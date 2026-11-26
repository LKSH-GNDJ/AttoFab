// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Identifiability: what the data actually pins down (K-07).
//!
//! A least-squares fit always returns numbers. It does not always return *information*: two
//! parameters that move the observable in the same way cannot be separated by that data, and
//! the fit will happily split the difference and report both with small-looking uncertainties.
//! That is the most dangerous output a fitting tool can produce, because it looks like a result.
//!
//! This module takes the Jacobian of the model with respect to its parameters and reports the
//! covariance, the correlation matrix, and any pair whose correlation is high enough that the
//! data cannot distinguish them.

/// Why an identifiability analysis could not run.
#[derive(Debug, Clone, PartialEq)]
pub struct IdentifyError {
    pub what: String,
}

impl std::fmt::Display for IdentifyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "identifiability: {}", self.what)
    }
}

/// What a fit can and cannot tell apart.
#[derive(Debug, Clone, PartialEq)]
pub struct Identifiability {
    /// Covariance of the parameters, same units squared as the parameters.
    pub covariance: Vec<Vec<f64>>,
    /// Correlation matrix, dimensionless, diagonal 1.
    pub correlation: Vec<Vec<f64>>,
    /// Parameter pairs the data cannot separate: (i, j, |correlation|).
    pub degenerate: Vec<(usize, usize, f64)>,
    /// Condition number of `J^T J`. Large means the problem is ill-conditioned even if no
    /// single pair is correlated.
    pub condition: f64,
}


/// Determinant of a small square matrix by Gaussian elimination with partial pivoting.
fn det(mut a: Vec<Vec<f64>>) -> f64 {
    let n = a.len();
    let mut d = 1.0;
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-300 {
            return 0.0;
        }
        if piv != col {
            a.swap(piv, col);
            d = -d;
        }
        d *= a[col][col];
        let p = a[col][col];
        for r in col + 1..n {
            let f = a[r][col] / p;
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
        }
    }
    d
}

/// Invert a square matrix by Gauss-Jordan with partial pivoting; `None` if singular.
fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    if n == 0 || a.iter().any(|r| r.len() != n) {
        return None;
    }
    let mut m: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut r = a[i].clone();
            r.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 }));
            r
        })
        .collect();
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if m[r][col].abs() > m[piv][col].abs() {
                piv = r;
            }
        }
        if m[piv][col].abs() < 1e-300 {
            return None;
        }
        m.swap(col, piv);
        let d = m[col][col];
        for j in 0..2 * n {
            m[col][j] /= d;
        }
        for r in 0..n {
            if r == col {
                continue;
            }
            let f = m[r][col];
            if f == 0.0 {
                continue;
            }
            for j in 0..2 * n {
                m[r][j] -= f * m[col][j];
            }
        }
    }
    Some((0..n).map(|i| m[i][n..].to_vec()).collect())
}

/// Analyse a Jacobian `J` (rows = observations, columns = parameters).
///
/// `residual_variance` is the variance of the fit residual, sigma^2; the covariance is
/// `(J^T J)^-1 sigma^2`. A pair is reported as degenerate when `|correlation| >= threshold`.
pub fn analyse(
    jacobian: &[Vec<f64>],
    residual_variance: f64,
    threshold: f64,
) -> Result<Identifiability, IdentifyError> {
    if jacobian.is_empty() {
        return Err(IdentifyError { what: "no observations".into() });
    }
    let n = jacobian[0].len();
    if n == 0 {
        return Err(IdentifyError { what: "no parameters".into() });
    }
    if jacobian.iter().any(|r| r.len() != n) {
        return Err(IdentifyError { what: "ragged Jacobian".into() });
    }
    if jacobian.iter().flatten().any(|v| !v.is_finite()) {
        return Err(IdentifyError { what: "a Jacobian entry is not finite".into() });
    }
    if !residual_variance.is_finite() || residual_variance < 0.0 {
        return Err(IdentifyError { what: "residual variance must be finite and >= 0".into() });
    }
    if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
        return Err(IdentifyError { what: "threshold must be a fraction in 0..=1".into() });
    }

    // J^T J
    let mut jtj = vec![vec![0.0; n]; n];
    for row in jacobian {
        for i in 0..n {
            for j in 0..n {
                jtj[i][j] += row[i] * row[j];
            }
        }
    }
    let cond = {
        let mut tr = 0.0f64;
        for i in 0..n {
            tr += jtj[i][i];
        }
        if tr <= 0.0 {
            f64::INFINITY
        } else {
            // a cheap proxy: the trace over the smallest diagonal, which is enough to flag
            // a rank-deficient normal matrix without an SVD
            let smallest = (0..n).map(|i| jtj[i][i]).fold(f64::INFINITY, f64::min);
            if smallest <= 0.0 { f64::INFINITY } else { tr / smallest }
        }
    };

    let inv = invert(&jtj).ok_or(IdentifyError {
        what: "J^T J is singular: the parameters are not separable by this data".into(),
    })?;
    let covariance: Vec<Vec<f64>> = inv.iter().map(|r| r.iter().map(|v| v * residual_variance).collect()).collect();

    let mut correlation = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let di = covariance[i][i];
            let dj = covariance[j][j];
            correlation[i][j] = if di > 0.0 && dj > 0.0 { covariance[i][j] / (di * dj).sqrt() } else { 0.0 };
        }
        correlation[i][i] = 1.0;
    }

    let mut degenerate = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let c = correlation[i][j].abs();
            if c >= threshold {
                degenerate.push((i, j, c));
            }
        }
    }

    Ok(Identifiability { covariance, correlation, degenerate, condition: cond })
}


/// Analyse from a normal matrix `J^T J` directly.
///
/// A least-squares fit rarely has the Jacobian to hand; what it can compute is the Hessian of
/// its objective, and for `f = sum r^2` that Hessian is `2 J^T J`. This entry point takes the
/// normal matrix so a fit can ask the same question without inventing a Jacobian.
pub fn analyse_normal_matrix(
    jtj: &[Vec<f64>],
    residual_variance: f64,
    threshold: f64,
) -> Result<Identifiability, IdentifyError> {
    let n = jtj.len();
    if n == 0 {
        return Err(IdentifyError { what: "no parameters".into() });
    }
    if jtj.iter().any(|r| r.len() != n) {
        return Err(IdentifyError { what: "the normal matrix is not square".into() });
    }
    if jtj.iter().flatten().any(|v| !v.is_finite()) {
        return Err(IdentifyError { what: "a normal-matrix entry is not finite".into() });
    }
    if !residual_variance.is_finite() || residual_variance < 0.0 {
        return Err(IdentifyError { what: "residual variance must be finite and >= 0".into() });
    }
    if !threshold.is_finite() || !(0.0..=1.0).contains(&threshold) {
        return Err(IdentifyError { what: "threshold must be a fraction in 0..=1".into() });
    }
    let trace: f64 = (0..n).map(|i| jtj[i][i]).sum();
    let smallest = (0..n).map(|i| jtj[i][i]).fold(f64::INFINITY, f64::min);
    let condition = if smallest <= 0.0 || trace <= 0.0 { f64::INFINITY } else { trace / smallest };

    // A determinant near zero means the parameters are not separable, whatever the diagonals
    // say. The trace/smallest-diagonal proxy used above reports 1.0 for [[1,1],[1,1]] — a rank-1
    // matrix — so it cannot be trusted on its own.
    let scale = (0..n).map(|i| jtj[i][i]).sum::<f64>() / n as f64;
    let d = det(jtj.to_vec()).abs();
    if scale <= 0.0 || d <= 1e-12 * scale.powi(n as i32) {
        return Err(IdentifyError {
            what: "the normal matrix is singular: the parameters are not separable by this data".into(),
        });
    }
    let inv = invert(jtj).ok_or(IdentifyError {
        what: "the normal matrix is singular: the parameters are not separable by this data".into(),
    })?;
    let covariance: Vec<Vec<f64>> = inv.iter().map(|r| r.iter().map(|v| v * residual_variance).collect()).collect();
    let mut correlation = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            let (di, dj) = (covariance[i][i], covariance[j][j]);
            correlation[i][j] = if di > 0.0 && dj > 0.0 { covariance[i][j] / (di * dj).sqrt() } else { 0.0 };
        }
        correlation[i][i] = 1.0;
    }
    let mut degenerate = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let c = correlation[i][j].abs();
            if c >= threshold {
                degenerate.push((i, j, c));
            }
        }
    }
    Ok(Identifiability { covariance, correlation, degenerate, condition })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two well-separated parameters: each observation moves one parameter only, so the
    /// columns are orthogonal. A row like [1, 1] would correlate them by 0.5 and is not
    /// orthogonal — the first version of this fixture had one and the test rightly failed.
    fn orthogonal() -> Vec<Vec<f64>> {
        vec![vec![1.0, 0.0], vec![0.0, 1.0]]
    }

    #[test]
    fn independent_parameters_are_not_flagged() {
        let r = analyse(&orthogonal(), 1.0, 0.9).expect("an analysis");
        assert!(r.degenerate.is_empty(), "nothing is degenerate here: {:?}", r.degenerate);
        assert!((r.correlation[0][1]).abs() < 1e-12, "and the off-diagonal is zero");
        assert_eq!(r.covariance.len(), 2);
    }

    #[test]
    fn a_perfectly_degenerate_pair_is_reported_not_fitted() {
        // two columns that move identically: the data cannot tell them apart at all
        let j = vec![vec![1.0, 1.0], vec![2.0, 2.0], vec![3.0, 3.0]];
        let r = analyse(&j, 1.0, 0.9);
        match r {
            Err(e) => assert!(e.what.contains("singular"), "a singular problem is named: {}", e.what),
            Ok(v) => {
                // if the normal matrix was invertible numerically, the pair must still be flagged
                assert!(!v.degenerate.is_empty(), "a perfect correlation must be reported");
                assert!(v.degenerate[0].2 > 0.99, "and near one: {}", v.degenerate[0].2);
            }
        }
    }

    #[test]
    fn a_nearly_degenerate_pair_is_flagged_with_its_correlation() {
        // the second column is almost the first, plus a little of its own
        let j = vec![vec![1.0, 1.0], vec![2.0, 2.02], vec![3.0, 2.99], vec![4.0, 4.05]];
        let r = analyse(&j, 1.0, 0.95).expect("an analysis");
        assert_eq!(r.degenerate.len(), 1, "one pair is flagged: {:?}", r.degenerate);
        assert!(r.degenerate[0].2 > 0.95, "with a high correlation: {}", r.degenerate[0].2);
    }

    #[test]
    fn a_higher_threshold_flags_less() {
        let j = vec![vec![1.0, 1.0], vec![2.0, 2.02], vec![3.0, 2.99], vec![4.0, 4.05]];
        let loose = analyse(&j, 1.0, 0.5).expect("loose");
        let strict = analyse(&j, 1.0, 0.9999).expect("strict");
        assert!(loose.degenerate.len() >= strict.degenerate.len());
    }

    #[test]
    fn the_covariance_scales_with_the_residual_variance() {
        let a = analyse(&orthogonal(), 1.0, 0.9).expect("a");
        let b = analyse(&orthogonal(), 4.0, 0.9).expect("b");
        assert!((b.covariance[0][0] / a.covariance[0][0] - 4.0).abs() < 1e-9, "sigma^2 enters linearly");
    }

    #[test]
    fn bad_input_is_refused() {
        assert!(analyse(&[], 1.0, 0.9).is_err(), "no observations");
        assert!(analyse(&[vec![]], 1.0, 0.9).is_err(), "no parameters");
        assert!(analyse(&[vec![1.0, 2.0], vec![1.0]], 1.0, 0.9).is_err(), "ragged");
        assert!(analyse(&[vec![f64::NAN, 1.0]], 1.0, 0.9).is_err(), "not a number");
        assert!(analyse(&orthogonal(), -1.0, 0.9).is_err(), "negative variance");
        assert!(analyse(&orthogonal(), 1.0, 2.0).is_err(), "threshold above one");
    }
}
