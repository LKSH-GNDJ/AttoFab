// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A response surface and a yield estimate (PRD W-03).
//!
//! A Monte Carlo gives samples; a surface gives the *shape* of the response, so
//! a process engineer can see which direction is safe rather than only how often
//! a point passed. The fit is linear in the declared inputs, solved from the
//! normal equations, and the yield is the fraction of samples inside the spec.

#[derive(Debug, Clone, PartialEq)]
/// A fitted linear response over the study knobs.
pub struct ResponseSurface {
    /// One coefficient per input, in the order the inputs were given.
    pub coefficients: Vec<f64>,
    pub intercept: f64,
    /// Coefficient of determination, 0..1.
    pub r_squared: f64,
}

impl ResponseSurface {
    /// The fitted response at a point.
    pub fn predict(&self, inputs: &[f64]) -> f64 {
        self.intercept
            + self
                .coefficients
                .iter()
                .zip(inputs.iter())
                .map(|(c, x)| c * x)
                .sum::<f64>()
    }
}

/// Fit `y = b0 + sum(bi xi)` by least squares.
///
/// Returns `None` when there are fewer samples than inputs or the normal
/// equations are singular — a fit that cannot be trusted is not returned as a
/// fit with large error bars.
pub fn fit(inputs: &[Vec<f64>], outputs: &[f64]) -> Option<ResponseSurface> {
    let n = outputs.len();
    if n == 0 || inputs.len() != n {
        return None;
    }
    let p = inputs[0].len();
    if p == 0 || n < p + 1 || inputs.iter().any(|r| r.len() != p) {
        return None;
    }
    // normal equations for the design matrix [1, x...]
    let m = p + 1;
    let mut ata = vec![vec![0.0; m]; m];
    let mut aty = vec![0.0; m];
    for i in 0..n {
        let mut row = Vec::with_capacity(m);
        row.push(1.0);
        row.extend_from_slice(&inputs[i]);
        for a in 0..m {
            for b in 0..m {
                ata[a][b] += row[a] * row[b];
            }
            aty[a] += row[a] * outputs[i];
        }
    }
    let beta = solve(m, ata, aty)?;
    // r^2
    let mean = outputs.iter().sum::<f64>() / n as f64;
    let mut ss_res = 0.0;
    let mut ss_tot = 0.0;
    let surface = ResponseSurface { coefficients: beta[1..].to_vec(), intercept: beta[0], r_squared: 0.0 };
    for i in 0..n {
        let pred = surface.predict(&inputs[i]);
        ss_res += (outputs[i] - pred) * (outputs[i] - pred);
        ss_tot += (outputs[i] - mean) * (outputs[i] - mean);
    }
    let r_squared = if ss_tot > 0.0 { 1.0 - ss_res / ss_tot } else { 1.0 };
    Some(ResponseSurface { r_squared, ..surface })
}

/// Gaussian elimination with partial pivoting. `None` if singular.
fn solve(n: usize, mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-12 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);
        for r in col + 1..n {
            let f = a[r][col] / a[col][col];
            for c in col..n {
                a[r][c] -= f * a[col][c];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for c in r + 1..n {
            s -= a[r][c] * x[c];
        }
        x[r] = s / a[r][r];
    }
    Some(x)
}

/// Fraction of samples inside `[low, high]`.
pub fn yield_estimate(outputs: &[f64], low: f64, high: f64) -> Option<f64> {
    if outputs.is_empty() || !(low <= high) {
        return None;
    }
    let inside = outputs.iter().filter(|v| **v >= low && **v <= high).count();
    Some(inside as f64 / outputs.len() as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_linear_response_is_fitted_exactly() {
        // y = 3 + 2*x1 - 5*x2
        let inputs: Vec<Vec<f64>> = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.0, 1.0],
            vec![1.0, 1.0],
            vec![2.0, 3.0],
        ];
        let outputs: Vec<f64> = inputs.iter().map(|r| 3.0 + 2.0 * r[0] - 5.0 * r[1]).collect();
        let s = fit(&inputs, &outputs).expect("fits");
        assert!((s.intercept - 3.0).abs() < 1e-9, "intercept {:.6}", s.intercept);
        assert!((s.coefficients[0] - 2.0).abs() < 1e-9);
        assert!((s.coefficients[1] + 5.0).abs() < 1e-9);
        assert!((s.r_squared - 1.0).abs() < 1e-9, "r2 {:.6}", s.r_squared);
        assert!((s.predict(&[4.0, 2.0]) - (3.0 + 8.0 - 10.0)).abs() < 1e-9);
    }

    #[test]
    fn a_noisy_response_reports_a_less_than_perfect_fit() {
        let inputs: Vec<Vec<f64>> = (0..10).map(|i| vec![i as f64]).collect();
        let outputs: Vec<f64> = (0..10).map(|i| i as f64 + if i % 2 == 0 { 0.5 } else { -0.5 }).collect();
        let s = fit(&inputs, &outputs).expect("fits");
        assert!(s.r_squared > 0.94 && s.r_squared < 1.0, "r2 {:.4}", s.r_squared);
    }

    #[test]
    fn a_fit_that_cannot_be_trusted_is_refused() {
        assert!(fit(&[], &[]).is_none());
        assert!(fit(&[vec![1.0]], &[1.0]).is_none(), "one sample, one input");
        // two identical rows: the normal equations are singular
        let same = vec![vec![1.0], vec![1.0]];
        assert!(fit(&same, &[1.0, 2.0]).is_none(), "singular design");
        // ragged input
        assert!(fit(&[vec![1.0, 2.0], vec![1.0]], &[1.0, 2.0]).is_none());
    }

    #[test]
    fn the_yield_is_the_fraction_inside_the_spec() {
        let out = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        assert_eq!(yield_estimate(&out, 2.0, 4.0), Some(0.6));
        assert_eq!(yield_estimate(&out, -1.0, 10.0), Some(1.0));
        assert_eq!(yield_estimate(&out, 10.0, 20.0), Some(0.0));
        assert_eq!(yield_estimate(&[], 0.0, 1.0), None);
        assert_eq!(yield_estimate(&out, 5.0, 1.0), None);
    }
}
