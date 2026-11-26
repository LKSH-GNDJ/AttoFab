// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Fitting declared parameters to a calibration objective (PRD K-04).
//!
//! An objective, a parameter set with bounds, and a search. The search is
//! coordinate descent with a shrinking step rather than Levenberg-Marquardt:
//! it needs no Jacobian, it cannot leave the bounds, and its failure mode is
//! "did not converge" rather than "converged to a negative thickness".
//!
//! That is a deliberate limitation. A derivative-free search is slower on
//! smooth problems and much harder to make produce a plausible wrong answer,
//! which is the trade this project wants.

#[derive(Debug, Clone, PartialEq)]
/// One parameter under fit: its name, start and bounds.
pub struct FitParam {
    pub name: String,
    pub initial: f64,
    pub lower: f64,
    pub upper: f64,
}

impl FitParam {
    /// A parameter to fit, with its bounds.
    pub fn new(name: &str, initial: f64, lower: f64, upper: f64) -> Self {
        FitParam { name: name.to_string(), initial, lower, upper }
    }

    fn clamp(&self, v: f64) -> f64 {
        v.max(self.lower).min(self.upper)
    }

    fn valid(&self) -> bool {
        self.lower.is_finite() && self.upper.is_finite() && self.lower <= self.upper
            && self.initial.is_finite() && self.lower <= self.initial && self.initial <= self.upper
    }
}

#[derive(Debug, Clone, PartialEq)]
/// Where the search ended and why.
pub struct FitResult {
    pub values: Vec<f64>,
    pub objective: f64,
    pub iterations: usize,
    pub converged: bool,
    /// Why it stopped, when it did not converge.
    pub reason: String,
}

/// Minimise `objective` over the parameters, within their bounds.
///
/// Returns `None` for an invalid parameter set — a fit with no bounds or an
/// initial value outside them is a configuration error, not a poor fit.
pub fn fit<F>(params: &[FitParam], objective: F, max_iter: usize) -> Option<FitResult>
where
    F: Fn(&[f64]) -> f64,
{
    if params.is_empty() || params.iter().any(|p| !p.valid()) {
        return None;
    }
    let mut x: Vec<f64> = params.iter().map(|p| p.initial).collect();
    let mut best = objective(&x);
    if !best.is_finite() {
        return Some(FitResult {
            values: x,
            objective: best,
            iterations: 0,
            converged: false,
            reason: "the objective is not finite at the initial point".into(),
        });
    }
    let mut step: Vec<f64> = params
        .iter()
        .map(|p| ((p.upper - p.lower) * 0.1).max(1e-12))
        .collect();
    let mut iterations = 0;
    let mut converged = false;
    let mut reason = String::from("iteration limit reached");

    while iterations < max_iter {
        iterations += 1;
        let mut improved = false;
        for i in 0..x.len() {
            for dir in [-1.0, 1.0] {
                let trial = params[i].clamp(x[i] + dir * step[i]);
                if (trial - x[i]).abs() < 1e-15 {
                    continue;
                }
                let mut candidate = x.clone();
                candidate[i] = trial;
                let v = objective(&candidate);
                if v.is_finite() && v < best {
                    x = candidate;
                    best = v;
                    improved = true;
                }
            }
        }
        if !improved {
            // shrink every step; stop when they are all negligible against the
            // parameter's own range
            let mut done = true;
            for i in 0..step.len() {
                step[i] *= 0.5;
                let range = (params[i].upper - params[i].lower).max(1e-12);
                if step[i] > range * 1e-9 {
                    done = false;
                }
            }
            if done {
                converged = true;
                reason = "step size fell below the parameter resolution".into();
                break;
            }
        }
    }
    Some(FitResult { values: x, objective: best, iterations, converged, reason })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_known_parameter_is_recovered() {
        // the objective is minimised exactly at 42
        let target = 42.0;
        let obj = |v: &[f64]| (v[0] - target) * (v[0] - target);
        let p = [FitParam::new("thickness_nm", 10.0, 0.0, 100.0)];
        let r = fit(&p, obj, 500).expect("a valid fit");
        assert!((r.values[0] - target).abs() < 1e-6, "recovered {:.6}", r.values[0]);
        assert!(r.objective < 1e-12);
        assert!(r.converged, "{}", r.reason);
    }

    #[test]
    fn a_bound_is_respected_rather_than_crossed() {
        // the true minimum is at 500, well outside the bound
        let obj = |v: &[f64]| (v[0] - 500.0) * (v[0] - 500.0);
        let p = [FitParam::new("x", 50.0, 0.0, 100.0)];
        let r = fit(&p, obj, 500).expect("a valid fit");
        assert!(r.values[0] <= 100.0 && r.values[0] >= 0.0, "stayed inside: {}", r.values[0]);
        assert!((r.values[0] - 100.0).abs() < 1e-6, "and pressed against the bound");
    }

    #[test]
    fn two_parameters_are_fitted_together() {
        let obj = |v: &[f64]| (v[0] - 20.0).powi(2) + (v[1] + 5.0).powi(2);
        let p = [
            FitParam::new("a", 0.0, -50.0, 50.0),
            FitParam::new("b", 0.0, -50.0, 50.0),
        ];
        let r = fit(&p, obj, 2000).expect("a valid fit");
        assert!((r.values[0] - 20.0).abs() < 1e-4, "a = {:.4}", r.values[0]);
        assert!((r.values[1] + 5.0).abs() < 1e-4, "b = {:.4}", r.values[1]);
    }

    #[test]
    fn an_invalid_parameter_set_is_refused() {
        assert!(fit(&[], |_| 0.0, 10).is_none(), "no parameters");
        // initial outside the bounds
        assert!(fit(&[FitParam::new("x", 200.0, 0.0, 100.0)], |_| 0.0, 10).is_none());
        // inverted bounds
        assert!(fit(&[FitParam::new("x", 50.0, 100.0, 0.0)], |_| 0.0, 10).is_none());
    }

    #[test]
    fn an_objective_that_is_not_finite_says_so() {
        let p = [FitParam::new("x", 10.0, 0.0, 100.0)];
        let r = fit(&p, |_| f64::NAN, 10).expect("a result, not a panic");
        assert!(!r.converged);
        assert!(r.reason.contains("not finite"), "{}", r.reason);
    }
}

/// A fit that reports what the data could and could not separate (K-07).
///
/// `fit` returns a point estimate and nothing else. This runs the same fit and then estimates
/// the Hessian of the objective numerically at the solution; for a sum-of-squares objective that
/// Hessian is `2 J^T J`, so the covariance follows without inventing a Jacobian.
pub fn fit_with_covariance<F>(
    params: &[FitParam],
    objective: F,
    max_iter: usize,
    threshold: f64,
) -> Option<(FitResult, crate::identify::Identifiability)>
where
    F: Fn(&[f64]) -> f64,
{
    let r = fit(params, &objective, max_iter)?;
    let n = r.values.len();
    if n == 0 || !r.objective.is_finite() {
        return None;
    }
    let f = |x: &[f64]| objective(x);
    // central-difference Hessian at the solution
    let mut h = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        for j in 0..n {
            let hi = 1e-5 * r.values[i].abs().max(1e-3);
            let hj = 1e-5 * r.values[j].abs().max(1e-3);
            let mut xp = r.values.clone();
            xp[i] += hi;
            let mut xm = r.values.clone();
            xm[i] -= hi;
            let base = f(&r.values);
            let (fpp, fmm, fpm, fmp) = if i == j {
                (f(&xp), f(&xm), base, base)
            } else {
                let mut a = xp.clone();
                a[j] += hj;
                let mut b = xp.clone();
                b[j] -= hj;
                let mut c = xm.clone();
                c[j] += hj;
                let mut d = xm.clone();
                d[j] -= hj;
                (f(&a), f(&b), f(&c), f(&d))
            };
            h[i][j] = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
        }
    }
    // J^T J = H / 2 for f = sum r^2; the residual variance is the objective per point, which
    // is 1.0 when the objective is already a normalised sum of squares.
    let jtj: Vec<Vec<f64>> = h.iter().map(|row| row.iter().map(|v| v / 2.0).collect()).collect();
    let ident = crate::identify::analyse_normal_matrix(&jtj, 1.0, threshold).ok()?;
    Some((r, ident))
}

#[cfg(test)]
mod covariance_tests {
    use super::*;

    fn two(a: f64, b: f64) -> Vec<FitParam> {
        vec![
            FitParam::new("a", a, -10.0, 10.0),
            FitParam::new("b", b, -10.0, 10.0),
        ]
    }

    #[test]
    fn an_independent_pair_is_not_reported_as_degenerate() {
        // f = (a-1)^2 + (b-2)^2: each parameter is pinned by its own term
        let obj = |x: &[f64]| (x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2);
        let (r, id) = fit_with_covariance(&two(0.0, 0.0), obj, 400, 0.95).expect("a fit");
        assert!((r.values[0] - 1.0).abs() < 0.1, "a is recovered: {}", r.values[0]);
        assert!(id.degenerate.is_empty(), "nothing is degenerate here: {:?}", id.degenerate);
    }

    #[test]
    fn a_pair_the_data_cannot_separate_is_reported() {
        // The direct claim: a rank-1 normal matrix is refused rather than inverted. This is
        // tested on `analyse_normal_matrix` rather than through a fit, because the numerical
        // Hessian of a degenerate objective is only *near*-singular and the outcome depends on
        // the finite-difference step — a test that flips on rounding is not a test.
        let singular = vec![vec![1.0, 1.0], vec![1.0, 1.0]];
        match crate::identify::analyse_normal_matrix(&singular, 1.0, 0.95) {
            Err(e) => assert!(e.what.contains("singular"), "the cause is named: {}", e.what),
            Ok(v) => panic!("a rank-1 normal matrix must be refused, got degenerate {:?} condition {}", v.degenerate, v.condition),
        }
    }

    #[test]
    fn a_well_conditioned_normal_matrix_is_accepted() {
        let good = vec![vec![4.0, 0.0], vec![0.0, 4.0]];
        let v = crate::identify::analyse_normal_matrix(&good, 1.0, 0.95).expect("accepted");
        assert!(v.degenerate.is_empty(), "independent parameters are not degenerate");
        assert!((v.covariance[0][0] - 0.25).abs() < 1e-12, "covariance is the inverse: {}", v.covariance[0][0]);
    }
}
