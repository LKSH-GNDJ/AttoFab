// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Pass/fail under measurement uncertainty (build PRD M-04).
//!
//! A measurement near a limit cannot be judged pass or fail from the number
//! alone. This decides with a guard band: the value is accepted only when it is
//! inside the tolerance by more than `k` standard uncertainties, and the
//! probability of a false pass is reported alongside the verdict rather than
//! hidden behind the word "pass".

use crate::{Measurement, Status};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The three-state outcome of a guarded check.
pub enum Verdict {
    Pass,
    /// Inside the tolerance, but not by enough to call it: within the guard band.
    Indeterminate,
    Fail,
}

impl Verdict {
    /// The word a report uses.
    pub fn word(self) -> &'static str {
        match self {
            Verdict::Pass => "PASS",
            Verdict::Indeterminate => "INDETERMINATE",
            Verdict::Fail => "FAIL",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
/// A verdict with the risk that produced it.
pub struct Decision {
    pub verdict: Verdict,
    /// Distance to the tolerance limit in standard uncertainties. Negative when
    /// the value is outside the limit.
    pub margin_sigma: f64,
    /// Probability that the true value lies outside the tolerance.
    pub risk: f64,
    pub reason: String,
}

/// Decide one measurement, with `k` as the coverage factor (2 is ~95%).
pub fn decide(m: &Measurement, k: f64) -> Decision {
    let Some(v) = m.value else {
        return Decision {
            verdict: Verdict::Indeterminate,
            margin_sigma: f64::NEG_INFINITY,
            risk: 0.5,
            reason: match &m.reason {
                Some(r) => format!("not measured — {r}"),
                None => "not measured, and no reason was recorded".into(),
            },
        };
    };
    if m.status == Status::NotMeasured {
        return Decision {
            verdict: Verdict::Indeterminate,
            margin_sigma: f64::NEG_INFINITY,
            risk: 0.5,
            reason: "the status says not measured but a value was present".into(),
        };
    }
    let margin = m.tolerance - (v - m.nominal).abs();
    let sigma = if m.sigma > 0.0 { m.sigma } else { 0.0 };
    if sigma <= 0.0 {
        // No uncertainty declared: the number is all there is, and saying
        // "indeterminate" forever would be useless.
        return Decision {
            verdict: if margin >= 0.0 { Verdict::Pass } else { Verdict::Fail },
            margin_sigma: if margin >= 0.0 { f64::INFINITY } else { f64::NEG_INFINITY },
            risk: if margin >= 0.0 { 0.0 } else { 1.0 },
            reason: if margin >= 0.0 {
                "inside tolerance; no uncertainty declared".into()
            } else {
                "outside tolerance; no uncertainty declared".into()
            },
        };
    }
    let ms = margin / sigma;
    // one-sided: probability the true value is beyond the limit
    let risk = 0.5 * erfc(ms / std::f64::consts::SQRT_2);
    let (verdict, reason) = if margin < 0.0 {
        (Verdict::Fail, format!("{:.2} sigma outside the limit", -ms))
    } else if margin < k * sigma {
        (
            Verdict::Indeterminate,
            format!("inside the limit by {:.2} sigma, less than the {:.1} sigma guard band", ms, k),
        )
    } else {
        (Verdict::Pass, format!("inside the limit by {:.2} sigma", ms))
    };
    Decision { verdict, margin_sigma: ms, risk, reason }
}

/// Complementary error function, Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7).
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * z);
    let poly = t * (0.254829592
        + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let e = poly * (-z * z).exp();
    if x >= 0.0 { e } else { 2.0 - e }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(value: f64, nominal: f64, tol: f64, sigma: f64) -> Measurement {
        Measurement::against("CD", "nm", value, nominal, tol, sigma)
    }

    #[test]
    fn a_value_well_inside_passes_with_a_small_risk() {
        let d = decide(&m(100.0, 100.0, 10.0, 1.0), 2.0);
        assert_eq!(d.verdict, Verdict::Pass);
        assert!(d.risk < 1e-6, "risk {:.2e}", d.risk);
        assert!(d.margin_sigma > 9.0);
    }

    #[test]
    fn a_value_outside_the_tolerance_fails() {
        let d = decide(&m(112.0, 100.0, 10.0, 1.0), 2.0);
        assert_eq!(d.verdict, Verdict::Fail);
        assert!(d.risk > 0.9, "risk {:.3}", d.risk);
        assert!(d.margin_sigma < 0.0);
    }

    #[test]
    fn a_value_inside_but_within_the_guard_band_is_indeterminate() {
        let d = decide(&m(109.0, 100.0, 10.0, 1.0), 2.0);
        assert_eq!(d.verdict, Verdict::Indeterminate, "{}", d.reason);
        assert!(d.risk > 0.10 && d.risk < 0.17, "risk {:.3}", d.risk);
        assert!(d.reason.contains("guard band"));
    }

    #[test]
    fn the_indeterminate_band_widens_with_the_coverage_factor() {
        let mm = m(108.0, 100.0, 10.0, 1.0);
        assert_eq!(decide(&mm, 1.0).verdict, Verdict::Pass, "2 sigma in, k=1");
        assert_eq!(decide(&mm, 3.0).verdict, Verdict::Indeterminate, "2 sigma in, k=3");
    }

    #[test]
    fn a_missing_measurement_is_indeterminate_and_says_why() {
        let d = decide(&Measurement::missing("CD", "nm", "no resist in the window"), 2.0);
        assert_eq!(d.verdict, Verdict::Indeterminate);
        assert!(d.reason.contains("no resist"), "{}", d.reason);
    }

    #[test]
    fn no_declared_uncertainty_still_gives_a_verdict() {
        assert_eq!(decide(&m(105.0, 100.0, 10.0, 0.0), 2.0).verdict, Verdict::Pass);
        assert_eq!(decide(&m(115.0, 100.0, 10.0, 0.0), 2.0).verdict, Verdict::Fail);
    }

    #[test]
    fn the_error_function_matches_published_values() {
        assert!((erfc(0.0) - 1.0).abs() < 1e-7);
        assert!((erfc(1.0) - 0.157299).abs() < 1e-6);
        assert!((erfc(2.0) - 0.0046777).abs() < 1e-6);
        assert!((erfc(-1.0) - 1.842701).abs() < 1e-6);
    }

    #[test]
    fn the_risk_is_one_sided_and_monotone_in_the_margin() {
        let a = decide(&m(104.0, 100.0, 10.0, 1.0), 2.0).risk;
        let b = decide(&m(108.0, 100.0, 10.0, 1.0), 2.0).risk;
        assert!(b > a, "risk must rise as the margin shrinks: {a:.3e} -> {b:.3e}");
        assert!(b < 0.5);
    }
}
