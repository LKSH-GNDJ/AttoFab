// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Sidewall passivation and polymer (E-09).
//!
//! A chemical etch with no ion guidance is isotropic: the lateral rate equals the vertical rate.
//! Anisotropy in a real plasma etch comes mostly from a polymer film that deposits on the
//! sidewall, where the ions arrive at a grazing angle and do not remove it, while the ions at
//! the feature bottom strike at normal incidence and do. The film shields the wall from the
//! chemistry. This module is that balance for the sidewall alone:
//!
//! `d(theta)/dt = k_dep * f * (1 - theta) - k_rem * theta`
//!
//! where `theta` is the fraction of the wall covered, `f` is the passivant fraction in the gas,
//! `k_dep` is the deposition rate constant at `f = 1` and `k_rem` is the removal constant. The
//! lateral etch rate is `L * (1 - theta)`. With `k = k_dep * f + k_rem` and
//! `theta_ss = k_dep * f / k` the solution is `theta(t) = theta_ss * (1 - exp(-k t))`, and the
//! undercut (the integral of the lateral rate) is
//! `L * (t - theta_ss * (t - (1 - exp(-k t)) / k))`.
//!
//! Coverage takes time to build, so a short etch undercuts more than a long one in proportion:
//! the profile steepens with depth. The bottom is taken as unpassivated, so the vertical rate is
//! the input vertical rate.
//!
//! **Maturity: L1.** The mechanism is the published qualitative one (sidewall passivation, Coburn
//! and Winters style ion-enhanced etching), but `k_dep` and `k_rem` are registered as *assumed*:
//! the tree holds no measured sidewall profile to fit them to.

use atto_units::{Fraction, Minutes, Nm, NmPerMin};

/// Deposition rate constant of the polymer at a passivant fraction of 1, per minute. Assumed.
pub const K_DEP_PER_MIN: f64 = 5.0;
/// Ion-assisted removal constant of the polymer from the wall, per minute. Assumed.
pub const K_REM_PER_MIN: f64 = 0.2;

/// The two rate constants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PassivationParams {
    /// Deposition constant at a passivant fraction of 1, per minute.
    pub k_dep_per_min: f64,
    /// Removal constant, per minute.
    pub k_rem_per_min: f64,
}

impl PassivationParams {
    /// The bundled, assumed constants.
    pub fn bundled() -> Self {
        PassivationParams { k_dep_per_min: K_DEP_PER_MIN, k_rem_per_min: K_REM_PER_MIN }
    }
}

/// Why a passivation could not be evaluated.
#[derive(Debug, Clone, PartialEq)]
pub struct PassivationError {
    /// What was wrong.
    pub what: String,
}

impl std::fmt::Display for PassivationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.what)
    }
}

impl std::error::Error for PassivationError {}

/// What the sidewall looks like after an etch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SidewallResult {
    /// Fraction of the wall covered by polymer at the end of the etch.
    pub coverage: f64,
    /// Lateral etch on each side.
    pub undercut: Nm,
    /// Vertical etch depth.
    pub depth: Nm,
    /// Sidewall angle from the horizontal, degrees: 90 is vertical, 45 is a fully isotropic etch.
    pub angle_deg: f64,
}

impl SidewallResult {
    /// One line for a report.
    pub fn label(&self) -> String {
        format!("wall {:.1} deg, undercut {:.1} nm at depth {:.0} nm, coverage {:.2}", self.angle_deg, self.undercut.0, self.depth.0, self.coverage)
    }
}

/// Etch for `time` with a chemistry whose bare lateral and vertical rates are given, and a gas
/// carrying `passivant` as the fraction of polymer-forming species.
pub fn sidewall(
    p: &PassivationParams,
    passivant: Fraction,
    lateral: NmPerMin,
    vertical: NmPerMin,
    time: Minutes,
) -> Result<SidewallResult, PassivationError> {
    let f = passivant.0;
    if !(f.is_finite() && (0.0..=1.0).contains(&f)) {
        return Err(PassivationError { what: "the passivant fraction must lie in [0, 1]".into() });
    }
    for (v, what) in [(lateral.0, "lateral rate"), (vertical.0, "vertical rate"), (time.0, "time")] {
        if !(v.is_finite() && v >= 0.0) {
            return Err(PassivationError { what: format!("the {what} must be finite and not negative") });
        }
    }
    if !(p.k_dep_per_min.is_finite() && p.k_dep_per_min >= 0.0 && p.k_rem_per_min.is_finite() && p.k_rem_per_min >= 0.0) {
        return Err(PassivationError { what: "the rate constants must be finite and not negative".into() });
    }
    let t = time.0;
    let k = p.k_dep_per_min * f + p.k_rem_per_min;
    let (theta_ss, decay) = if k > 0.0 { (p.k_dep_per_min * f / k, (1.0 - (-k * t).exp()) / k) } else { (0.0, t) };
    let coverage = if k > 0.0 { theta_ss * (1.0 - (-k * t).exp()) } else { 0.0 };
    let undercut = lateral.0 * (t - theta_ss * (t - decay));
    let depth = vertical.0 * t;
    let angle_deg = depth.atan2(undercut).to_degrees();
    Ok(SidewallResult { coverage, undercut: Nm(undercut), depth: Nm(depth), angle_deg })
}

#[cfg(test)]
mod tests {
    use super::*;

    const L: NmPerMin = NmPerMin(300.0);
    const V: NmPerMin = NmPerMin(300.0);

    #[test]
    fn without_passivation_an_isotropic_chemistry_etches_at_forty_five_degrees() {
        let r = sidewall(&PassivationParams::bundled(), Fraction(0.0), L, V, Minutes(5.0)).unwrap();
        assert!((r.undercut.0 - 1500.0).abs() < 1e-9, "undercut {}", r.undercut.0);
        assert!((r.angle_deg - 45.0).abs() < 0.5, "angle {}", r.angle_deg);
    }

    #[test]
    fn passivation_makes_an_isotropic_chemistry_vertical() {
        let r = sidewall(&PassivationParams::bundled(), Fraction(1.0), L, V, Minutes(5.0)).unwrap();
        assert!(r.angle_deg >= 85.0, "angle {}", r.angle_deg);
        assert!(r.coverage > 0.9);
    }

    #[test]
    fn more_passivant_gives_a_steeper_wall() {
        let mut prev = 0.0;
        for f in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let a = sidewall(&PassivationParams::bundled(), Fraction(f), L, V, Minutes(5.0)).unwrap().angle_deg;
            assert!(a > prev, "angle {a} at f {f} did not exceed {prev}");
            prev = a;
        }
    }

    #[test]
    fn the_closed_form_matches_a_numerical_integration() {
        // An independent check: step theta forward and integrate the lateral rate, no formula shared.
        let p = PassivationParams::bundled();
        let f = 0.6;
        let t_end = 3.0;
        let n = 200_000;
        let dt = t_end / n as f64;
        let (mut theta, mut undercut) = (0.0f64, 0.0f64);
        for _ in 0..n {
            undercut += L.0 * (1.0 - theta) * dt;
            theta += (p.k_dep_per_min * f * (1.0 - theta) - p.k_rem_per_min * theta) * dt;
        }
        let r = sidewall(&p, Fraction(f), L, V, Minutes(t_end)).unwrap();
        assert!((r.undercut.0 - undercut).abs() / undercut < 0.01, "{} vs {undercut}", r.undercut.0);
        assert!((r.coverage - theta).abs() < 0.01);
    }

    #[test]
    fn bad_inputs_are_refused() {
        let p = PassivationParams::bundled();
        assert!(sidewall(&p, Fraction(1.5), L, V, Minutes(1.0)).is_err());
        assert!(sidewall(&p, Fraction(0.5), NmPerMin(-1.0), V, Minutes(1.0)).is_err());
        assert!(sidewall(&p, Fraction(0.5), L, V, Minutes(f64::NAN)).is_err());
    }
}
