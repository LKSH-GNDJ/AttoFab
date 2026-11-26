// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Point-defect pair diffusion and transient enhanced diffusion, D-13.
//!
//! An implant leaves the silicon supersaturated with self-interstitials. Boron diffuses by
//! swapping with those interstitials, so while the supersaturation lasts boron moves far faster
//! than its equilibrium diffusivity — and then it stops. That is **transient** enhanced diffusion,
//! and it is why an implant anneal broadens a boron spike in the first seconds and barely moves it
//! afterwards.
//!
//! The interstitial population recombines with vacancies, which is second order in the
//! supersaturation, so it decays as `1/(1 + t/tau)` rather than exponentially. A model with an
//! exponential decay gets the shape of TED wrong: too fast at the start, too slow at the end.
use atto_units::{Minutes};

/// A point-defect population relaxing by recombination.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairConfig {
    /// Initial interstitial supersaturation, `I/I*`. 1.0 is equilibrium, above 1.0 is an implant.
    pub initial_supersaturation: f64,
    /// Recombination rate constant, per minute.
    pub k_recombine: f64,
}

/// Why a pair-diffusion calculation could not run.
#[derive(Debug, Clone, PartialEq)]
pub struct PairError {
    pub what: String,
}

impl std::fmt::Display for PairError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "pair diffusion: {}", self.what)
    }
}

impl PairConfig {
    /// A typical implant supersaturation.
    ///
    /// **Assumed.** `10x` supersaturation and a recombination constant that makes the transient
    /// last a few minutes are the observed order; neither is fitted here. Registered as
    /// `pair.initial.supersaturation` and `pair.k.recombine`.
    pub fn implant() -> Self {
        PairConfig { initial_supersaturation: 10.0, k_recombine: 0.5 }
    }

    fn check(&self) -> Result<(), PairError> {
        if !self.initial_supersaturation.is_finite() || self.initial_supersaturation < 1.0 {
            return Err(PairError { what: "the initial supersaturation must be finite and >= 1".into() });
        }
        if !self.k_recombine.is_finite() || self.k_recombine <= 0.0 {
            return Err(PairError { what: "the recombination constant must be finite and > 0".into() });
        }
        Ok(())
    }
}

/// The supersaturation `I/I*` at a time, for second-order recombination.
///
/// `dS/dt = -k S^2` gives `S(t) = S0 / (1 + S0 k t)` — a hyperbolic decay, not an exponential.
pub fn supersaturation(t_min: Minutes, c: &PairConfig) -> Result<f64, PairError> {
    let t_min = t_min.0;
    c.check()?;
    if !t_min.is_finite() || t_min < 0.0 {
        return Err(PairError { what: "time must be finite and >= 0 min".into() });
    }
    let s0 = c.initial_supersaturation;
    // Recombination can only bring the population back to equilibrium, never below it: the
    // uncapped 1/(1+t/tau) keeps falling past 1 and would make the enhancement negative.
    Ok((s0 / (1.0 + s0 * c.k_recombine * t_min)).max(1.0))
}

/// The diffusion enhancement factor at a time.
///
/// The dopant diffusivity scales with the interstitial supersaturation, so this is the
/// supersaturation itself. It returns to 1.0 as the defects recombine, which is what makes the
/// enhancement transient.
pub fn enhancement(t_min: Minutes, c: &PairConfig) -> Result<f64, PairError> {
    let t_min = t_min.0;
    supersaturation(atto_units::Minutes(t_min), c)
}

/// The integrated enhancement — how much extra diffusion the transient delivers in total.
///
/// `integral of (S(t) - 1) dt`, which is the quantity that decides whether a boron spike moves at
/// all. It grows as `ln(1 + S0 k t)/k` minus the equilibrium part.
pub fn integrated_enhancement(t_min: Minutes, c: &PairConfig) -> Result<f64, PairError> {
    let t_min = t_min.0;
    c.check()?;
    if !t_min.is_finite() || t_min < 0.0 {
        return Err(PairError { what: "time must be finite and >= 0 min".into() });
    }
    let s0 = c.initial_supersaturation;
    let k = c.k_recombine;
    // The transient ends when the supersaturation has fallen back to 1.0, at
    // t_end = (S0 - 1)/(S0 k). Integrating S - 1 up to that point gives
    // (ln S0 - (S0 - 1)/S0)/k, and nothing is added afterwards.
    let t_end = (s0 - 1.0) / (s0 * k);
    let t = t_min.min(t_end);
    let total = (1.0 + s0 * k * t).ln() / k;
    Ok((total - t).max(0.0))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cfg() -> PairConfig {
        PairConfig::implant()
    }

    #[test]
    fn the_supersaturation_decays_to_equilibrium() {
        let c = cfg();
        let early = supersaturation(atto_units::Minutes(0.0), &c).expect("early");
        let late = supersaturation(atto_units::Minutes(1000.0), &c).expect("late");
        assert!((early - 10.0).abs() < 1e-12, "it starts at the implant value: {early}");
        assert!(late < 1.01, "and relaxes to equilibrium: {late}");
    }

    #[test]
    fn the_decay_is_hyperbolic_not_exponential() {
        // The signature of 1/(1+bt) is that 1/S is LINEAR in time; for an exponential it is ln S
        // that is linear. The first version of this test used a ratio identity that holds for
        // neither, and failed for that reason rather than because the model was wrong.
        let c = cfg();
        let s = |t: f64| supersaturation(atto_units::Minutes(t), &c).expect("s");
        let (t0, t1, t2) = (0.05, 0.10, 0.15);
        let inv = [1.0 / s(t0), 1.0 / s(t1), 1.0 / s(t2)];
        let d1 = inv[1] - inv[0];
        let d2 = inv[2] - inv[1];
        assert!((d1 - d2).abs() < 1e-9, "1/S is linear: {d1} vs {d2}");
        let ln = [s(t0).ln(), s(t1).ln(), s(t2).ln()];
        let e1 = ln[1] - ln[0];
        let e2 = ln[2] - ln[1];
        assert!((e1 - e2).abs() > 1e-3, "and ln S is not: {e1} vs {e2}");
    }

    #[test]
    fn the_enhancement_is_transient() {
        // D-13's point: enhanced early, gone later
        let c = cfg();
        let early = enhancement(atto_units::Minutes(0.1), &c).expect("early");
        let late = enhancement(atto_units::Minutes(100.0), &c).expect("late");
        assert!(early > 2.0, "a strong enhancement at first: {early}");
        assert!(late < 1.05, "and none afterwards: {late}");
        assert!(early > late, "it decays: {early} vs {late}");
    }

    #[test]
    fn a_larger_implant_gives_a_larger_integrated_enhancement() {
        let small = PairConfig { initial_supersaturation: 2.0, k_recombine: 0.5 };
        let large = PairConfig { initial_supersaturation: 20.0, k_recombine: 0.5 };
        let a = integrated_enhancement(atto_units::Minutes(10.0), &small).expect("a");
        let b = integrated_enhancement(atto_units::Minutes(10.0), &large).expect("b");
        assert!(b > a, "more supersaturation, more total diffusion: {b} vs {a}");
    }

    #[test]
    fn the_integrated_enhancement_saturates() {
        // The transient ends when the supersaturation reaches 1.0 — at (S0-1)/(S0 k) = 1.8 min
        // for these defaults — and past that point more time adds nothing. The first version of
        // this test compared 20 and 40 minutes and expected growth; both are past the end, so it
        // correctly found them equal.
        let c = cfg();
        let early = integrated_enhancement(atto_units::Minutes(0.5), &c).expect("early");
        let late = integrated_enhancement(atto_units::Minutes(1.5), &c).expect("late");
        let after = integrated_enhancement(atto_units::Minutes(100.0), &c).expect("after");
        let later = integrated_enhancement(atto_units::Minutes(200.0), &c).expect("later");
        assert!(late > early, "still growing during the transient: {late} vs {early}");
        assert!(after > late, "and more by the end: {after} vs {late}");
        assert!((after - later).abs() < 1e-12, "past the end it is spent: {after} vs {later}");
    }

    #[test]
    fn impossible_input_is_refused() {
        let c = cfg();
        assert!(supersaturation(atto_units::Minutes(-1.0), &c).is_err(), "negative time");
        assert!(supersaturation(atto_units::Minutes(f64::NAN), &c).is_err(), "not a number");
        assert!(supersaturation(atto_units::Minutes(1.0), &PairConfig { initial_supersaturation: 0.5, k_recombine: 1.0 }).is_err(), "below equilibrium");
        assert!(supersaturation(atto_units::Minutes(1.0), &PairConfig { initial_supersaturation: 10.0, k_recombine: 0.0 }).is_err(), "no recombination");
    }
}
