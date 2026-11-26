// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Atomic layer deposition: self-limiting growth (S-07).
//!
//! ALD works because each half-reaction saturates: once the surface is covered, more precursor
//! changes nothing. That gives two properties no other deposition has, and they are worth
//! stating as the model:
//!
//! * **Thickness is linear in cycles.** One cycle, one atomic layer.
//! * **Above the saturation dose, dose does not matter.** Doubling it changes nothing.
//!
//! Below saturation the growth is dose-limited instead, which is how a starved ALD process
//! behaves and why a real one runs well above the knee.
use atto_units::{Nm};

/// An ALD process.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AldProcess {
    /// Thickness per saturated cycle, nm.
    pub nm_per_cycle: f64,
    /// Precursor dose at which the surface is fully covered, arbitrary units.
    pub saturation_dose: f64,
}

impl AldProcess {
    /// Defaults for an oxide ALD.
    ///
    /// **Assumed.** `0.1 nm/cycle` is a typical saturated rate for a silica ALD; the saturation
    /// dose is in whatever units the process uses and is chosen so the knee is visible in a
    /// normal range. Registered as `ald.nm.per.cycle` and `ald.saturation.dose`.
    pub fn oxide() -> Self {
        AldProcess { nm_per_cycle: 0.1, saturation_dose: 1.0 }
    }
}

/// Why an ALD thickness could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub struct AldError {
    pub what: String,
}

impl std::fmt::Display for AldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "ALD: {}", self.what)
    }
}

/// Thickness after a number of cycles at a dose, nm.
///
/// At or above the saturation dose this is exactly `cycles x nm_per_cycle`; below it the growth
/// is scaled linearly by `dose / saturation_dose`, which is the starved regime.
pub fn thickness_nm(cycles: u32, dose: f64, p: &AldProcess) -> Result<f64, AldError> {
    if !p.nm_per_cycle.is_finite() || p.nm_per_cycle <= 0.0 {
        return Err(AldError { what: "nm per cycle must be finite and > 0".into() });
    }
    if !p.saturation_dose.is_finite() || p.saturation_dose <= 0.0 {
        return Err(AldError { what: "the saturation dose must be finite and > 0".into() });
    }
    if !dose.is_finite() || dose < 0.0 {
        return Err(AldError { what: "dose must be finite and >= 0".into() });
    }
    if cycles == 0 {
        return Ok(0.0);
    }
    let coverage = if dose >= p.saturation_dose { 1.0 } else { dose / p.saturation_dose };
    Ok(cycles as f64 * p.nm_per_cycle * coverage)
}

/// How many cycles are needed for a target thickness, at a saturated dose.
pub fn cycles_for(target_nm: Nm, p: &AldProcess) -> Result<u32, AldError> {
    let target_nm = target_nm.0;
    if !target_nm.is_finite() || target_nm < 0.0 {
        return Err(AldError { what: "the target thickness must be finite and >= 0".into() });
    }
    if !p.nm_per_cycle.is_finite() || p.nm_per_cycle <= 0.0 {
        return Err(AldError { what: "nm per cycle must be finite and > 0".into() });
    }
    Ok((target_nm / p.nm_per_cycle).ceil() as u32)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thickness_is_linear_in_cycles() {
        // S-07's first half, asserted as an identity
        let p = AldProcess::oxide();
        let a = thickness_nm(100, 5.0, &p).expect("a");
        let b = thickness_nm(200, 5.0, &p).expect("b");
        assert!((a - 10.0).abs() < 1e-12, "100 cycles at 0.1 nm is 10 nm: {a}");
        assert!((b / a - 2.0).abs() < 1e-12, "twice the cycles, twice the thickness: {}", b / a);
    }

    #[test]
    fn doubling_the_dose_above_saturation_changes_nothing() {
        // S-07's second half, and the property that defines ALD
        let p = AldProcess::oxide();
        let a = thickness_nm(100, 1.0, &p).expect("a");
        let b = thickness_nm(100, 2.0, &p).expect("b");
        let c = thickness_nm(100, 100.0, &p).expect("c");
        assert_eq!(a, b, "saturated: {a} vs {b}");
        assert_eq!(a, c, "and more dose is still saturated: {a} vs {c}");
    }

    #[test]
    fn below_saturation_the_growth_is_dose_limited() {
        let p = AldProcess::oxide();
        let starved = thickness_nm(100, 0.5, &p).expect("starved");
        let saturated = thickness_nm(100, 1.0, &p).expect("saturated");
        assert!(starved < saturated, "half the dose, less growth: {starved} vs {saturated}");
        assert!((starved / saturated - 0.5).abs() < 1e-12, "and linearly so: {}", starved / saturated);
    }

    #[test]
    fn no_cycles_is_no_film() {
        let p = AldProcess::oxide();
        assert_eq!(thickness_nm(0, 100.0, &p).expect("none"), 0.0);
    }

    #[test]
    fn the_cycle_count_for_a_target_rounds_up() {
        let p = AldProcess::oxide();
        assert_eq!(cycles_for(atto_units::Nm(10.0), &p).expect("cycles"), 100);
        assert_eq!(cycles_for(atto_units::Nm(10.05), &p).expect("cycles"), 101, "a partial cycle is still a cycle");
    }

    #[test]
    fn an_impossible_process_is_refused() {
        let p = AldProcess::oxide();
        assert!(thickness_nm(10, -1.0, &p).is_err(), "negative dose");
        assert!(thickness_nm(10, f64::NAN, &p).is_err(), "not a number");
        assert!(thickness_nm(10, 1.0, &AldProcess { nm_per_cycle: 0.0, saturation_dose: 1.0 }).is_err(), "no growth per cycle");
        assert!(cycles_for(atto_units::Nm(-1.0), &p).is_err(), "negative target");
    }
}
