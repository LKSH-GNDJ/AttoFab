// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Oxidation-enhanced and -retarded diffusion, D-16.
//!
//! A growing oxide injects silicon interstitials into the silicon below it. Boron diffuses
//! mainly through interstitials, so it moves faster under an oxidation than under an inert
//! anneal at the same temperature; antimony, which moves through vacancies, moves slower. That
//! asymmetry is what makes LOCOS work and what makes a boron well widen under a field oxide.
//!
//! The model here is the simplest one that carries the effect: the diffusivity is multiplied by
//! `1 + k x (oxidation rate)`, with the sign of `k` deciding enhancement or retardation.
use atto_units::{Celsius, NmPerMin};

use atto_params::Species;

/// Supersaturation per unit oxidation rate, per (nm/min).
///
/// **Assumed.** D-16 asks that boron under a growing oxide diffuse faster than under an inert
/// anneal; the magnitude here is chosen to make that true and plausible, not fitted to a
/// measured OED curve. Registered as `oed.k`.
pub const OED_K: f64 = 2.0;

/// Why an OED factor could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub struct OedError {
    pub what: String,
}

impl std::fmt::Display for OedError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "OED: {}", self.what)
    }
}

/// Whether a species is interstitial- or vacancy-mediated.
///
/// Boron and phosphorus diffuse mainly through interstitials and are **enhanced** by oxidation;
/// antimony goes through vacancies and is **retarded**. That sign is the whole content of D-16.
pub fn sign(species: Species) -> f64 {
    match species {
        Species::Boron => 1.0,
        Species::Phosphorus => 1.0,
        Species::Arsenic => -0.5,
    }
}

/// The diffusivity multiplier from a growing oxide.
///
/// `1 + sign x k x rate`. Zero oxidation rate returns exactly 1.0, so an inert anneal is
/// untouched — which is the regression that matters.
pub fn oed_factor(species: Species, oxidation_rate_nm_per_min: NmPerMin) -> Result<f64, OedError> {
    let oxidation_rate_nm_per_min = oxidation_rate_nm_per_min.0;
    if !oxidation_rate_nm_per_min.is_finite() || oxidation_rate_nm_per_min < 0.0 {
        return Err(OedError { what: "the oxidation rate must be finite and >= 0 nm/min".into() });
    }
    if !OED_K.is_finite() || OED_K < 0.0 {
        return Err(OedError { what: "the OED coefficient must be finite and >= 0".into() });
    }
    let f = 1.0 + sign(species) * OED_K * oxidation_rate_nm_per_min;
    if !f.is_finite() || f <= 0.0 {
        return Err(OedError { what: "the factor left its domain; the rate is too large for this model".into() });
    }
    Ok(f)
}

/// The diffusivity with oxidation enhancement applied.
pub fn enhanced_diffusivity(
    species: Species,
    temp_c: Celsius,
    oxidation_rate_nm_per_min: NmPerMin,
) -> Result<f64, OedError> {
    let oxidation_rate_nm_per_min = oxidation_rate_nm_per_min.0;
    let temp_c = temp_c.0;
    let base = crate::anneal::ThermalHistory::diffusivity_at(species, temp_c);
    if !base.is_finite() || base <= 0.0 {
        return Err(OedError { what: "the base diffusivity is not usable at this temperature".into() });
    }
    Ok(base * oed_factor(species, atto_units::NmPerMin(oxidation_rate_nm_per_min))?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn an_inert_anneal_is_untouched() {
        // the regression that matters: no oxidation, no change
        let f = oed_factor(Species::Boron, atto_units::NmPerMin(0.0)).expect("a factor");
        assert!((f - 1.0).abs() < 1e-12, "zero rate gives exactly one: {f}");
    }

    #[test]
    fn boron_diffuses_faster_under_a_growing_oxide() {
        // D-16's acceptance
        let f = oed_factor(Species::Boron, atto_units::NmPerMin(1.0)).expect("a factor");
        assert!(f > 1.0, "boron is enhanced: {f}");
        let inert = enhanced_diffusivity(Species::Boron, atto_units::Celsius(1000.0), atto_units::NmPerMin(0.0)).expect("inert");
        let growing = enhanced_diffusivity(Species::Boron, atto_units::Celsius(1000.0), atto_units::NmPerMin(1.0)).expect("growing");
        assert!(growing > inert, "and its diffusivity is higher: {growing} vs {inert}");
    }

    #[test]
    fn arsenic_is_retarded_not_enhanced() {
        // the other half of the asymmetry: a vacancy-mediated species slows down.
        // At rate 1.0 with sign -0.5 and k = 2 the factor is exactly zero, which the model
        // refuses — so the test uses a rate inside its range and checks the sign.
        let f = oed_factor(Species::Arsenic, atto_units::NmPerMin(0.5)).expect("a factor");
        assert!(f < 1.0 && f > 0.0, "arsenic is retarded in this model: {f}");
    }

    #[test]
    fn the_enhancement_grows_with_the_oxidation_rate() {
        let a = oed_factor(Species::Boron, atto_units::NmPerMin(0.5)).expect("a");
        let b = oed_factor(Species::Boron, atto_units::NmPerMin(2.0)).expect("b");
        assert!(b > a, "a faster oxide enhances more: {b} vs {a}");
    }

    #[test]
    fn a_rate_the_model_cannot_hold_is_refused() {
        // at a large enough rate the retardation term drives the factor negative
        assert!(oed_factor(Species::Arsenic, atto_units::NmPerMin(10.0)).is_err(), "the factor would be negative");
        assert!(oed_factor(Species::Boron, atto_units::NmPerMin(-1.0)).is_err(), "negative rate");
        assert!(oed_factor(Species::Boron, atto_units::NmPerMin(f64::NAN)).is_err(), "not a number");
    }
}
