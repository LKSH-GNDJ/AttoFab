// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Activation and solid solubility (D-14).
//!
//! A doping step puts atoms in the silicon. An anneal decides how many of them sit on lattice
//! sites and donate a carrier. Above the solid solubility limit the excess precipitates and
//! does nothing electrically, so **total dose and active dose are different numbers** — and it
//! is the active one that sheet resistance and every device metric actually see.
//!
//! A model that reports the implant dose as though it were active is wrong by however much of
//! the profile sits above the limit, which for a 1e20 implant is most of the peak.

use atto_params::Species;
use atto_units::{Kelvin, PerCm3};

/// Boltzmann constant, eV/K.
const KB_EV: f64 = 8.617_333_262e-5;

/// A solubility limit expressed as `N0 exp(-Ea / kT)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Solubility {
    /// Pre-exponential, cm^-3.
    pub n0_cm3: f64,
    /// Activation energy, eV.
    pub ea_ev: f64,
}

impl Solubility {
    /// The limit for a species.
    ///
    /// **These are assumed, not fitted.** D-14 asks for a 1e20 cm^-3 implant annealed at
    /// 1000 C to report active < total; the pairs below are chosen to reproduce the accepted
    /// ordering (B < As < P at a given temperature) and the right order of magnitude, not a
    /// particular published table. Calibrating them against one is a separate block.
    pub fn of(species: Species) -> Self {
        match species {
            // boron is the least soluble of the three at 1000 C
            Species::Boron => Solubility { n0_cm3: 1.2e22, ea_ev: 1.05 },
            Species::Arsenic => Solubility { n0_cm3: 3.0e22, ea_ev: 1.02 },
            Species::Phosphorus => Solubility { n0_cm3: 1.5e23, ea_ev: 1.10 },
        }
    }

    /// The limit at a temperature.
    pub fn at(&self, temp: Kelvin) -> Result<PerCm3, ActivationError> {
        let t = temp.0;
        if !t.is_finite() || t <= 0.0 {
            return Err(ActivationError::BadTemperature);
        }
        if !self.n0_cm3.is_finite() || self.n0_cm3 <= 0.0 || !self.ea_ev.is_finite() {
            return Err(ActivationError::BadSolubility);
        }
        let n = self.n0_cm3 * (-self.ea_ev / (KB_EV * t)).exp();
        if !n.is_finite() || n <= 0.0 {
            return Err(ActivationError::BadSolubility);
        }
        Ok(PerCm3(n))
    }
}

/// Why an activation calculation could not run.
#[derive(Debug, Clone, PartialEq)]
pub enum ActivationError {
    /// The temperature was not a positive absolute temperature.
    BadTemperature,
    /// A solubility pair was not usable at this temperature.
    BadSolubility,
    /// The total concentration was not a positive finite number.
    BadConcentration,
}

impl std::fmt::Display for ActivationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ActivationError::BadTemperature => write!(f, "temperature must be a positive absolute temperature"),
            ActivationError::BadSolubility => write!(f, "the solubility limit is not usable at this temperature"),
            ActivationError::BadConcentration => write!(f, "total concentration must be finite and > 0"),
        }
    }
}

/// What an anneal makes electrically active.
#[derive(Debug, Clone, PartialEq)]
pub struct Activation {
    /// Total concentration present.
    pub total: PerCm3,
    /// Concentration on lattice sites, at most the solubility limit.
    pub active: PerCm3,
    /// The limit itself, so a report can say what capped it.
    pub limit: PerCm3,
    /// Active fraction, 0..1.
    pub fraction: f64,
    /// True when the total exceeded the limit and the excess precipitated.
    pub capped: bool,
}

/// The active concentration after an anneal at `temp`.
///
/// Everything above the solubility limit is assumed electrically dead. That is the standard
/// first-order treatment and it is wrong in detail — clustering and deactivation kinetics mean
/// the active fraction depends on the anneal history, not only its peak temperature — but it
/// captures the one thing D-14 asks for: active is not total.
pub fn activate(total: PerCm3, species: Species, temp: Kelvin) -> Result<Activation, ActivationError> {
    if !total.0.is_finite() || total.0 <= 0.0 {
        return Err(ActivationError::BadConcentration);
    }
    let limit = Solubility::of(species).at(temp)?;
    let active = total.0.min(limit.0);
    Ok(Activation {
        total,
        active: PerCm3(active),
        limit,
        fraction: active / total.0,
        capped: total.0 > limit.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn at_1000c() -> Kelvin {
        Kelvin::new(1273.15)
    }

    #[test]
    fn a_heavy_implant_reports_less_active_than_total() {
        // D-14's acceptance: 1e20 cm^-3 annealed at 1000 C
        let a = activate(PerCm3(1.0e20), Species::Boron, at_1000c()).expect("an activation");
        assert!(a.active.0 < a.total.0, "active must be below total: {} vs {}", a.active.0, a.total.0);
        assert!(a.capped, "the excess precipitated");
        assert!(a.fraction < 1.0 && a.fraction > 0.0, "the fraction is a fraction: {}", a.fraction);
        assert!((a.active.0 - a.limit.0).abs() < 1e-9, "and it sits at the limit");
    }

    #[test]
    fn a_light_implant_is_fully_active() {
        let a = activate(PerCm3(1.0e17), Species::Boron, at_1000c()).expect("an activation");
        assert!(!a.capped, "nothing precipitates below the limit");
        assert!((a.active.0 - 1.0e17).abs() < 1.0, "all of it is active");
        assert!((a.fraction - 1.0).abs() < 1e-12, "fraction one");
    }

    #[test]
    fn solubility_rises_with_temperature() {
        let cold = Solubility::of(Species::Arsenic).at(Kelvin::new(1073.15)).expect("1000 K").0;
        let hot = Solubility::of(Species::Arsenic).at(Kelvin::new(1373.15)).expect("1100 C").0;
        assert!(hot > cold, "a hotter anneal dissolves more: {} vs {}", hot, cold);
    }

    #[test]
    fn the_species_order_is_boron_below_arsenic_below_phosphorus() {
        // the accepted ordering at a given temperature, which is the one thing the assumed
        // pairs are chosen to reproduce
        let t = at_1000c();
        let b = Solubility::of(Species::Boron).at(t).expect("B").0;
        let a = Solubility::of(Species::Arsenic).at(t).expect("As").0;
        let p = Solubility::of(Species::Phosphorus).at(t).expect("P").0;
        assert!(b < a, "boron below arsenic: {b} vs {a}");
        assert!(a < p, "arsenic below phosphorus: {a} vs {p}");
    }

    #[test]
    fn a_recipe_that_cannot_exist_is_refused() {
        assert!(activate(PerCm3(0.0), Species::Boron, at_1000c()).is_err(), "no dopant");
        assert!(activate(PerCm3(f64::NAN), Species::Boron, at_1000c()).is_err(), "not a number");
        assert!(activate(PerCm3(1.0e18), Species::Boron, Kelvin::new(0.0)).is_err(), "zero kelvin");
        assert!(activate(PerCm3(1.0e18), Species::Boron, Kelvin::new(-5.0)).is_err(), "negative kelvin");
    }
}
