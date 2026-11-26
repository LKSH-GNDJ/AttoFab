// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Stress coupling to diffusivity, oxidation and mobility (S-11).
//!
//! A film under stress is not a film under no stress with the same numbers. Compressive stress
//! slows diffusion and oxidation — the activation volume is positive, so squeezing the lattice
//! makes the jump harder — and it lowers carrier mobility. Tensile stress does the opposite.
//!
//! The coupling is exponential in the stress: `f = exp(-sigma x V / (k T))`, with `V` an
//! activation volume. That makes it a small effect at 100 MPa and a large one at 1 GPa, which is
//! the behaviour S-11 asks for and the reason capping layers work.
use atto_units::{Celsius};

/// Boltzmann constant in units that make `MPa x cm^3/mol / (R T)` dimensionless.
///
/// The gas constant in J/(mol K); stress in MPa is N/mm^2, and the activation volume is in
/// cm^3/mol, so the product needs a factor of 1e-6 to reach joules per mole.
const R_J_PER_MOL_K: f64 = 8.314_462_618;

/// Activation volumes for the three coupled quantities, cm^3/mol.
///
/// **Assumed.** S-11 asks that compressive capping retard oxidation and lower diffusivity by the
/// cited factor; the magnitudes here give a few tens of percent at 500 MPa, which is the observed
/// order. Registered as `stress.v.*`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ActivationVolumes {
    /// For diffusion, cm^3/mol.
    pub diffusion: f64,
    /// For the oxidation rate, cm^3/mol.
    pub oxidation: f64,
    /// For mobility, cm^3/mol.
    pub mobility: f64,
}

impl ActivationVolumes {
    /// The assumed activation volumes for silicon, cm^3/mol.
    pub fn silicon() -> Self {
        ActivationVolumes { diffusion: 5.0, oxidation: 8.0, mobility: 3.0 }
    }
}

/// Why a stress factor could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub struct StressCouplingError {
    pub what: String,
}

impl std::fmt::Display for StressCouplingError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "stress coupling: {}", self.what)
    }
}

fn factor(stress_mpa: f64, volume_cm3_per_mol: f64, temp_c: f64) -> Result<f64, StressCouplingError> {
    if !stress_mpa.is_finite() {
        return Err(StressCouplingError { what: "stress must be finite".into() });
    }
    if !volume_cm3_per_mol.is_finite() || volume_cm3_per_mol <= 0.0 {
        return Err(StressCouplingError { what: "the activation volume must be finite and > 0".into() });
    }
    let t_k = temp_c + 273.15;
    if !t_k.is_finite() || t_k <= 0.0 {
        return Err(StressCouplingError { what: "temperature must be above absolute zero".into() });
    }
    // MPa -> Pa is 1e6; cm^3/mol -> m^3/mol is 1e-6; the two cancel.
    //
    // The sign: compression is a NEGATIVE stress here, and it must SLOW the rate. The barrier
    // rises by pV with p = -sigma, so the factor is exp(-pV/RT) = exp(sigma V / RT), which is
    // below one for negative sigma. The first version had a minus in front and therefore made
    // compression accelerate diffusion — the test caught it.
    let exponent = stress_mpa * volume_cm3_per_mol / (R_J_PER_MOL_K * t_k);
    let f = exponent.exp();
    if !f.is_finite() || f <= 0.0 {
        return Err(StressCouplingError { what: "the factor left its domain; the stress is too large".into() });
    }
    Ok(f)
}

/// Multiplier on diffusivity from a stress, dimensionless.
///
/// Compressive stress is negative and gives a factor below one.
pub fn diffusion_factor(stress_mpa: f64, v: &ActivationVolumes, temp_c: Celsius) -> Result<f64, StressCouplingError> {
    let temp_c = temp_c.0;
    factor(stress_mpa, v.diffusion, temp_c)
}

/// Multiplier on the oxidation rate from a stress.
pub fn oxidation_factor(stress_mpa: f64, v: &ActivationVolumes, temp_c: Celsius) -> Result<f64, StressCouplingError> {
    let temp_c = temp_c.0;
    factor(stress_mpa, v.oxidation, temp_c)
}

/// Multiplier on carrier mobility from a stress.
pub fn mobility_factor(stress_mpa: f64, v: &ActivationVolumes, temp_c: Celsius) -> Result<f64, StressCouplingError> {
    let temp_c = temp_c.0;
    factor(stress_mpa, v.mobility, temp_c)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v() -> ActivationVolumes {
        ActivationVolumes::silicon()
    }

    #[test]
    fn no_stress_changes_nothing() {
        // the regression that matters: an unstressed film keeps its numbers exactly
        assert!((diffusion_factor(0.0, &v(), atto_units::Celsius(1000.0)).expect("f") - 1.0).abs() < 1e-15);
        assert!((oxidation_factor(0.0, &v(), atto_units::Celsius(1000.0)).expect("f") - 1.0).abs() < 1e-15);
        assert!((mobility_factor(0.0, &v(), atto_units::Celsius(300.0)).expect("f") - 1.0).abs() < 1e-15);
    }

    #[test]
    fn compressive_stress_retards_oxidation_and_diffusion() {
        // S-11's acceptance
        let compressive = -500.0;
        let d = diffusion_factor(compressive, &v(), atto_units::Celsius(1000.0)).expect("d");
        let o = oxidation_factor(compressive, &v(), atto_units::Celsius(1000.0)).expect("o");
        assert!(d < 1.0, "diffusion is slower: {d}");
        assert!(o < 1.0, "oxidation is slower: {o}");
        assert!(o < d, "and oxidation is retarded more, its volume being larger: {o} vs {d}");
    }

    #[test]
    fn tensile_stress_accelerates_them() {
        let tensile = 500.0;
        assert!(diffusion_factor(tensile, &v(), atto_units::Celsius(1000.0)).expect("d") > 1.0);
        assert!(oxidation_factor(tensile, &v(), atto_units::Celsius(1000.0)).expect("o") > 1.0);
    }

    #[test]
    fn the_effect_grows_with_stress() {
        let a = diffusion_factor(-200.0, &v(), atto_units::Celsius(1000.0)).expect("a");
        let b = diffusion_factor(-800.0, &v(), atto_units::Celsius(1000.0)).expect("b");
        assert!(b < a, "more compression, more retardation: {b} vs {a}");
    }

    #[test]
    fn a_hotter_wafer_couples_less() {
        // the exponent is divided by kT, so the same stress matters less when hot
        let cold = diffusion_factor(-500.0, &v(), atto_units::Celsius(700.0)).expect("cold");
        let hot = diffusion_factor(-500.0, &v(), atto_units::Celsius(1100.0)).expect("hot");
        assert!(hot > cold, "the same stress couples less at temperature: {hot} vs {cold}");
    }

    #[test]
    fn impossible_input_is_refused() {
        assert!(diffusion_factor(f64::NAN, &v(), atto_units::Celsius(1000.0)).is_err(), "not a number");
        assert!(diffusion_factor(0.0, &v(), atto_units::Celsius(-300.0)).is_err(), "below absolute zero");
        assert!(diffusion_factor(0.0, &ActivationVolumes { diffusion: 0.0, oxidation: 1.0, mobility: 1.0 }, atto_units::Celsius(1000.0)).is_err(), "no volume");
    }
}
