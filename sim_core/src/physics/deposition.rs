//! Stage 2: Thin-Film Deposition — Mechanical Stress & Thermodynamics Principle.
//!
//! Governing law: Stoney Equation
//!     sigma = E / (1 - nu) * t_s^2 / (6 * t_f * R)

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DepositionInputs {
    /// Substrate Young's modulus, Pa.
    pub e_substrate: f64,
    /// Substrate Poisson's ratio (dimensionless).
    pub nu_substrate: f64,
    /// Substrate thickness, meters.
    pub t_substrate: f64,
    /// Deposited film thickness, meters.
    pub t_film: f64,
    /// Measured/predicted wafer bow radius of curvature, meters.
    pub bow_radius: f64,
    /// Film material yield strength, Pa — used as a delamination threshold.
    pub yield_strength: f64,
}

impl Default for DepositionInputs {
    fn default() -> Self {
        Self {
            e_substrate: 130e9, // Si <100>, Pa
            nu_substrate: 0.28,
            t_substrate: 500e-6, // 500 um wafer
            t_film: 1e-6,        // 1 um film
            bow_radius: 10.0,    // meters
            yield_strength: 1e9,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DepositionOutputs {
    /// Biaxial film stress, Pa (positive = tensile, negative = compressive by convention here).
    pub stress_pa: f64,
    /// True if predicted stress exceeds the film's yield strength.
    pub delamination_risk: bool,
}

#[derive(Debug, thiserror::Error)]
pub enum DepositionError {
    #[error("film thickness must be > 0, got {0}")]
    InvalidFilmThickness(f64),
    #[error("bow radius must be nonzero, got {0}")]
    InvalidBowRadius(f64),
    #[error("poisson ratio nu = {0} must be finite and in (-1, 1)")]
    InvalidPoissonRatio(f64),
    #[error("input {0} must be finite and physically valid")]
    InvalidInput(&'static str),
}

/// Solve the Stoney Equation for biaxial film stress.
pub fn run(inputs: &DepositionInputs) -> Result<DepositionOutputs, DepositionError> {
    if !inputs.t_film.is_finite() || inputs.t_film <= 0.0 {
        return Err(DepositionError::InvalidFilmThickness(inputs.t_film));
    }
    if !inputs.bow_radius.is_finite() || inputs.bow_radius == 0.0 {
        return Err(DepositionError::InvalidBowRadius(inputs.bow_radius));
    }
    if !inputs.nu_substrate.is_finite() || !(-1.0..1.0).contains(&inputs.nu_substrate) {
        return Err(DepositionError::InvalidPoissonRatio(inputs.nu_substrate));
    }
    for (name, value) in [("substrate modulus", inputs.e_substrate), ("substrate thickness", inputs.t_substrate), ("yield strength", inputs.yield_strength)] {
        if !value.is_finite() || value <= 0.0 {
            return Err(DepositionError::InvalidInput(name));
        }
    }

    let stress_pa = (inputs.e_substrate / (1.0 - inputs.nu_substrate))
        * (inputs.t_substrate.powi(2) / (6.0 * inputs.t_film * inputs.bow_radius));

    let delamination_risk = stress_pa.abs() > inputs.yield_strength;

    Ok(DepositionOutputs {
        stress_pa,
        delamination_risk,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thicker_film_reduces_stress() {
        let mut inputs = DepositionInputs::default();
        let thin = run(&inputs).unwrap().stress_pa;
        inputs.t_film *= 2.0;
        let thick = run(&inputs).unwrap().stress_pa;
        assert!(thick.abs() < thin.abs());
    }

    #[test]
    fn non_finite_inputs_are_rejected() {
        let mut inputs = DepositionInputs::default();
        inputs.bow_radius = f64::NAN;
        assert!(matches!(run(&inputs), Err(DepositionError::InvalidBowRadius(_))));
    }

    #[test]
    fn zero_film_thickness_is_rejected() {
        let mut inputs = DepositionInputs::default();
        inputs.t_film = 0.0;
        assert!(matches!(run(&inputs), Err(DepositionError::InvalidFilmThickness(_))));
    }
}
