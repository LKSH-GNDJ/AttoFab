// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Stage 1: Substrate Preparation — Surface Energy & Wettability Principle.
//!
//! Governing law: Young's Equation
//!     gamma_SG = gamma_SL + gamma_LG * cos(theta)
//!
//! Rearranged to solve for the contact angle theta, which predicts whether a
//! liquid chemical or photoresist will wet the substrate without pinholes.

use serde::{Deserialize, Serialize};

/// Inputs to the Young's Equation wettability model.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PrepInputs {
    /// Solid-gas interfacial energy (mJ/m^2 or mN/m).
    pub gamma_sg: f64,
    /// Solid-liquid interfacial energy.
    pub gamma_sl: f64,
    /// Liquid-gas interfacial energy (surface tension of the liquid).
    pub gamma_lg: f64,
    /// Initial substrate RMS roughness Ra (nm) — used for a pinhole-risk heuristic.
    pub ra_nm: f64,
    /// RCA clean chemical ratio, e.g. NH4OH:H2O2:H2O, encoded as a scalar dilution factor.
    pub rca_dilution: f64,
}

impl Default for PrepInputs {
    fn default() -> Self {
        Self {
            gamma_sg: 60.0,
            gamma_sl: 20.0,
            gamma_lg: 72.8, // water at room temperature, mN/m
            ra_nm: 0.2,
            rca_dilution: 5.0,
        }
    }
}

/// Outputs of the substrate-prep stage.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PrepOutputs {
    /// Predicted contact angle theta, in radians.
    pub theta_rad: f64,
    /// Contact angle in degrees, for display.
    pub theta_deg: f64,
    /// True if the surface is considered wettable (theta < 90 deg) and
    /// downstream coating should be permitted.
    pub wettable: bool,
    /// Validation flag: blocks downstream deposition if wettability is
    /// violated (theta too large) or roughness is excessive for pinhole-free
    /// coating.
    pub validation_passed: bool,
}

#[derive(Debug, thiserror::Error)]
pub enum PrepError {
    #[error("gamma_lg (liquid-gas surface tension) must be > 0, got {0}")]
    InvalidSurfaceTension(f64),
    #[error("cos(theta) = {0} is out of the physical range [-1, 1]; inputs are thermodynamically inconsistent")]
    OutOfRange(f64),
    #[error("input {0} must be finite and non-negative")]
    InvalidInput(&'static str),
}

/// Solve Young's Equation for the equilibrium contact angle and evaluate a
/// pinhole-risk validation gate.
///
/// gamma_SG = gamma_SL + gamma_LG * cos(theta)
/// => cos(theta) = (gamma_SG - gamma_SL) / gamma_LG
pub fn run(inputs: &PrepInputs) -> Result<PrepOutputs, PrepError> {
    if !inputs.gamma_lg.is_finite() || inputs.gamma_lg <= 0.0 {
        return Err(PrepError::InvalidSurfaceTension(inputs.gamma_lg));
    }
    for (name, value) in [("solid-gas surface energy", inputs.gamma_sg), ("solid-liquid surface energy", inputs.gamma_sl), ("roughness", inputs.ra_nm), ("RCA dilution", inputs.rca_dilution)] {
        if !value.is_finite() || value < 0.0 {
            return Err(PrepError::InvalidInput(name));
        }
    }

    let cos_theta = (inputs.gamma_sg - inputs.gamma_sl) / inputs.gamma_lg;

    if !(-1.0..=1.0).contains(&cos_theta) {
        return Err(PrepError::OutOfRange(cos_theta));
    }

    let theta_rad = cos_theta.acos();
    let theta_deg = theta_rad.to_degrees();

    // Wettable (hydrophilic-style) behavior below 90 degrees.
    let wettable = theta_deg < 90.0;

    // Pinhole-risk heuristic: even a wettable surface can pinhole if
    // roughness Ra exceeds ~1 nm for sub-10nm resist films; this is a
    // simplified guard, not a full stochastic pinhole model.
    let roughness_ok = inputs.ra_nm < 1.0;

    Ok(PrepOutputs {
        theta_rad,
        theta_deg,
        wettable,
        validation_passed: wettable && roughness_ok,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn perfectly_wetting_surface_gives_zero_contact_angle() {
        let inputs = PrepInputs {
            gamma_sg: 92.8,
            gamma_sl: 20.0,
            gamma_lg: 72.8,
            ra_nm: 0.1,
            rca_dilution: 5.0,
        };
        let out = run(&inputs).unwrap();
        assert_relative_eq!(out.theta_deg, 0.0, epsilon = 1e-6);
        assert!(out.wettable);
        assert!(out.validation_passed);
    }

    #[test]
    fn out_of_range_cos_theta_is_rejected() {
        let inputs = PrepInputs {
            gamma_sg: 500.0,
            gamma_sl: 0.0,
            gamma_lg: 72.8,
            ra_nm: 0.1,
            rca_dilution: 5.0,
        };
        assert!(matches!(run(&inputs), Err(PrepError::OutOfRange(_))));
    }
}
