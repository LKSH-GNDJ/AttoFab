// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Stage 6: Metrology — Wave-Particle Interaction & Ellipsometry Principle.
//!
//!   rho = r_p / r_s = tan(Psi) * exp(i * Delta)
//!
//! Represented here with real-valued amplitude ratio (r_p, r_s) inputs and
//! a complex rho, from which Psi and Delta are recovered.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct MetrologyInputs {
    /// Real part of the p-polarized reflection coefficient r_p.
    pub r_p_re: f64,
    /// Imaginary part of r_p.
    pub r_p_im: f64,
    /// Real part of the s-polarized reflection coefficient r_s.
    pub r_s_re: f64,
    /// Imaginary part of r_s.
    pub r_s_im: f64,
    /// Angle of incidence, degrees (recorded for provenance).
    pub incidence_deg: f64,
}

impl Default for MetrologyInputs {
    fn default() -> Self {
        Self {
            r_p_re: 0.30,
            r_p_im: 0.10,
            r_s_re: 0.50,
            r_s_im: -0.05,
            incidence_deg: 70.0,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct MetrologyOutputs {
    /// Amplitude ratio angle Psi, degrees.
    pub psi_deg: f64,
    /// Phase difference Delta, degrees.
    pub delta_deg: f64,
    /// |rho|, magnitude of the complex ratio.
    pub rho_magnitude: f64,
}

#[derive(Debug, thiserror::Error)]
pub enum MetrologyError {
    #[error("r_s must be nonzero (division by zero in rho = r_p / r_s)")]
    ZeroRs,
    #[error("reflection coefficients must all be finite")]
    NonFiniteInput,
}

/// Solve the fundamental ellipsometry equation for Psi and Delta given
/// complex reflection coefficients r_p and r_s.
pub fn run(inputs: &MetrologyInputs) -> Result<MetrologyOutputs, MetrologyError> {
    let (rp_re, rp_im) = (inputs.r_p_re, inputs.r_p_im);
    let (rs_re, rs_im) = (inputs.r_s_re, inputs.r_s_im);
    if [rp_re, rp_im, rs_re, rs_im, inputs.incidence_deg].iter().any(|v| !v.is_finite()) {
        return Err(MetrologyError::NonFiniteInput);
    }

    let rs_mag_sq = rs_re * rs_re + rs_im * rs_im;
    if rs_mag_sq == 0.0 {
        return Err(MetrologyError::ZeroRs);
    }

    // Complex division: rho = r_p / r_s
    let rho_re = (rp_re * rs_re + rp_im * rs_im) / rs_mag_sq;
    let rho_im = (rp_im * rs_re - rp_re * rs_im) / rs_mag_sq;

    let rho_magnitude = (rho_re * rho_re + rho_im * rho_im).sqrt();
    let delta_rad = rho_im.atan2(rho_re);
    let psi_rad = rho_magnitude.atan(); // tan(Psi) = |rho|

    Ok(MetrologyOutputs {
        psi_deg: psi_rad.to_degrees(),
        delta_deg: delta_rad.to_degrees(),
        rho_magnitude,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_p_and_s_gives_45_degree_psi_and_zero_delta() {
        let inputs = MetrologyInputs {
            r_p_re: 0.4,
            r_p_im: 0.0,
            r_s_re: 0.4,
            r_s_im: 0.0,
            incidence_deg: 70.0,
        };
        let out = run(&inputs).unwrap();
        assert!((out.psi_deg - 45.0).abs() < 1e-6);
        assert!(out.delta_deg.abs() < 1e-6);
    }

    #[test]
    fn non_finite_reflection_coefficients_are_rejected() {
        let inputs = MetrologyInputs { r_p_re: f64::NAN, ..MetrologyInputs::default() };
        assert!(matches!(run(&inputs), Err(MetrologyError::NonFiniteInput)));
    }

    #[test]
    fn zero_rs_is_rejected() {
        let inputs = MetrologyInputs { r_s_re: 0.0, r_s_im: 0.0, ..MetrologyInputs::default() };
        assert!(matches!(run(&inputs), Err(MetrologyError::ZeroRs)));
    }
}
