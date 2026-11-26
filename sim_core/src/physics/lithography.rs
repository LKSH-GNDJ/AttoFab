// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Stage 3: Lithography — Optics, Particle Scattering & Viscoelasticity Principles.
//!
//! Three interchangeable exposure modalities are modeled:
//!   - Photolithography: Rayleigh Criterion, R = k1 * lambda / NA
//!   - Electron-Beam Lithography: Monte Carlo electron trajectory kinematics
//!     with a Gaussian Proximity Effect Correction (PEC) point-spread function
//!   - Nanoimprint Lithography (NIL): viscoelastic polymer flow under pressure

use serde::{Deserialize, Serialize};

/// Which patterning modality this stage uses.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LithographyMode {
    Photo,
    EBeam,
    Nanoimprint,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PhotoInputs {
    /// Process constant k1 (typically 0.25-0.4 for advanced nodes).
    pub k1: f64,
    /// Exposure wavelength, meters (e.g. 13.5e-9 for EUV).
    pub lambda: f64,
    /// Numerical aperture of the projection optics.
    pub na: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EBeamInputs {
    /// Beam accelerating voltage, volts.
    pub voltage_v: f64,
    /// Beam current, amps.
    pub current_a: f64,
    /// Nominal beam spot / forward-scatter radius at this voltage & resist
    /// stack, meters (a simplified stand-in for a full Monte Carlo cascade).
    pub forward_scatter_radius: f64,
    /// Backscatter radius (long-range proximity blur), meters.
    pub backscatter_radius: f64,
    /// Backscatter energy fraction (eta), typically 0.1-0.5 depending on substrate Z.
    pub backscatter_eta: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct NilInputs {
    /// Imprint pressure, Pa.
    pub pressure_pa: f64,
    /// Polymer (resist) viscosity, Pa*s.
    pub viscosity: f64,
    /// Mold feature gap / channel half-height, meters.
    pub gap: f64,
    /// Imprint (dwell) time, seconds.
    pub time_s: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct LithographyOutputs {
    /// Minimum printable half-pitch feature size, meters (photo mode).
    pub rayleigh_limit_m: Option<f64>,
    /// Proximity Effect Correction dose-blur point-spread parameters (e-beam mode).
    pub pec_forward_radius_m: Option<f64>,
    pub pec_backscatter_radius_m: Option<f64>,
    pub pec_eta: Option<f64>,
    /// Predicted polymer flow front displacement, meters (NIL mode), from
    /// squeeze-flow (Stefan) lubrication approximation:
    ///   x(t) = gap * sqrt(1 + (2 * P * t) / (3 * mu) )  [simplified 1D form]
    pub nil_flow_front_m: Option<f64>,
}

#[derive(Debug, thiserror::Error)]
pub enum LithographyError {
    #[error("numerical aperture NA must be > 0, got {0}")]
    InvalidNa(f64),
    #[error("viscosity must be > 0, got {0}")]
    InvalidViscosity(f64),
    #[error("gap must be > 0, got {0}")]
    InvalidGap(f64),
    #[error("input {0} must be finite and physically valid")]
    InvalidInput(&'static str),
}

/// Rayleigh Criterion: R = k1 * lambda / NA
pub fn run_photo(inputs: &PhotoInputs) -> Result<f64, LithographyError> {
    if !inputs.na.is_finite() || inputs.na <= 0.0 {
        return Err(LithographyError::InvalidNa(inputs.na));
    }
    if !inputs.k1.is_finite() || inputs.k1 < 0.0 || !inputs.lambda.is_finite() || inputs.lambda < 0.0 {
        return Err(LithographyError::InvalidInput("Rayleigh k1/wavelength"));
    }
    Ok(inputs.k1 * inputs.lambda / inputs.na)
}

/// Simplified Monte-Carlo-informed PEC blur model. A full Monte Carlo
/// cascade traces individual electron trajectories stochastically; here we
/// expose the two-Gaussian point-spread function (forward + backscatter)
/// that such a cascade is typically fit to, so downstream dose-correction
/// code has a deterministic, reproducible model to consume.
pub fn run_ebeam(inputs: &EBeamInputs) -> (f64, f64, f64) {
    (
        inputs.forward_scatter_radius,
        inputs.backscatter_radius,
        inputs.backscatter_eta,
    )
}

/// Simplified viscoelastic squeeze-flow model for nanoimprint polymer
/// displacement under a constant applied pressure (1D lubrication
/// approximation of Stefan's equation for parallel-plate squeeze flow).
pub fn run_nil(inputs: &NilInputs) -> Result<f64, LithographyError> {
    if !inputs.viscosity.is_finite() || inputs.viscosity <= 0.0 {
        return Err(LithographyError::InvalidViscosity(inputs.viscosity));
    }
    if !inputs.gap.is_finite() || inputs.gap <= 0.0 {
        return Err(LithographyError::InvalidGap(inputs.gap));
    }
    if !inputs.pressure_pa.is_finite() || inputs.pressure_pa < 0.0 || !inputs.time_s.is_finite() || inputs.time_s < 0.0 {
        return Err(LithographyError::InvalidInput("NIL pressure/time"));
    }
    let ratio = (2.0 * inputs.pressure_pa * inputs.time_s) / (3.0 * inputs.viscosity);
    let displacement = inputs.gap * (1.0 + ratio).sqrt();
    if !displacement.is_finite() {
        return Err(LithographyError::InvalidGap(displacement));
    }
    Ok(displacement)
}

/// Run whichever lithography sub-model applies and package results into the
/// unified LithographyOutputs struct for the dependency graph.
pub fn run(
    mode: LithographyMode,
    photo: Option<&PhotoInputs>,
    ebeam: Option<&EBeamInputs>,
    nil: Option<&NilInputs>,
) -> Result<LithographyOutputs, LithographyError> {
    let mut out = LithographyOutputs {
        rayleigh_limit_m: None,
        pec_forward_radius_m: None,
        pec_backscatter_radius_m: None,
        pec_eta: None,
        nil_flow_front_m: None,
    };

    match mode {
        LithographyMode::Photo => {
            if let Some(p) = photo {
                out.rayleigh_limit_m = Some(run_photo(p)?);
            }
        }
        LithographyMode::EBeam => {
            if let Some(e) = ebeam {
                let (fwd, back, eta) = run_ebeam(e);
                out.pec_forward_radius_m = Some(fwd);
                out.pec_backscatter_radius_m = Some(back);
                out.pec_eta = Some(eta);
            }
        }
        LithographyMode::Nanoimprint => {
            if let Some(n) = nil {
                out.nil_flow_front_m = Some(run_nil(n)?);
            }
        }
    }

    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euv_gives_tighter_pitch_than_duv() {
        let euv = PhotoInputs { k1: 0.3, lambda: 13.5e-9, na: 0.33 };
        let duv = PhotoInputs { k1: 0.3, lambda: 193e-9, na: 1.35 };
        let euv_r = run_photo(&euv).unwrap();
        let duv_r = run_photo(&duv).unwrap();
        assert!(euv_r < duv_r);
    }

    #[test]
    fn nil_flow_increases_with_time() {
        let mut n = NilInputs { pressure_pa: 5e6, viscosity: 1e3, gap: 50e-9, time_s: 1.0 };
        let early = run_nil(&n).unwrap();
        n.time_s = 10.0;
        let later = run_nil(&n).unwrap();
        assert!(later > early);
    }
}
