//! Stage 4: Etching — Directional Anisotropy & Material Selectivity Principles.
//!
//!   Anisotropy:  A = 1 - R_lat / R_vert
//!   Selectivity: S = ER_target / ER_mask

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EtchingInputs {
    /// Lateral etch rate, nm/min.
    pub r_lateral: f64,
    /// Vertical etch rate, nm/min.
    pub r_vertical: f64,
    /// Target film etch rate, nm/min.
    pub er_target: f64,
    /// Mask etch rate, nm/min.
    pub er_mask: f64,
    /// RF plasma power, watts.
    pub rf_power_w: f64,
    /// Chamber pressure, mTorr.
    pub pressure_mtorr: f64,
    /// Target trench depth, nm — used to check whether the mask survives
    /// the full etch given its selectivity.
    pub target_depth_nm: f64,
    /// Mask thickness available before pattern transfer, nm.
    pub mask_thickness_nm: f64,
}

impl Default for EtchingInputs {
    fn default() -> Self {
        Self {
            r_lateral: 5.0,
            r_vertical: 50.0,
            er_target: 50.0,
            er_mask: 2.0,
            rf_power_w: 300.0,
            pressure_mtorr: 10.0,
            target_depth_nm: 500.0,
            mask_thickness_nm: 150.0,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct EtchingOutputs {
    /// Anisotropy factor A in [0, 1]; 1.0 is perfectly vertical (RIE-like),
    /// 0.0 is fully isotropic (wet-etch-like).
    pub anisotropy: f64,
    /// Selectivity ratio S = target etch rate / mask etch rate.
    pub selectivity: f64,
    /// Mask consumed during the full target-depth etch, nm.
    pub mask_consumed_nm: f64,
    /// True if the mask is fully consumed before reaching target depth.
    pub mask_erosion_failure: bool,
}

#[derive(Debug, thiserror::Error)]
pub enum EtchingError {
    #[error("vertical etch rate must be > 0, got {0}")]
    InvalidVerticalRate(f64),
    #[error("mask etch rate must be > 0, got {0}")]
    InvalidMaskRate(f64),
    #[error("input {0} must be finite and physically valid")]
    InvalidInput(&'static str),
}

pub fn run(inputs: &EtchingInputs) -> Result<EtchingOutputs, EtchingError> {
    if !inputs.r_vertical.is_finite() || inputs.r_vertical <= 0.0 {
        return Err(EtchingError::InvalidVerticalRate(inputs.r_vertical));
    }
    if !inputs.er_mask.is_finite() || inputs.er_mask <= 0.0 {
        return Err(EtchingError::InvalidMaskRate(inputs.er_mask));
    }
    for (name, value) in [("lateral etch rate", inputs.r_lateral), ("target etch rate", inputs.er_target), ("target depth", inputs.target_depth_nm), ("mask thickness", inputs.mask_thickness_nm), ("RF power", inputs.rf_power_w), ("chamber pressure", inputs.pressure_mtorr)] {
        if !value.is_finite() || value < 0.0 {
            return Err(EtchingError::InvalidInput(name));
        }
    }
    if inputs.er_target == 0.0 {
        return Err(EtchingError::InvalidInput("target etch rate"));
    }

    let anisotropy = (1.0 - (inputs.r_lateral / inputs.r_vertical)).clamp(0.0, 1.0);
    let selectivity = inputs.er_target / inputs.er_mask;

    // Time to etch to target depth, then how much mask that consumes.
    let etch_time_min = inputs.target_depth_nm / inputs.er_target;
    let mask_consumed_nm = inputs.er_mask * etch_time_min;
    let mask_erosion_failure = mask_consumed_nm >= inputs.mask_thickness_nm;

    Ok(EtchingOutputs {
        anisotropy,
        selectivity,
        mask_consumed_nm,
        mask_erosion_failure,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn purely_vertical_etch_has_unity_anisotropy() {
        let inputs = EtchingInputs { r_lateral: 0.0, ..EtchingInputs::default() };
        let out = run(&inputs).unwrap();
        assert_eq!(out.anisotropy, 1.0);
    }

    #[test]
    fn non_finite_rate_is_rejected() {
        let inputs = EtchingInputs { r_lateral: f64::NAN, ..EtchingInputs::default() };
        assert!(matches!(run(&inputs), Err(EtchingError::InvalidInput(_))));
    }

    #[test]
    fn low_selectivity_triggers_mask_erosion_failure() {
        let inputs = EtchingInputs { er_mask: 40.0, mask_thickness_nm: 100.0, ..EtchingInputs::default() };
        let out = run(&inputs).unwrap();
        assert!(out.mask_erosion_failure);
    }
}
