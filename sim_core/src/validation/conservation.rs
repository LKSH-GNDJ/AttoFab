//! Stage-level material-volume diagnostics.
//!
//! These checks compare geometric volumes only. They are not mass conservation
//! checks unless density, stoichiometry, and material conversion are explicitly
//! included by the caller. Tolerances must be selected for each model/grid.

use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VolumeSnapshot {
    pub material: String,
    pub volume_nm3: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct VolumeDelta {
    pub before_nm3: f64,
    pub after_nm3: f64,
    pub delta_nm3: f64,
    pub relative_change: f64,
}

#[derive(Debug, Error, PartialEq)]
pub enum ConservationError {
    #[error("volumes must be finite and non-negative")]
    InvalidVolume,
    #[error("relative tolerance must be finite and non-negative")]
    InvalidTolerance,
    #[error("{material} volume changed by {relative_change:.6}% which exceeds tolerance {tolerance_percent:.6}%")]
    ToleranceExceeded {
        material: String,
        relative_change: f64,
        tolerance_percent: f64,
    },
}

/// Compute a relative volume change; if initial volume is zero, absolute change
/// is normalized to 1 nm³ to avoid division by zero while remaining explicit.
pub fn volume_delta(before_nm3: f64, after_nm3: f64) -> Result<VolumeDelta, ConservationError> {
    if !before_nm3.is_finite() || !after_nm3.is_finite() || before_nm3 < 0.0 || after_nm3 < 0.0 {
        return Err(ConservationError::InvalidVolume);
    }
    let delta_nm3 = after_nm3 - before_nm3;
    let relative_change = delta_nm3.abs() / before_nm3.abs().max(1.0) * 100.0;
    Ok(VolumeDelta {
        before_nm3,
        after_nm3,
        delta_nm3,
        relative_change,
    })
}

/// Assert a material-specific relative volume tolerance. This is an explicit
/// acceptance check, not a `debug_assert!`, so release builds behave identically.
pub fn check_volume_tolerance(
    material: impl Into<String>,
    before_nm3: f64,
    after_nm3: f64,
    tolerance_percent: f64,
) -> Result<VolumeDelta, ConservationError> {
    if !tolerance_percent.is_finite() || tolerance_percent < 0.0 {
        return Err(ConservationError::InvalidTolerance);
    }
    let material = material.into();
    let delta = volume_delta(before_nm3, after_nm3)?;
    if delta.relative_change > tolerance_percent {
        return Err(ConservationError::ToleranceExceeded {
            material,
            relative_change: delta.relative_change,
            tolerance_percent,
        });
    }
    Ok(delta)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn reports_volume_delta_and_accepts_tolerance() {
        let d = check_volume_tolerance("SiO2", 100.0, 100.5, 1.0).unwrap();
        assert!((d.delta_nm3 - 0.5).abs() < 1e-12);
        assert!((d.relative_change - 0.5).abs() < 1e-12);
    }
    #[test]
    fn rejects_change_outside_tolerance() {
        assert!(matches!(
            check_volume_tolerance("Si", 100.0, 110.0, 1.0),
            Err(ConservationError::ToleranceExceeded { .. })
        ));
    }
    #[test]
    fn rejects_nonfinite_values_and_tolerance() {
        assert_eq!(
            volume_delta(f64::NAN, 1.0),
            Err(ConservationError::InvalidVolume)
        );
        assert_eq!(
            check_volume_tolerance("Si", 1.0, 1.0, f64::NAN),
            Err(ConservationError::InvalidTolerance)
        );
    }
}
