// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Grid Convergence Index (GCI) calculations for a scalar quantity of interest.
//!
//! Uses the three-grid Richardson extrapolation form commonly associated with
//! ASME V&V 20 practice. The caller must provide monotonically refined grids and
//! comparable results from the same physical problem.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// One grid level and its computed quantity of interest.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct GridResult {
    /// Representative grid spacing (same units for every level).
    pub spacing: f64,
    /// Quantity of interest from the simulation.
    pub value: f64,
}

/// Three-grid convergence report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GciReport {
    /// Estimated observed order of accuracy.
    pub observed_order: f64,
    /// Richardson-extrapolated value on a zero-spacing grid.
    pub extrapolated_value: f64,
    /// Fine-grid GCI as a percentage.
    pub fine_gci_percent: f64,
    /// Medium-grid GCI as a percentage.
    pub medium_gci_percent: f64,
    /// Whether the input sequence was monotonic and usable.
    pub monotonic_convergence: bool,
}

/// Errors from invalid GCI inputs.
#[derive(Debug, Error, PartialEq)]
pub enum GciError {
    /// Grid spacing or value was not finite.
    #[error("grid spacing and values must be finite")]
    NonFinite,
    /// Grid spacings must be positive and ordered coarse to fine.
    #[error("grid spacings must be positive and strictly ordered coarse to fine")]
    InvalidSpacing,
    /// Values do not support a stable three-grid estimate.
    #[error("three-grid values do not support a stable observed-order estimate")]
    UnstableValues,
    /// Safety factor must be positive.
    #[error("safety factor must be finite and positive")]
    InvalidSafetyFactor,
}

/// Compute a three-grid GCI report. `coarse`, `medium`, and `fine` must be
/// ordered from largest to smallest grid spacing. `safety_factor` is commonly
/// 1.25 when the model and refinement study justify that choice.
pub fn calculate_gci(
    coarse: GridResult,
    medium: GridResult,
    fine: GridResult,
    safety_factor: f64,
) -> Result<GciReport, GciError> {
    let grids = [coarse, medium, fine];
    if !safety_factor.is_finite() || safety_factor <= 0.0 {
        return Err(GciError::InvalidSafetyFactor);
    }
    if grids
        .iter()
        .any(|g| !g.spacing.is_finite() || !g.value.is_finite())
    {
        return Err(GciError::NonFinite);
    }
    if coarse.spacing <= medium.spacing || medium.spacing <= fine.spacing || fine.spacing <= 0.0 {
        return Err(GciError::InvalidSpacing);
    }

    let r21 = coarse.spacing / medium.spacing;
    let r32 = medium.spacing / fine.spacing;
    if (r21 - r32).abs() > 1.0e-6 * r21.max(r32).max(1.0) {
        return Err(GciError::UnstableValues);
    }

    let e21 = medium.value - coarse.value;
    let e32 = fine.value - medium.value;
    if e21.abs() < f64::EPSILON || e32.abs() < f64::EPSILON || e21.signum() != e32.signum() {
        return Err(GciError::UnstableValues);
    }
    let ratio = e21.abs() / e32.abs();
    let p = ratio.ln() / r21.ln();
    if !p.is_finite() || p <= 0.0 {
        return Err(GciError::UnstableValues);
    }

    let extrapolated = fine.value + (fine.value - medium.value) / (r32.powf(p) - 1.0);
    let fine_denominator = fine.value.abs().max(f64::EPSILON);
    let medium_denominator = medium.value.abs().max(f64::EPSILON);
    let gci_fine = safety_factor * ((fine.value - medium.value) / fine_denominator).abs()
        / (r32.powf(p) - 1.0)
        * 100.0;
    let gci_medium = safety_factor * ((medium.value - coarse.value) / medium_denominator).abs()
        / (r21.powf(p) - 1.0)
        * 100.0;

    Ok(GciReport {
        observed_order: p,
        extrapolated_value: extrapolated,
        fine_gci_percent: gci_fine,
        medium_gci_percent: gci_medium,
        monotonic_convergence: (fine.value - extrapolated).abs()
            < (medium.value - extrapolated).abs(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_second_order_sequence() {
        // Exact value 1, with error proportional to h^2.
        let report = calculate_gci(
            GridResult {
                spacing: 4.0,
                value: 1.16,
            },
            GridResult {
                spacing: 2.0,
                value: 1.04,
            },
            GridResult {
                spacing: 1.0,
                value: 1.01,
            },
            1.25,
        )
        .unwrap();
        assert!((report.observed_order - 2.0).abs() < 1.0e-8);
        assert!((report.extrapolated_value - 1.0).abs() < 1.0e-8);
    }

    #[test]
    fn rejects_unordered_spacing() {
        let g = GridResult {
            spacing: 1.0,
            value: 1.0,
        };
        assert_eq!(
            calculate_gci(g, g, g, 1.25).unwrap_err(),
            GciError::InvalidSpacing
        );
    }
}
