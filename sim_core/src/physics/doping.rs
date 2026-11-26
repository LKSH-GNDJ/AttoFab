// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Stage 5: Doping — Kinetic Diffusion & Gaussian Distribution Principle.
//!
//!   N(x) = Phi / (sqrt(2*pi) * dRp) * exp( -(x - Rp)^2 / (2 * dRp^2) )

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DopingInputs {
    /// Implant dose Phi, ions/cm^2.
    pub dose_phi: f64,
    /// Projected range Rp, nm.
    pub projected_range_nm: f64,
    /// Longitudinal straggle dRp, nm.
    pub straggle_nm: f64,
    /// Accelerating voltage, kV (recorded for provenance / calibration lookups).
    pub voltage_kv: f64,
}

impl Default for DopingInputs {
    fn default() -> Self {
        Self {
            dose_phi: 1e15,
            projected_range_nm: 100.0,
            straggle_nm: 30.0,
            voltage_kv: 80.0,
        }
    }
}

#[derive(Debug, thiserror::Error)]
pub enum DopingError {
    #[error("straggle dRp must be finite and > 0, got {0}")]
    InvalidStraggle(f64),
    #[error("sparse field cell size must be finite and > 0, got {0}")]
    InvalidCellSize(f64),
    #[error("dopant concentration must be finite and >= 0, got {0}")]
    InvalidConcentration(f64),
    #[error("implant dose must be finite and >= 0, got {0}")]
    InvalidDose(f64),
    #[error("depth coordinate must be finite, got {0}")]
    InvalidDepth(f64),
    #[error("field origin coordinates must be finite")]
    InvalidOrigin,
}

/// Sparse sidecar concentration field. Geometry/material identity stays in
/// the SVO; this map stores only cells with non-zero dopant concentration.
/// Coordinates are nanometres and concentrations are ions/cm^3.
#[derive(Debug, Clone)]
pub struct SparseDopantField {
    origin_nm: (f64, f64, f64),
    cell_size_nm: f64,
    values_ions_cm3: BTreeMap<(i64, i64, i64), f64>,
}

impl SparseDopantField {
    pub fn new(origin_nm: (f64, f64, f64), cell_size_nm: f64) -> Result<Self, DopingError> {
        if !cell_size_nm.is_finite() || cell_size_nm <= 0.0 {
            return Err(DopingError::InvalidCellSize(cell_size_nm));
        }
        if !origin_nm.0.is_finite() || !origin_nm.1.is_finite() || !origin_nm.2.is_finite() {
            return Err(DopingError::InvalidOrigin);
        }
        Ok(Self { origin_nm, cell_size_nm, values_ions_cm3: BTreeMap::new() })
    }

    fn key_for(&self, p_nm: (f64, f64, f64)) -> (i64, i64, i64) {
        let idx = |p: f64, origin: f64| ((p - origin) / self.cell_size_nm).floor() as i64;
        (idx(p_nm.0, self.origin_nm.0), idx(p_nm.1, self.origin_nm.1), idx(p_nm.2, self.origin_nm.2))
    }

    pub fn set_at_nm(&mut self, p_nm: (f64, f64, f64), concentration_ions_cm3: f64) -> Result<(), DopingError> {
        if !concentration_ions_cm3.is_finite() || concentration_ions_cm3 < 0.0 {
            return Err(DopingError::InvalidConcentration(concentration_ions_cm3));
        }
        let key = self.key_for(p_nm);
        if concentration_ions_cm3 == 0.0 {
            self.values_ions_cm3.remove(&key);
        } else {
            self.values_ions_cm3.insert(key, concentration_ions_cm3);
        }
        Ok(())
    }

    /// Piecewise-constant nearest-cell lookup at a nanometre coordinate.
    pub fn concentration_at_nm(&self, p_nm: (f64, f64, f64)) -> f64 {
        self.values_ions_cm3.get(&self.key_for(p_nm)).copied().unwrap_or(0.0)
    }

    pub fn nonzero_cells(&self) -> usize { self.values_ions_cm3.len() }

    /// Approximate number of dopant atoms represented by stored cells.
    pub fn integrated_atoms(&self) -> f64 {
        let cell_edge_cm = self.cell_size_nm * 1.0e-7;
        let cell_volume_cm3 = cell_edge_cm.powi(3);
        self.values_ions_cm3.values().sum::<f64>() * cell_volume_cm3
    }
}

/// Evaluate a Gaussian implant profile at depth `x_nm` (nm).
///
/// `dose_phi` is ions/cm^2 and the returned concentration is ions/cm^3.
/// Convert the straggle to centimetres in the normalization term; the
/// Gaussian exponent remains in nm so the public depth API stays convenient.
pub fn concentration_at(inputs: &DopingInputs, x_nm: f64) -> Result<f64, DopingError> {
    if !inputs.straggle_nm.is_finite() || inputs.straggle_nm <= 0.0 {
        return Err(DopingError::InvalidStraggle(inputs.straggle_nm));
    }
    if !inputs.dose_phi.is_finite() || inputs.dose_phi < 0.0 {
        return Err(DopingError::InvalidDose(inputs.dose_phi));
    }
    if !x_nm.is_finite() || !inputs.projected_range_nm.is_finite() || !inputs.voltage_kv.is_finite() {
        return Err(DopingError::InvalidDepth(x_nm));
    }
    let sigma_nm = inputs.straggle_nm;
    let sigma_cm = sigma_nm * 1.0e-7;
    let r_p_nm = inputs.projected_range_nm;
    let prefactor = inputs.dose_phi / ((2.0 * PI).sqrt() * sigma_cm);
    let exponent = -((x_nm - r_p_nm).powi(2)) / (2.0 * sigma_nm.powi(2));
    Ok(prefactor * exponent.exp())
}

/// Sample the Gaussian profile across a depth range into a 1D array,
/// suitable for populating a column of charge-carrier density voxels.
pub fn sample_profile(
    inputs: &DopingInputs,
    x_start_nm: f64,
    x_end_nm: f64,
    n_samples: usize,
) -> Result<Vec<(f64, f64)>, DopingError> {
    if n_samples == 0 {
        return Ok(Vec::new());
    }
    let step = (x_end_nm - x_start_nm) / (n_samples.max(1) as f64 - 1.0).max(1.0);
    let mut out = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let x = x_start_nm + step * i as f64;
        out.push((x, concentration_at(inputs, x)?));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn sparse_field_stores_only_nonzero_cells() {
        let mut field = SparseDopantField::new((0.0, 0.0, 0.0), 10.0).unwrap();
        field.set_at_nm((5.0, 5.0, 5.0), 1.0e20).unwrap();
        field.set_at_nm((6.0, 6.0, 6.0), 2.0e20).unwrap();
        assert_eq!(field.nonzero_cells(), 1);
        assert_eq!(field.concentration_at_nm((7.0, 7.0, 7.0)), 2.0e20);
        field.set_at_nm((5.0, 5.0, 5.0), 0.0).unwrap();
        assert_eq!(field.nonzero_cells(), 0);
    }

    #[test]
    fn peak_concentration_occurs_at_projected_range() {
        let inputs = DopingInputs::default();
        let at_peak = concentration_at(&inputs, inputs.projected_range_nm).unwrap();
        let off_peak = concentration_at(&inputs, inputs.projected_range_nm + 50.0).unwrap();
        assert!(at_peak > off_peak);
    }

    #[test]
    fn peak_has_physical_cm3_units() {
        let inputs = DopingInputs::default();
        let expected = inputs.dose_phi / ((2.0 * PI).sqrt() * inputs.straggle_nm * 1.0e-7);
        let peak = concentration_at(&inputs, inputs.projected_range_nm).unwrap();
        assert_relative_eq!(peak, expected, max_relative = 1e-12);
        assert!(peak > 1.0e19 && peak < 1.0e21);
    }

    #[test]
    fn integrated_profile_recovers_implanted_dose() {
        let inputs = DopingInputs::default();
        let start = inputs.projected_range_nm - 8.0 * inputs.straggle_nm;
        let end = inputs.projected_range_nm + 8.0 * inputs.straggle_nm;
        let steps = 20_000usize;
        let dx_nm = (end - start) / steps as f64;
        let mut integral_ions_cm2 = 0.0;
        for i in 0..=steps {
            let x = start + i as f64 * dx_nm;
            let weight = if i == 0 || i == steps { 0.5 } else { 1.0 };
            // Convert the nm step to cm before integrating concentration.
            integral_ions_cm2 += weight * concentration_at(&inputs, x).unwrap() * dx_nm * 1.0e-7;
        }
        assert_relative_eq!(integral_ions_cm2, inputs.dose_phi, max_relative = 1e-5);
    }

    #[test]
    fn profile_is_symmetric_about_rp() {
        let inputs = DopingInputs::default();
        let left = concentration_at(&inputs, inputs.projected_range_nm - 20.0).unwrap();
        let right = concentration_at(&inputs, inputs.projected_range_nm + 20.0).unwrap();
        assert_relative_eq!(left, right, epsilon = 1e-6);
    }
}
