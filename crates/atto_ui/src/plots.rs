// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Plottable series for the device tab (PRD G-14).
//!
//! A plot is data plus labels. Keeping it that way means the curve can be tested
//! without a canvas, and a shell only has to draw points it is handed.

use atto_units::{Volts};

/// One named series with its axes.
#[derive(Debug, Clone, PartialEq)]
pub struct Series {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub points: Vec<(f64, f64)>,
}

impl Series {
    /// Whether the series has any points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }
    /// The x extent, for scaling the axis.
    pub fn x_range(&self) -> Option<(f64, f64)> {
        let mut it = self.points.iter().map(|p| p.0);
        let first = it.next()?;
        Some(it.fold((first, first), |(lo, hi), v| (lo.min(v), hi.max(v))))
    }
    /// The y extent, for scaling the axis.
    pub fn y_range(&self) -> Option<(f64, f64)> {
        let mut it = self.points.iter().map(|p| p.1);
        let first = it.next()?;
        Some(it.fold((first, first), |(lo, hi), v| (lo.min(v), hi.max(v))))
    }
}

/// A C-V curve as a plottable series.
pub fn cv_series(col: &atto_device::poisson::PoissonColumn, from: Volts, to: Volts, steps: usize) -> Option<Series> {
    let (from_v, to_v) = (from.value(), to.value());
    let pts = atto_device::poisson::cv_curve(col, from_v, to_v, steps).ok()?;
    Some(Series {
        title: "C-V".into(),
        x_label: "gate voltage (V)".into(),
        y_label: "capacitance (F/cm2)".into(),
        points: pts.iter().map(|p| (p.vg, p.c)).collect(),
    })
}

/// An I-V curve as a plottable series.
pub fn iv_series(diode: &atto_device::newton::SmallDiode, biases: &[f64]) -> Option<Series> {
    let mut points = Vec::with_capacity(biases.len());
    let mut prev = None;
    for v in biases {
        let s = diode.solve(*v, prev.as_ref()).ok()?;
        points.push((*v, s.j));
        prev = Some(s);
    }
    Some(Series { title: "I-V".into(), x_label: "bias (V)".into(), y_label: "current density (A/cm2)".into(), points })
}

#[cfg(test)]
mod tests {
    use super::*;
    use atto_device::newton::SmallDiode;
    use atto_device::poisson::PoissonColumn;

    fn col() -> PoissonColumn {
        PoissonColumn { doping: vec![-1.0e15; 4000], h_nm: 1.0, tox_nm: 32.0, qf_cm2: 0.0, phi_ms: -0.55, temp_k: 300.15 }
    }

    #[test]
    fn the_cv_series_is_plottable_and_has_the_right_shape() {
        let s = cv_series(&col(), Volts::new(-1.5), Volts::new(2.0), 15).expect("a curve");
        assert_eq!(s.points.len(), 15);
        assert!(s.points.iter().all(|p| p.1.is_finite() && p.1 > 0.0), "capacitance is positive");
        assert_eq!(s.x_label, "gate voltage (V)");
        let (x0, x1) = s.x_range().unwrap();
        assert!(x0 < x1);
        // accumulation capacitance exceeds the inversion minimum
        let first = s.points[0].1;
        let min = s.points.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        assert!(first > min, "the curve dips: {first:.3e} vs {min:.3e}");
    }

    #[test]
    fn the_iv_series_rises_with_bias() {
        let s = iv_series(&SmallDiode::small(), &[0.1, 0.2, 0.3]).expect("a curve");
        assert_eq!(s.points.len(), 3);
        assert!(s.points.windows(2).all(|w| w[1].1 > w[0].1), "current rises");
        assert!(s.y_label.contains("A/cm2"));
    }

    #[test]
    fn an_empty_series_reports_no_range() {
        let s = Series { title: "x".into(), x_label: "a".into(), y_label: "b".into(), points: vec![] };
        assert!(s.is_empty());
        assert_eq!(s.x_range(), None);
        assert_eq!(s.y_range(), None);
    }
}
