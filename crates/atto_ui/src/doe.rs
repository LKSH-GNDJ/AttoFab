// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The DOE tab (PRD G-18): define a sweep, see what it did.
//!
//! The panel is a presentation of `atto_study`'s sensitivity table and Monte
//! Carlo attribution — one row per knob, ordered by how much it matters, with the
//! direction of its effect. A table sorted by influence is the thing an engineer
//! reads; the raw samples are not.

use atto_study::{sensitivity, Knobs};

#[derive(Debug, Clone, PartialEq)]
/// One knob and how much it moves the target.
pub struct SensitivityRow {
    pub knob: String,
    /// Relative sensitivity of the junction depth to this knob.
    pub xj: f64,
    /// Relative sensitivity of the threshold voltage.
    pub vt: f64,
    /// The larger of the two magnitudes, used for ordering.
    pub weight: f64,
}

/// One row per knob, strongest first.
pub fn sensitivity_rows(base: &Knobs) -> Vec<SensitivityRow> {
    let mut rows: Vec<SensitivityRow> = sensitivity(base)
        .into_iter()
        .map(|(k, o)| SensitivityRow {
            knob: k.label().to_string(),
            xj: o.xj_nm,
            vt: o.vt_v,
            weight: o.xj_nm.abs().max(o.vt_v.abs()),
        })
        .collect();
    rows.sort_by(|a, b| b.weight.partial_cmp(&a.weight).unwrap_or(std::cmp::Ordering::Equal));
    rows
}

/// The DOE tab as a shell draws it: a title and the rows in order.
#[derive(Debug, Clone, PartialEq)]
pub struct DoePanel {
    pub title: String,
    pub rows: Vec<SensitivityRow>,
    pub empty_reason: Option<String>,
}

impl DoePanel {
    /// A panel seeded from the current knob values.
    pub fn from_base(base: &Knobs) -> Self {
        let rows = sensitivity_rows(base);
        let empty_reason = if rows.is_empty() {
            Some("the base case does not run, so no sensitivity can be computed".to_string())
        } else {
            None
        };
        DoePanel { title: "Sensitivity".into(), rows, empty_reason }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_knob_gets_a_row_ordered_by_influence() {
        let panel = DoePanel::from_base(&Knobs::default());
        assert!(panel.rows.len() >= 6, "the default sweep has several knobs");
        for w in panel.rows.windows(2) {
            assert!(w[0].weight >= w[1].weight, "rows must descend: {} then {}", w[0].weight, w[1].weight);
        }
        assert!(panel.rows.iter().all(|r| r.weight.is_finite()));
        assert!(panel.empty_reason.is_none());
    }

    #[test]
    fn the_strongest_knob_is_named_and_is_not_the_weakest() {
        let panel = DoePanel::from_base(&Knobs::default());
        let top = &panel.rows[0];
        let bottom = &panel.rows[panel.rows.len() - 1];
        assert!(!top.knob.is_empty());
        assert!(top.weight > bottom.weight, "{} {:.2} vs {} {:.2}", top.knob, top.weight, bottom.knob, bottom.weight);
    }

    #[test]
    fn the_rows_carry_both_outputs() {
        let panel = DoePanel::from_base(&Knobs::default());
        // the anneal temperature moves the junction and the oxide moves Vt, so at
        // least one row must show an effect on each
        assert!(panel.rows.iter().any(|r| r.xj.abs() > 1.0), "something moves xj strongly");
        assert!(panel.rows.iter().any(|r| r.vt.abs() > 1.0), "something moves Vt strongly");
    }
}
