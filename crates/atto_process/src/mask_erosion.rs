// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Mask erosion and selectivity (E-12).
//!
//! An etch attacks the mask as well as the film. The mask survives if it is thick enough for the
//! depth required at the selectivity the chemistry gives — and "thick enough" is not a matter of
//! opinion: it is `depth / selectivity`.
//!
//! Two numbers decide whether a process is viable, and they are often confused:
//!
//! * **Selectivity** — film rate over mask rate. A property of the chemistry and the materials.
//! * **Depth budget** — how deep the film can be etched before the mask is gone. A property of the
//!   selectivity *and* the mask thickness.
//!
//! A process can have excellent selectivity and still fail, because the mask was too thin.
use atto_units::{Minutes, Nm, NmPerMin};

/// What an etch does to its mask.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MaskErosion {
    /// Mask thickness lost, nm.
    pub lost_nm: f64,
    /// Mask remaining, nm. Negative means the mask failed.
    pub remaining_nm: f64,
    /// Film rate over mask rate, dimensionless.
    pub selectivity: f64,
    /// Film depth this mask can support before it is gone, nm.
    pub depth_budget_nm: f64,
}

impl MaskErosion {
    /// Whether the mask survived the etch.
    pub fn survived(&self) -> bool {
        self.remaining_nm > 0.0
    }

    /// One line for a report.
    pub fn label(&self) -> String {
        format!(
            "mask lost {:.0} nm of {:.0} (selectivity {:.1}:1, budget {:.0} nm)",
            self.lost_nm,
            self.lost_nm + self.remaining_nm,
            self.selectivity,
            self.depth_budget_nm
        )
    }
}

/// Mask erosion for a given chemistry, time and starting thickness.
///
/// `film_rate` and `mask_rate` are nm/min. The loss is `mask_rate x time` — E-12's acceptance is
/// exactly that identity, and it is asserted directly rather than through a tolerance.
pub fn mask_erosion(
    mask_rate_nm_per_min: NmPerMin,
    film_rate_nm_per_min: NmPerMin,
    time_min: Minutes,
    initial_nm: Nm,
) -> Result<MaskErosion, String> {
    let film_rate_nm_per_min = film_rate_nm_per_min.0;
    let initial_nm = initial_nm.0;
    let mask_rate_nm_per_min = mask_rate_nm_per_min.0;
    let time_min = time_min.0;
    if !mask_rate_nm_per_min.is_finite() || mask_rate_nm_per_min < 0.0 {
        return Err("the mask rate must be finite and >= 0 nm/min".into());
    }
    if !film_rate_nm_per_min.is_finite() || film_rate_nm_per_min < 0.0 {
        return Err("the film rate must be finite and >= 0 nm/min".into());
    }
    if !time_min.is_finite() || time_min < 0.0 {
        return Err("the etch time must be finite and >= 0 min".into());
    }
    if !initial_nm.is_finite() || initial_nm <= 0.0 {
        return Err("the mask must start with a positive thickness".into());
    }
    let lost = mask_rate_nm_per_min * time_min;
    // selectivity is film over mask; an infinite mask rate means nothing is etched selectively,
    // and an infinite film rate means the mask is never touched — both are reported as they are
    let selectivity = if mask_rate_nm_per_min > 0.0 { film_rate_nm_per_min / mask_rate_nm_per_min } else { f64::INFINITY };
    let budget = if selectivity.is_finite() { initial_nm * selectivity } else { f64::INFINITY };
    Ok(MaskErosion {
        lost_nm: lost,
        remaining_nm: initial_nm - lost,
        selectivity,
        depth_budget_nm: budget,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_loss_is_rate_times_time() {
        // E-12's acceptance, asserted as an identity rather than a tolerance
        let e = mask_erosion(atto_units::NmPerMin(2.0), atto_units::NmPerMin(40.0), atto_units::Minutes(3.0), atto_units::Nm(100.0)).expect("erosion");
        assert!((e.lost_nm - 6.0).abs() < 1e-12, "2 nm/min for 3 min is 6 nm: {}", e.lost_nm);
        assert!((e.remaining_nm - 94.0).abs() < 1e-12, "and 94 remains: {}", e.remaining_nm);
    }

    #[test]
    fn selectivity_is_film_over_mask() {
        let e = mask_erosion(atto_units::NmPerMin(2.0), atto_units::NmPerMin(40.0), atto_units::Minutes(1.0), atto_units::Nm(100.0)).expect("erosion");
        assert!((e.selectivity - 20.0).abs() < 1e-12, "40 over 2 is 20:1: {}", e.selectivity);
    }

    #[test]
    fn the_depth_budget_is_the_mask_times_the_selectivity() {
        // a good selectivity does not save a thin mask
        let e = mask_erosion(atto_units::NmPerMin(2.0), atto_units::NmPerMin(40.0), atto_units::Minutes(1.0), atto_units::Nm(50.0)).expect("erosion");
        assert!((e.depth_budget_nm - 1000.0).abs() < 1e-9, "50 nm at 20:1 budgets 1000 nm: {}", e.depth_budget_nm);
    }

    #[test]
    fn a_mask_that_is_too_thin_fails_and_says_so() {
        // 3 nm/min for 40 min is 120 nm from a 100 nm mask
        let e = mask_erosion(atto_units::NmPerMin(3.0), atto_units::NmPerMin(30.0), atto_units::Minutes(40.0), atto_units::Nm(100.0)).expect("erosion");
        assert!(!e.survived(), "the mask is gone: {} nm left", e.remaining_nm);
        assert!(e.remaining_nm < 0.0, "and the report says by how much: {}", e.remaining_nm);
    }

    #[test]
    fn a_mask_that_is_never_attacked_has_no_budget_limit() {
        let e = mask_erosion(atto_units::NmPerMin(0.0), atto_units::NmPerMin(40.0), atto_units::Minutes(10.0), atto_units::Nm(100.0)).expect("erosion");
        assert_eq!(e.lost_nm, 0.0, "an inert mask loses nothing");
        assert!(e.selectivity.is_infinite(), "and the selectivity is unbounded");
        assert!(e.survived());
    }

    #[test]
    fn a_recipe_that_cannot_exist_is_refused() {
        assert!(mask_erosion(atto_units::NmPerMin(-1.0), atto_units::NmPerMin(40.0), atto_units::Minutes(1.0), atto_units::Nm(100.0)).is_err(), "negative mask rate");
        assert!(mask_erosion(atto_units::NmPerMin(2.0), atto_units::NmPerMin(-1.0), atto_units::Minutes(1.0), atto_units::Nm(100.0)).is_err(), "negative film rate");
        assert!(mask_erosion(atto_units::NmPerMin(2.0), atto_units::NmPerMin(40.0), atto_units::Minutes(-1.0), atto_units::Nm(100.0)).is_err(), "negative time");
        assert!(mask_erosion(atto_units::NmPerMin(2.0), atto_units::NmPerMin(40.0), atto_units::Minutes(1.0), atto_units::Nm(0.0)).is_err(), "no mask");
        assert!(mask_erosion(atto_units::NmPerMin(f64::NAN), atto_units::NmPerMin(40.0), atto_units::Minutes(1.0), atto_units::Nm(100.0)).is_err(), "not a number");
    }
}
