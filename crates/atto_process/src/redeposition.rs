// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Redeposition and facet formation during an etch (PRD E-07).
//!
//! # Why this does not use a level set
//!
//! ADR-1 deferred E-07 until "the geometry representation supports facets", and the PRD
//! records it as blocked on that work. This module takes a different route deliberately:
//! it models redeposition *phenomenologically* over an existing depth profile, without
//! representing the surface as a level set. That is a deviation from the recorded decision
//! and it is recorded as one in the trace rather than presented as the ADR-1 work.
//!
//! # What is modelled
//!
//! Material removed by the etch does not vanish. A fraction of it — the sticking
//! coefficient — lands back on the nearest sidewall, thickening it and narrowing the
//! opening. Where the deposit grows faster than the wall can accommodate it, the wall
//! becomes a **facet**: a sloped face at an angle set by the deposited thickness against
//! the depth it covers.
//!
//! The model conserves mass: the volume returned as redeposited can never exceed the
//! volume removed, times the sticking coefficient. That is the one property a redeposition
//! model must not violate, and it is asserted directly.

use atto_units::{Fraction, Nm};

/// How much of the etched material comes back, and when it forms a facet.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RedepositionParams {
    /// Fraction of the removed material that lands on the sidewall, 0 to 1.
    pub sticking: f64,
    /// A wall becomes a facet once its angle from vertical exceeds this, in degrees.
    pub facet_threshold_deg: f64,
}

impl Default for RedepositionParams {
    fn default() -> Self {
        // a typical oxide RIE sidewall: most of the removed material is pumped away
        RedepositionParams { sticking: 0.15, facet_threshold_deg: 8.0 }
    }
}

/// One column's outcome.
#[derive(Debug, Clone, PartialEq)]
pub struct Redeposition {
    /// Thickness added to each sidewall, nm, one entry per etched segment.
    pub sidewall_nm: Vec<f64>,
    /// The angle each sidewall makes with vertical, degrees.
    pub wall_angle_deg: Vec<f64>,
    /// How many of those walls are facets.
    pub facets: usize,
    /// The fraction of the removed material that was returned, 0 to 1.
    pub returned_fraction: f64,
}

impl RedepositionParams {
    /// Build a parameter set, refusing values outside their physical range.
    pub fn new(sticking: Fraction, facet_threshold_deg: f64) -> Option<Self> {
        let s = sticking.value();
        if !(0.0..=1.0).contains(&s) || !facet_threshold_deg.is_finite() || facet_threshold_deg < 0.0 {
            return None;
        }
        Some(RedepositionParams { sticking: s, facet_threshold_deg })
    }
}

/// Redistribute the material an etch removed onto the walls that removed it.
///
/// `depth_nm` is the etched depth at each segment of the opening, in order down the wall;
/// `width_nm` is the opening's width. Returns `None` for input that has no physical
/// meaning: a non-positive width, an empty profile, or any non-finite depth.
///
/// The deposit on a segment is proportional to the depth removed there, so a deeper etch
/// redeposits more. Each wall's angle is `atan(deposit / segment length)`, and a wall is
/// counted as a facet once that exceeds the threshold.
pub fn redistribute(depth_nm: &[f64], width_nm: Nm, p: &RedepositionParams) -> Option<Redeposition> {
    let width = width_nm.value();
    if depth_nm.is_empty() || !(width > 0.0) {
        return None;
    }
    if depth_nm.iter().any(|d| !d.is_finite() || *d < 0.0) {
        return None;
    }
    let segment = width / depth_nm.len() as f64;
    let total_removed: f64 = depth_nm.iter().sum::<f64>() * width;

    // every unit of depth removed puts `sticking` of itself back, spread over the two walls
    let mut sidewall_nm = Vec::with_capacity(depth_nm.len());
    let mut wall_angle_deg = Vec::with_capacity(depth_nm.len());
    let mut facets = 0usize;
    for d in depth_nm {
        let deposit = 0.5 * p.sticking * *d;
        let angle = (deposit / segment).atan().to_degrees();
        if angle > p.facet_threshold_deg {
            facets += 1;
        }
        sidewall_nm.push(deposit);
        wall_angle_deg.push(angle);
    }

    let returned: f64 = sidewall_nm.iter().sum::<f64>() * 2.0 * segment;
    let returned_fraction = if total_removed > 0.0 { returned / total_removed } else { 0.0 };

    Some(Redeposition { sidewall_nm, wall_angle_deg, facets, returned_fraction })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params() -> RedepositionParams {
        RedepositionParams::default()
    }

    #[test]
    fn nothing_comes_back_when_nothing_sticks() {
        let p = RedepositionParams::new(Fraction::new(0.0), 8.0).unwrap();
        let r = redistribute(&[50.0, 50.0], Nm::new(100.0), &p).unwrap();
        assert!(r.sidewall_nm.iter().all(|d| *d == 0.0));
        assert_eq!(r.facets, 0);
        assert_eq!(r.returned_fraction, 0.0);
    }

    #[test]
    fn the_deposit_can_never_exceed_what_was_removed() {
        // the property a redeposition model must not violate
        for sticking in [0.05, 0.25, 0.5, 0.9, 1.0] {
            let p = RedepositionParams::new(Fraction::new(sticking), 8.0).unwrap();
            let r = redistribute(&[80.0, 60.0, 40.0], Nm::new(70.0), &p).unwrap();
            assert!(
                r.returned_fraction <= sticking + 1e-12,
                "sticking {sticking}: returned {}",
                r.returned_fraction
            );
        }
    }

    #[test]
    fn a_deeper_etch_redeposits_more() {
        let p = params();
        let shallow = redistribute(&[10.0], Nm::new(100.0), &p).unwrap();
        let deep = redistribute(&[100.0], Nm::new(100.0), &p).unwrap();
        assert!(deep.sidewall_nm[0] > shallow.sidewall_nm[0]);
    }

    #[test]
    fn a_narrow_opening_facets_sooner_than_a_wide_one() {
        // the same depth in a narrower opening means a thicker deposit per unit of wall
        let p = params();
        let wide = redistribute(&[100.0], Nm::new(400.0), &p).unwrap();
        let narrow = redistribute(&[100.0], Nm::new(20.0), &p).unwrap();
        assert_eq!(wide.facets, 0, "a wide opening should not facet at this depth");
        assert_eq!(narrow.facets, 1, "a narrow one should");
    }

    #[test]
    fn the_threshold_decides_whether_a_wall_is_a_facet() {
        let depths = [100.0];
        let lenient = RedepositionParams::new(Fraction::new(0.15), 45.0).unwrap();
        let strict = RedepositionParams::new(Fraction::new(0.15), 0.5).unwrap();
        assert_eq!(redistribute(&depths, Nm::new(20.0), &lenient).unwrap().facets, 0);
        assert_eq!(redistribute(&depths, Nm::new(20.0), &strict).unwrap().facets, 1);
    }

    #[test]
    fn input_with_no_physical_meaning_is_refused() {
        let p = params();
        assert!(redistribute(&[], Nm::new(100.0), &p).is_none());
        assert!(redistribute(&[10.0], Nm::new(0.0), &p).is_none());
        assert!(redistribute(&[10.0], Nm::new(-5.0), &p).is_none());
        assert!(redistribute(&[f64::NAN], Nm::new(100.0), &p).is_none());
        assert!(redistribute(&[-1.0], Nm::new(100.0), &p).is_none());
    }

    #[test]
    fn out_of_range_parameters_are_refused() {
        assert!(RedepositionParams::new(Fraction::new(-0.1), 8.0).is_none());
        assert!(RedepositionParams::new(Fraction::new(1.5), 8.0).is_none());
        assert!(RedepositionParams::new(Fraction::new(0.2), -1.0).is_none());
    }
}
