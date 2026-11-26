// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Via and contact fill: does the metal close the mouth before it reaches the bottom?
//!
//! A via is a hole through dielectric, filled with metal. Whether it fills cleanly is not a
//! matter of "enough metal": a conformal film grows on the sidewalls as fast as on the floor,
//! so a narrow deep hole pinches shut at the mouth while the bottom is still bare. What is
//! left is a void — a sealed cavity that later opens in CMP and ruins the contact.
//!
//! The model is deliberately one-dimensional: it tracks how far the floor has risen and how
//! far the two sidewalls have closed, and asks which reaches its limit first.

/// How the metal arrives.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FillMode {
    /// PVD / CVD: film grows on every surface it can see. Pinches at the mouth.
    Conformal,
    /// Electroplated or CVD with inhibitors: the floor grows, the sidewalls barely.
    BottomUp,
}

/// Geometry and process for one via.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViaFillSpec {
    /// Opening diameter, nm.
    pub diameter_nm: f64,
    /// Hole depth, nm. Aspect ratio is `depth_nm / diameter_nm`.
    pub depth_nm: f64,
    pub mode: FillMode,
    /// Sidewall growth as a fraction of the floor growth, 0..1. Conformal is 1.0.
    pub sidewall_fraction: f64,
    /// Film growth rate on a flat surface, nm/min.
    pub rate_nm_per_min: f64,
    /// Time the deposition runs, minutes.
    pub time_min: f64,
}

/// What the fill produced.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ViaFill {
    /// Fraction of the hole volume that is metal, 0..1.
    pub filled_fraction: f64,
    /// True when the mouth closed before the floor reached it.
    pub void: bool,
    /// Height of the sealed cavity, nm. Zero when there is no void.
    pub void_height_nm: f64,
    /// True when the two sidewalls met without a cavity — a seam, weaker but not hollow.
    pub seam: bool,
    /// Time at which the mouth closed, minutes. `None` when it never closed.
    pub closed_at_min: Option<f64>,
}

impl ViaFillSpec {
    /// Aspect ratio of the hole.
    pub fn aspect_ratio(&self) -> f64 {
        if self.diameter_nm > 0.0 { self.depth_nm / self.diameter_nm } else { f64::INFINITY }
    }
}

/// Fill one via.
///
/// Returns `Err` for a spec that cannot describe a hole: non-finite or non-positive
/// diameter, depth or rate, a sidewall fraction outside 0..1, or negative time. Those are
/// caller mistakes, not physics, and a number invented for them would be read as a result.
pub fn fill(s: &ViaFillSpec) -> Result<ViaFill, String> {
    let ok = |v: f64| v.is_finite() && v > 0.0;
    if !ok(s.diameter_nm) || !ok(s.depth_nm) || !ok(s.rate_nm_per_min) {
        return Err("diameter, depth and rate must be finite and > 0".into());
    }
    if !s.sidewall_fraction.is_finite() || !(0.0..=1.0).contains(&s.sidewall_fraction) {
        return Err("sidewall fraction must be within [0, 1]".into());
    }
    if !s.time_min.is_finite() || s.time_min < 0.0 {
        return Err("time must be finite and >= 0".into());
    }

    let floor_rate = s.rate_nm_per_min;
    let wall_rate = s.rate_nm_per_min * s.sidewall_fraction;
    // the mouth closes when the two sidewalls have grown half the diameter each
    let close_time = if wall_rate > 0.0 { Some(s.diameter_nm / (2.0 * wall_rate)) } else { None };

    let floor_at_close = close_time.map(|t| (floor_rate * t).min(s.depth_nm));
    let closed = close_time.filter(|t| *t <= s.time_min);

    match closed {
        // the mouth closed while the deposition was still running
        Some(t) => {
            // lint:allow no floor has closed yet, so nothing has risen: zero is the physical state
            let risen = floor_at_close.unwrap_or(0.0);
            let cavity = (s.depth_nm - risen).max(0.0);
            let volume = std::f64::consts::PI * (s.diameter_nm * 0.5).powi(2) * s.depth_nm;
            let metal = std::f64::consts::PI * (s.diameter_nm * 0.5).powi(2) * risen;
            Ok(ViaFill {
                filled_fraction: if volume > 0.0 { (metal / volume).clamp(0.0, 1.0) } else { 1.0 },
                void: cavity > 1e-9,
                void_height_nm: cavity,
                seam: cavity <= 1e-9,
                closed_at_min: Some(t),
            })
        }
        // the deposition ended before the mouth closed
        None => {
            let risen = (floor_rate * s.time_min).min(s.depth_nm);
            let gap = (s.diameter_nm - 2.0 * wall_rate * s.time_min).max(0.0);
            let volume = std::f64::consts::PI * (s.diameter_nm * 0.5).powi(2) * s.depth_nm;
            let metal = std::f64::consts::PI * (s.diameter_nm * 0.5).powi(2) * risen;
            Ok(ViaFill {
                filled_fraction: if volume > 0.0 { (metal / volume).clamp(0.0, 1.0) } else { 1.0 },
                void: false,
                void_height_nm: 0.0,
                seam: gap <= 1e-9 && risen >= s.depth_nm,
                closed_at_min: None,
            })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec(diameter_nm: f64, depth_nm: f64, mode: FillMode) -> ViaFillSpec {
        ViaFillSpec {
            diameter_nm,
            depth_nm,
            mode,
            sidewall_fraction: match mode { FillMode::Conformal => 1.0, FillMode::BottomUp => 0.05 },
            rate_nm_per_min: 100.0,
            time_min: 60.0,
        }
    }

    #[test]
    fn a_wide_shallow_via_fills_with_no_void() {
        // 200 nm across, 100 nm deep: the floor reaches the top long before the walls close
        let r = fill(&spec(200.0, 100.0, FillMode::Conformal)).expect("a fill");
        assert!(!r.void, "no cavity in a shallow hole");
        assert!((r.filled_fraction - 1.0).abs() < 1e-9, "full: {}", r.filled_fraction);
    }

    #[test]
    fn a_narrow_deep_via_pinches_and_leaves_a_void() {
        // 100 nm across, 1000 nm deep: the mouth closes at 50 nm of sidewall growth,
        // by which time the floor has risen only 50 nm of the 1000
        let r = fill(&spec(100.0, 1000.0, FillMode::Conformal)).expect("a fill");
        assert!(r.void, "a 10:1 hole must pinch off");
        assert!(r.void_height_nm > 900.0, "almost the whole hole is empty: {}", r.void_height_nm);
        assert!(r.closed_at_min.is_some(), "the mouth closed");
        assert!(r.filled_fraction < 0.1, "very little metal went in: {}", r.filled_fraction);
    }

    #[test]
    fn bottom_up_fills_the_same_hole_without_a_void() {
        // the same 10:1 geometry, filled from the floor
        let r = fill(&spec(100.0, 1000.0, FillMode::BottomUp)).expect("a fill");
        assert!(!r.void, "bottom-up does not pinch: {:?}", r);
        assert!((r.filled_fraction - 1.0).abs() < 1e-9, "full: {}", r.filled_fraction);
    }

    #[test]
    fn a_short_deposition_leaves_the_hole_part_filled() {
        let mut s = spec(200.0, 1000.0, FillMode::BottomUp);
        s.time_min = 2.0; // 200 nm of 1000
        let r = fill(&s).expect("a fill");
        assert!(!r.void, "no cavity, just unfinished");
        assert!((r.filled_fraction - 0.2).abs() < 1e-6, "a fifth full: {}", r.filled_fraction);
    }

    #[test]
    fn the_aspect_ratio_is_reported() {
        assert!((spec(100.0, 1000.0, FillMode::Conformal).aspect_ratio() - 10.0).abs() < 1e-9);
    }

    #[test]
    fn a_spec_that_cannot_describe_a_hole_is_refused() {
        let mut s = spec(100.0, 1000.0, FillMode::Conformal);
        s.diameter_nm = 0.0;
        assert!(fill(&s).is_err(), "zero diameter");
        let mut s = spec(100.0, 1000.0, FillMode::Conformal);
        s.rate_nm_per_min = f64::NAN;
        assert!(fill(&s).is_err(), "not a number");
        let mut s = spec(100.0, 1000.0, FillMode::Conformal);
        s.sidewall_fraction = 1.5;
        assert!(fill(&s).is_err(), "fraction above one");
        let mut s = spec(100.0, 1000.0, FillMode::Conformal);
        s.time_min = -1.0;
        assert!(fill(&s).is_err(), "negative time");
    }
}
