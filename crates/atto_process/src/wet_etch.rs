// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Crystallographic wet etch: KOH and TMAH (E-11).
//!
//! These etches are not isotropic and not rate-limited in the usual way: they stop on {111}
//! planes, because those planes etch roughly a hundred times slower than {100}. A square opening
//! in (100) silicon therefore does not deepen as a box — it becomes a pyramidal pit whose walls
//! are fixed at 54.74 degrees by the crystal, whatever the chemistry does.
//!
//! That angle is exact: the {111} plane makes `atan(sqrt(2))` with the (100) surface. It is the
//! one thing about this process that does not depend on a coefficient, which is why it is the
//! acceptance test.
use atto_units::{Nm};

use std::f64::consts::SQRT_2;

/// The angle a {111} sidewall makes with a (100) surface, degrees.
///
/// `atan(sqrt(2))` = 54.7356..., and it is a property of the diamond cubic lattice, not of the
/// etch. No constant is registered for it because it is not empirical.
pub fn sidewall_angle_deg() -> f64 {
    SQRT_2.atan().to_degrees()
}

/// Which crystal face the wafer presents.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Orientation {
    /// (100) — the usual CMOS surface; wet etches make pyramidal pits.
    Si100,
    /// (110) — vertical {111} walls are available, so a wet etch can cut straight down.
    Si110,
}

/// Why a wet-etch geometry could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub struct WetEtchError {
    pub what: String,
}

impl std::fmt::Display for WetEtchError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "wet etch: {}", self.what)
    }
}

/// A pit made by a mask-aligned square opening.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pit {
    /// Depth at which the {111} planes meet, nm.
    pub depth_nm: f64,
    /// Sidewall angle from the surface, degrees.
    pub angle_deg: f64,
    /// Sideways undercut at the surface, nm. Zero for (110).
    pub undercut_nm: f64,
    /// Whether the pit has come to a point (the planes met) or is still a flat-bottomed box.
    pub pointed: bool,
}

impl Pit {
    /// One line for a report.
    pub fn label(&self) -> String {
        format!(
            "pit {:.0} nm deep, walls {:.1} deg, undercut {:.0} nm{}",
            self.depth_nm,
            self.angle_deg,
            self.undercut_nm,
            if self.pointed { ", pointed" } else { ", flat bottomed" }
        )
    }
}

/// The pit a square mask opening makes in (100) silicon.
///
/// The {111} planes rise from the mask corners at 54.74 degrees, so a square of side `w` comes to
/// a point at `w/2 x sqrt(2)` below the surface. Etching longer than that changes nothing —
/// which is the property that makes KOH useful for membranes and dangerous for shallow features.
pub fn square_pit(mask_width_nm: Nm, orientation: Orientation) -> Result<Pit, WetEtchError> {
    let mask_width_nm = mask_width_nm.0;
    if !mask_width_nm.is_finite() || mask_width_nm <= 0.0 {
        return Err(WetEtchError { what: "the mask opening must be finite and > 0 nm".into() });
    }
    match orientation {
        Orientation::Si100 => {
            let depth = mask_width_nm * 0.5 * SQRT_2;
            Ok(Pit {
                depth_nm: depth,
                angle_deg: sidewall_angle_deg(),
                undercut_nm: depth / SQRT_2,
                pointed: true,
            })
        }
        Orientation::Si110 => {
            // vertical {111} walls: the etch goes straight down and never closes
            Ok(Pit { depth_nm: 0.0, angle_deg: 90.0, undercut_nm: 0.0, pointed: false })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_sidewall_angle_is_the_crystal_angle() {
        // E-11's acceptance: the correct sidewall angle, which is exact
        let a = sidewall_angle_deg();
        assert!((a - 54.7356).abs() < 1e-3, "atan(sqrt 2) in degrees: {a}");
    }

    #[test]
    fn a_square_makes_a_pyramid_of_the_right_depth() {
        // a 10 um square comes to a point at 7.071 um
        let p = square_pit(atto_units::Nm(10_000.0), Orientation::Si100).expect("a pit");
        assert!((p.depth_nm - 7_071.0678).abs() < 1e-2, "10 um square, 7.071 um deep: {}", p.depth_nm);
        assert!(p.pointed, "the planes meet");
        assert!((p.angle_deg - 54.7356).abs() < 1e-3, "at the crystal angle");
    }

    #[test]
    fn the_undercut_matches_the_angle() {
        // the undercut is depth / sqrt(2); check it against the angle directly
        let p = square_pit(atto_units::Nm(2_000.0), Orientation::Si100).expect("a pit");
        let from_angle = p.depth_nm / (p.angle_deg.to_radians().tan());
        assert!((p.undercut_nm - from_angle).abs() < 1e-6, "undercut and angle agree: {} vs {}", p.undercut_nm, from_angle);
    }

    #[test]
    fn the_depth_scales_with_the_mask_width() {
        let a = square_pit(atto_units::Nm(1_000.0), Orientation::Si100).expect("a");
        let b = square_pit(atto_units::Nm(4_000.0), Orientation::Si100).expect("b");
        assert!((b.depth_nm / a.depth_nm - 4.0).abs() < 1e-12, "four times the width, four times the depth");
    }

    #[test]
    fn a_110_wafer_can_cut_straight_down() {
        let p = square_pit(atto_units::Nm(1_000.0), Orientation::Si110).expect("a pit");
        assert_eq!(p.angle_deg, 90.0, "vertical walls are available");
        assert_eq!(p.undercut_nm, 0.0, "and there is no undercut");
        assert!(!p.pointed, "the pit does not close");
    }

    #[test]
    fn a_mask_that_cannot_exist_is_refused() {
        assert!(square_pit(atto_units::Nm(0.0), Orientation::Si100).is_err(), "no opening");
        assert!(square_pit(atto_units::Nm(-1.0), Orientation::Si100).is_err(), "negative width");
        assert!(square_pit(atto_units::Nm(f64::NAN), Orientation::Si100).is_err(), "not a number");
    }
}
