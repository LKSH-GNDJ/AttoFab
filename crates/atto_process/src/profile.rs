// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Etch profile defects: bowing, notching, microtrenching (E-10).
//!
//! A deep trench that etches "correctly" by depth alone can still be unusable. Three defects
//! decide whether it is:
//!
//! * **Bowing** — the sidewall bulges outward below the mask, because the ions arrive slightly
//!   divergent and the mask erodes at its edge. The feature is wider in the middle than at the
//!   top.
//! * **Notching** — the foot undercuts sideways where the sidewall charge builds up and bends the
//!   ion trajectories. It appears at the bottom and is invisible from above.
//! * **Microtrenching** — the corners of the floor etch faster than the centre, because ions
//!   scattered off the sidewalls land there. A trench bottom that is deeper at its edges than in
//!   the middle.
//!
//! All three grow with aspect ratio, which is why they matter at exactly the dimensions modern
//! processes use. The model here is phenomenological: each defect is a fraction of the depth,
//! scaled by how far the aspect ratio exceeds the threshold where it starts.
use atto_units::{Nm};

/// The three defects, in nanometres, for one feature.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileDefects {
    /// Sidewall bulge at its widest, nm.
    pub bow_nm: f64,
    /// Lateral undercut at the foot, nm.
    pub notch_nm: f64,
    /// Extra depth at the floor corners over the centre, nm.
    pub microtrench_nm: f64,
    /// The aspect ratio they were computed at.
    pub aspect_ratio: f64,
}

impl ProfileDefects {
    /// Whether any defect is large enough to matter at a given resolution.
    pub fn is_visible(&self, resolution_nm: Nm) -> bool {
    let resolution_nm = resolution_nm.0;
        let worst = self.bow_nm.max(self.notch_nm).max(self.microtrench_nm);
        worst > resolution_nm
    }

    /// One line for a report.
    pub fn label(&self) -> String {
        format!(
            "AR {:.1}: bow {:.0} nm, notch {:.0} nm, microtrench {:.0} nm",
            self.aspect_ratio, self.bow_nm, self.notch_nm, self.microtrench_nm
        )
    }
}

/// Coefficients for the defect model, as fractions of depth per unit aspect ratio above threshold.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DefectCoefficients {
    /// Aspect ratio below which no defect appears.
    pub threshold_ar: f64,
    /// Bow fraction per unit AR above threshold.
    pub bow_per_ar: f64,
    /// Notch fraction per unit AR above threshold.
    pub notch_per_ar: f64,
    /// Microtrench fraction per unit AR above threshold.
    pub microtrench_per_ar: f64,
}

impl DefectCoefficients {
    /// Defaults for a silicon etch.
    ///
    /// **Assumed, not fitted.** E-10 asks that a trench at AR 10 show bowing and a bottom notch;
    /// these reproduce that qualitative behaviour with plausible magnitudes and are registered as
    /// `SourceKind::Assumed`. They are not a calibration.
    pub fn silicon() -> Self {
        DefectCoefficients {
            threshold_ar: 4.0,
            bow_per_ar: 0.006,
            notch_per_ar: 0.004,
            microtrench_per_ar: 0.003,
        }
    }
}

/// The defects on a feature of a given width and depth.
pub fn defects(width_nm: Nm, depth_nm: Nm, c: &DefectCoefficients) -> Result<ProfileDefects, String> {
    let depth_nm = depth_nm.0;
    let width_nm = width_nm.0;
    if !width_nm.is_finite() || width_nm <= 0.0 {
        return Err("feature width must be finite and > 0".into());
    }
    if !depth_nm.is_finite() || depth_nm < 0.0 {
        return Err("etch depth must be finite and >= 0".into());
    }
    if !c.threshold_ar.is_finite() || c.threshold_ar < 0.0 {
        return Err("the aspect-ratio threshold must be finite and >= 0".into());
    }
    for (n, v) in [("bow", c.bow_per_ar), ("notch", c.notch_per_ar), ("microtrench", c.microtrench_per_ar)] {
        if !v.is_finite() || v < 0.0 {
            return Err(format!("the {n} coefficient must be finite and >= 0"));
        }
    }
    let ar = depth_nm / width_nm;
    let excess = (ar - c.threshold_ar).max(0.0);
    Ok(ProfileDefects {
        bow_nm: depth_nm * c.bow_per_ar * excess,
        notch_nm: depth_nm * c.notch_per_ar * excess,
        microtrench_nm: depth_nm * c.microtrench_per_ar * excess,
        aspect_ratio: ar,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_deep_trench_shows_bowing_and_a_notch() {
        // E-10's acceptance: a trench at AR 10 must show both, in nanometres
        let d = defects(atto_units::Nm(20.0), atto_units::Nm(200.0), &DefectCoefficients::silicon()).expect("defects");
        assert!((d.aspect_ratio - 10.0).abs() < 1e-12, "AR 10: {}", d.aspect_ratio);
        assert!(d.bow_nm > 0.0, "bowing is present: {}", d.bow_nm);
        assert!(d.notch_nm > 0.0, "a bottom notch is present: {}", d.notch_nm);
        assert!(d.is_visible(atto_units::Nm(1.0)), "and both are visible at a nanometre: {}", d.label());
    }

    #[test]
    fn a_shallow_feature_has_none() {
        let d = defects(atto_units::Nm(200.0), atto_units::Nm(100.0), &DefectCoefficients::silicon()).expect("defects");
        assert!(d.aspect_ratio < 4.0, "below the threshold: {}", d.aspect_ratio);
        assert_eq!(d.bow_nm, 0.0);
        assert_eq!(d.notch_nm, 0.0);
        assert_eq!(d.microtrench_nm, 0.0);
    }

    #[test]
    fn every_defect_grows_with_aspect_ratio() {
        let c = DefectCoefficients::silicon();
        let a = defects(atto_units::Nm(20.0), atto_units::Nm(100.0), &c).expect("a");
        let b = defects(atto_units::Nm(20.0), atto_units::Nm(300.0), &c).expect("b");
        assert!(b.bow_nm > a.bow_nm, "bow grows: {} vs {}", b.bow_nm, a.bow_nm);
        assert!(b.notch_nm > a.notch_nm, "notch grows: {} vs {}", b.notch_nm, a.notch_nm);
        assert!(b.microtrench_nm > a.microtrench_nm, "microtrench grows: {} vs {}", b.microtrench_nm, a.microtrench_nm);
    }

    #[test]
    fn the_defects_are_fractions_of_the_depth() {
        // doubling the depth at the same width doubles every defect, since each is depth x
        // fraction x excess-AR and the excess scales with depth
        let c = DefectCoefficients::silicon();
        let a = defects(atto_units::Nm(20.0), atto_units::Nm(200.0), &c).expect("a");
        let b = defects(atto_units::Nm(20.0), atto_units::Nm(400.0), &c).expect("b");
        assert!(b.bow_nm > 2.0 * a.bow_nm, "the excess-AR term makes it more than double: {} vs {}", b.bow_nm, a.bow_nm);
    }

    #[test]
    fn a_feature_that_cannot_exist_is_refused() {
        let c = DefectCoefficients::silicon();
        assert!(defects(atto_units::Nm(0.0), atto_units::Nm(200.0), &c).is_err(), "no width");
        assert!(defects(atto_units::Nm(20.0), atto_units::Nm(-1.0), &c).is_err(), "negative depth");
        assert!(defects(atto_units::Nm(f64::NAN), atto_units::Nm(200.0), &c).is_err(), "not a number");
        let mut bad = c;
        bad.notch_per_ar = -1.0;
        assert!(defects(atto_units::Nm(20.0), atto_units::Nm(200.0), &bad).is_err(), "negative coefficient");
    }
}
