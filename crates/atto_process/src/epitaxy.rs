// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Epitaxy: graded composition and selective growth (S-06).
//!
//! An epitaxial layer copies the substrate's crystal orientation because it grows from it. That
//! works on silicon and does not work on oxide: there is no lattice to continue, so the film is
//! polycrystalline or does not grow at all. The distinction matters more than the thickness —
//! a poly layer has grain boundaries, different mobility, and no strain.
//!
//! A graded SiGe layer adds a second issue: the lattice constant changes with depth, so beyond a
//! critical thickness the strain cannot be accommodated elastically and misfit dislocations form.
//! The critical thickness falls sharply with germanium fraction, which is why graded buffers are
//! used at all.
use atto_units::{Nm};

use atto_geom::svo::Material;

/// A graded Si(1-x)Ge(x) layer.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradedSiGe {
    /// Germanium fraction at the substrate interface.
    pub ge_start: f64,
    /// Germanium fraction at the surface.
    pub ge_end: f64,
    /// Layer thickness, nm.
    pub thickness_nm: f64,
}

/// What a growth on a given surface produces.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Growth {
    /// A single-crystal layer that continues the substrate lattice.
    Single,
    /// A polycrystalline layer: no lattice to continue, but the film still grows.
    Poly,
    /// No growth at all.
    None,
}

/// Why an epitaxy calculation could not run.
#[derive(Debug, Clone, PartialEq)]
pub struct EpitaxyError {
    pub what: String,
}

impl std::fmt::Display for EpitaxyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "epitaxy: {}", self.what)
    }
}

/// Whether epitaxy on this surface is single-crystal, poly, or nothing.
///
/// Silicon continues its own lattice. An oxide or a nitride has none, so the film is poly. A metal
/// is not a growth surface for a semiconductor at all.
pub fn growth_on(surface: Material) -> Growth {
    match surface {
        Material::Silicon => Growth::Single,
        Material::SiliconDioxide | Material::SiliconNitride => Growth::Poly,
        _ => Growth::None,
    }
}

/// Germanium fraction at a depth below the surface, 0..1.
pub fn composition_at(depth_nm: Nm, g: &GradedSiGe) -> Result<f64, EpitaxyError> {
    let depth_nm = depth_nm.0;
    if !g.thickness_nm.is_finite() || g.thickness_nm <= 0.0 {
        return Err(EpitaxyError { what: "the layer thickness must be finite and > 0 nm".into() });
    }
    if !(0.0..=1.0).contains(&g.ge_start) || !(0.0..=1.0).contains(&g.ge_end) {
        return Err(EpitaxyError { what: "germanium fractions must lie in 0..=1".into() });
    }
    if !depth_nm.is_finite() || depth_nm < 0.0 || depth_nm > g.thickness_nm {
        return Err(EpitaxyError { what: "the depth must lie inside the layer".into() });
    }
    // depth is measured from the surface down; the layer grows from the substrate up
    let t = 1.0 - depth_nm / g.thickness_nm;
    Ok(g.ge_start + (g.ge_end - g.ge_start) * t)
}

/// The composition profile through the layer, as (depth, fraction) pairs.
pub fn profile(g: &GradedSiGe, samples: usize) -> Result<Vec<(f64, f64)>, EpitaxyError> {
    if samples < 2 {
        return Err(EpitaxyError { what: "a profile needs at least two samples".into() });
    }
    let mut out = Vec::with_capacity(samples);
    for i in 0..samples {
        let d = g.thickness_nm * i as f64 / (samples - 1) as f64;
        out.push((d, composition_at(atto_units::Nm(d), g)?));
    }
    Ok(out)
}

/// Critical thickness for a uniform SiGe layer, nm.
///
/// `h_c = A / x^2` with `A` chosen so that 20% germanium gives about 10 nm, which is the
/// experimentally quoted order. **Assumed** — the real criterion is Matthews-Blakeslee with the
/// elastic constants, and this is its shape rather than its calculation.
pub const CRITICAL_A_NM: f64 = 0.4;

/// Critical thickness for a uniform SiGe layer, nm.
pub fn critical_thickness_nm(ge_fraction: f64) -> Result<f64, EpitaxyError> {
    if !ge_fraction.is_finite() || !(0.0..=1.0).contains(&ge_fraction) {
        return Err(EpitaxyError { what: "the germanium fraction must lie in 0..=1".into() });
    }
    if ge_fraction <= 0.0 {
        // pure silicon has no mismatch and no critical thickness
        return Ok(f64::INFINITY);
    }
    Ok(CRITICAL_A_NM / (ge_fraction * ge_fraction))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn graded() -> GradedSiGe {
        GradedSiGe { ge_start: 0.0, ge_end: 0.3, thickness_nm: 100.0 }
    }

    #[test]
    fn a_graded_layer_reports_a_composition_profile() {
        // S-06's acceptance: a profile, not a single number
        let p = profile(&graded(), 5).expect("a profile");
        assert_eq!(p.len(), 5);
        assert!((p[0].1 - 0.3).abs() < 1e-12, "0.30 at the surface: {}", p[0].1);
        assert!((p[4].1 - 0.0).abs() < 1e-12, "0.00 at the substrate: {}", p[4].1);
        for w in p.windows(2) {
            assert!(w[1].1 <= w[0].1 + 1e-12, "the grade is monotone: {:?}", w);
        }
    }

    #[test]
    fn a_uniform_layer_has_one_composition() {
        let g = GradedSiGe { ge_start: 0.2, ge_end: 0.2, thickness_nm: 50.0 };
        for d in [0.0, 25.0, 50.0] {
            assert!((composition_at(atto_units::Nm(d), &g).expect("c") - 0.2).abs() < 1e-12, "constant at {d} nm");
        }
    }

    #[test]
    fn growth_on_oxide_is_poly_not_single_crystal() {
        // S-06's second half
        assert_eq!(growth_on(Material::Silicon), Growth::Single, "silicon continues its lattice");
        assert_eq!(growth_on(Material::SiliconDioxide), Growth::Poly, "oxide has none");
        assert_eq!(growth_on(Material::SiliconNitride), Growth::Poly, "nor has nitride");
    }

    #[test]
    fn the_critical_thickness_falls_with_germanium() {
        let a = critical_thickness_nm(0.1).expect("a");
        let b = critical_thickness_nm(0.3).expect("b");
        assert!(b < a, "more germanium, less stable strain: {b} vs {a}");
        assert!((critical_thickness_nm(0.2).expect("c") - 10.0).abs() < 1e-9, "20% gives 10 nm");
        assert!(critical_thickness_nm(0.0).expect("pure").is_infinite(), "pure silicon never relaxes");
    }

    #[test]
    fn a_depth_outside_the_layer_is_refused() {
        let g = graded();
        assert!(composition_at(atto_units::Nm(-1.0), &g).is_err(), "above the surface");
        assert!(composition_at(atto_units::Nm(101.0), &g).is_err(), "below the substrate");
        assert!(profile(&g, 1).is_err(), "a profile needs two samples");
    }

    #[test]
    fn an_impossible_layer_is_refused() {
        let mut g = graded();
        g.ge_end = 1.5;
        assert!(composition_at(atto_units::Nm(10.0), &g).is_err(), "fraction above one");
        g = graded();
        g.thickness_nm = 0.0;
        assert!(composition_at(atto_units::Nm(0.0), &g).is_err(), "no thickness");
        assert!(critical_thickness_nm(-0.1).is_err(), "negative fraction");
    }
}
