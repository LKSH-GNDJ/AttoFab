// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Deal-Grove oxidation (build PRD O-01, O-02).
//!
//! x^2 + A x = B (t + tau), with A and B carried to temperature by their
//! activation energies. The result drives silicon consumption and the outward
//! surface rise in the geometry operators.

use atto_units::{Kelvin, Nm};
use atto_params::{Ambient, DealGrove};

#[derive(Debug, Clone, PartialEq)]
/// Why an oxidation was refused.
pub enum DealGroveError {
    NonFiniteInput,
    InvalidTemperature,
}

impl std::fmt::Display for DealGroveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DealGroveError::NonFiniteInput => write!(f, "oxidation inputs must be finite"),
            DealGroveError::InvalidTemperature => write!(f, "temperature must be above absolute zero"),
        }
    }
}
impl std::error::Error for DealGroveError {}

/// Fraction of the grown oxide thickness consumed from the silicon.
/// Re-exported from the manifest so there is one definition (build PRD K-01).
pub use atto_params::SI_CONSUMED_RATIO;

/// Solve Deal-Grove for the oxide thickness in nm.
///
/// `initial_nm` is the oxide already present; it sets the offset `tau` so that
/// growth continues from it rather than restarting.
pub fn oxidation_nm(ambient: Ambient, temp: Kelvin, minutes: f64, initial_nm: f64) -> Result<f64, DealGroveError> {
    let temp_c = temp.celsius();
    if ![temp_c, minutes, initial_nm].iter().all(|v| v.is_finite()) {
        return Err(DealGroveError::NonFiniteInput);
    }
    if minutes < 0.0 || initial_nm < 0.0 {
        return Err(DealGroveError::NonFiniteInput);
    }
    let t_k = temp_c + 273.15;
    if !(t_k > 0.0) {
        return Err(DealGroveError::InvalidTemperature);
    }
    let dg = DealGrove::bundled(ambient);
    // Work in micrometres and hours, which is the unit system of the coefficients.
    let a = dg.a_at(t_k);
    let b = dg.b_at(t_k);
    if !(a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0) {
        return Err(DealGroveError::InvalidTemperature);
    }
    let x0 = initial_nm / 1000.0;
    // tau is the fictitious time that brings the model to x0 at t = 0.
    let tau = (x0 * x0 + a * x0) / b;
    let t = minutes / 60.0 + tau;
    // x^2 + A x - B t = 0
    let disc = a * a + 4.0 * b * t;
    let x = (-a + disc.sqrt()) / 2.0;
    Ok(x * 1000.0)
}

/// The geometry change produced by growing `oxide_nm` of oxide.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OxidationGeometry {
    /// Total oxide thickness after growth, nm.
    pub oxide_nm: f64,
    /// Silicon consumed from the substrate, nm.
    pub silicon_consumed_nm: f64,
    /// Outward rise of the outer surface, nm.
    pub surface_rise_nm: f64,
}

/// Where the surfaces end up for a given oxide thickness.
pub fn geometry_for(total_oxide: Nm, previous_oxide_nm: f64) -> OxidationGeometry {
    let total_oxide_nm = total_oxide.value();
    let grown = (total_oxide_nm - previous_oxide_nm).max(0.0);
    OxidationGeometry {
        oxide_nm: total_oxide_nm,
        silicon_consumed_nm: grown * SI_CONSUMED_RATIO,
        surface_rise_nm: grown * (1.0 - SI_CONSUMED_RATIO),
    }
}

// ---------------------------------------------------------------------------
// Oxidation as a geometry operator (build PRD O-02, O-03)
// ---------------------------------------------------------------------------

use atto_geom::column::{set_interval, solid_top, ColumnGrid};
use atto_geom::svo::Material;

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// What an oxidation grew and consumed.
pub struct OxidationStats {
    pub columns: usize,
    /// Columns that grew no oxide because a mask covered them.
    pub masked_columns: usize,
    pub mean_oxide_nm: f64,
    /// Mean outward rise of the surface, nm.
    pub mean_rise_nm: f64,
    /// Mean silicon consumed below the original interface, nm.
    pub mean_consumed_nm: f64,
    /// Silicon volume removed, nm^3.
    pub silicon_removed_nm3: f64,
    /// Oxide volume added, nm^3.
    pub oxide_added_nm3: f64,
    /// Charge recorded at every silicon/oxide boundary this step created.
    pub interface: InterfaceCharge,
}

/// Fixed charge and interface states at an oxide/silicon boundary (build PRD
/// O-05). Typical numbers for a good thermal oxide; declared as assumptions in
/// the manifest rather than measured.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InterfaceCharge {
    /// Fixed oxide charge, cm^-2 (positive for a thermal oxide).
    pub qf_cm2: f64,
    /// Interface state density, cm^-2 eV^-1.
    pub dit_cm2ev: f64,
}

impl Default for InterfaceCharge {
    fn default() -> Self {
        InterfaceCharge { qf_cm2: 1.0e11, dit_cm2ev: 1.0e10 }
    }
}

/// Lateral encroachment length as a multiple of the grown oxide thickness.
/// A bird's beak runs roughly one oxide thickness under the mask edge.
pub const BEAK_RATIO: f64 = 1.0;

/// Which columns are protected by a mask that oxidation cannot pass.
fn masked_columns(g: &ColumnGrid) -> Vec<bool> {
    let mut m = vec![false; g.cols.len()];
    for (k, col) in g.cols.iter().enumerate() {
        if !g.active[k] {
            continue;
        }
        m[k] = matches!(solid_top(col), Some((_, Material::Photoresist)) | Some((_, Material::SiliconNitride)));
    }
    m
}

/// Distance in cells from each masked column to the nearest open one, by a
/// two-pass chamfer over the active grid. `u32::MAX` means "no opening nearby".
fn distance_to_open(g: &ColumnGrid, masked: &[bool]) -> Vec<u32> {
    const FAR: u32 = u32::MAX / 4;
    let mut d = vec![FAR; g.cols.len()];
    for k in 0..g.cols.len() {
        if g.active[k] && !masked[k] {
            d[k] = 0;
        }
    }
    // forward pass
    for iy in 0..g.ny {
        for ix in 0..g.nx {
            let k = g.idx(ix, iy);
            if !g.active[k] {
                continue;
            }
            let mut best = d[k];
            if ix > 0 {
                best = best.min(d[g.idx(ix - 1, iy)].saturating_add(1));
            }
            if iy > 0 {
                best = best.min(d[g.idx(ix, iy - 1)].saturating_add(1));
            }
            d[k] = best;
        }
    }
    // backward pass
    for iy in (0..g.ny).rev() {
        for ix in (0..g.nx).rev() {
            let k = g.idx(ix, iy);
            if !g.active[k] {
                continue;
            }
            let mut best = d[k];
            if ix + 1 < g.nx {
                best = best.min(d[g.idx(ix + 1, iy)].saturating_add(1));
            }
            if iy + 1 < g.ny {
                best = best.min(d[g.idx(ix, iy + 1)].saturating_add(1));
            }
            d[k] = best;
        }
    }
    d
}

/// Grow oxide on every column and return the transformed grid.
///
/// Each column is solved with Deal-Grove starting from the oxide it already
/// carries, then the stack is rewritten: the interface moves down by 0.44 of
/// the growth and the surface up by 0.56. A mask blocks the reaction, except
/// for the lateral encroachment of a bird's beak, which decays exponentially
/// away from the opening over one oxide thickness.
pub fn oxidise_geometry(
    g: &ColumnGrid,
    ambient: Ambient,
    temp_c: f64,
    minutes: f64,
    time_fraction: f64,
) -> Result<(ColumnGrid, OxidationStats), DealGroveError> {
    oxidise_geometry_with_charge(g, ambient, temp_c, minutes, time_fraction, InterfaceCharge::default())
}

/// As [`oxidise_geometry`], with the interface charge to record.
pub fn oxidise_geometry_with_charge(
    g: &ColumnGrid,
    ambient: Ambient,
    temp_c: f64,
    minutes: f64,
    time_fraction: f64,
    charge: InterfaceCharge,
) -> Result<(ColumnGrid, OxidationStats), DealGroveError> {
    if !time_fraction.is_finite() {
        return Err(DealGroveError::NonFiniteInput);
    }
    let t = minutes * time_fraction.clamp(0.0, 1.0);
    let masked = masked_columns(g);
    let dist = distance_to_open(g, &masked);
    let mut out = g.clone();
    let mut st = OxidationStats { interface: charge, ..Default::default() };
    let area = g.cell_area();
    let (mut ox_sum, mut rise_sum, mut cons_sum, mut n) = (0.0, 0.0, 0.0, 0usize);

    for k in 0..g.cols.len() {
        if !g.active[k] {
            continue;
        }
        let col = &g.cols[k];
        let Some(si) = col.iter().rev().find(|r| r.material == Material::Silicon) else {
            continue;
        };
        let interface = si.z1;
        let existing: f64 = col
            .iter()
            .filter(|r| r.material == Material::SiliconDioxide && r.z0 >= interface - 1e-9)
            .map(|r| r.thickness())
            .sum();

        let total = oxidation_nm(ambient, Kelvin::from_celsius(temp_c), t, existing)?;
        let grown_raw = (total - existing).max(0.0);

        // A mask blocks the reaction. The beak is a lateral encroachment whose
        // length scales with the growth, so it must be computed from `grown_raw`
        // rather than from whatever oxide was already there.
        // Count a masked column as masked whether or not a sliver of beak grows.
        if masked[k] {
            st.masked_columns += 1;
        }
        let scale = if masked[k] {
            let cells = dist[k];
            if cells >= u32::MAX / 8 {
                continue;
            }
            let beak_nm = grown_raw.max(1.0) * BEAK_RATIO;
            (-(cells as f64 * g.pitch) / beak_nm).exp()
        } else {
            1.0
        };
        if scale <= 1e-6 {
            continue;
        }
        let grown = grown_raw * scale;
        if grown <= 1e-9 {
            continue;
        }
        let consumed = grown * SI_CONSUMED_RATIO;
        let rise = grown * (1.0 - SI_CONSUMED_RATIO);
        let interface_new = g.snap_z(interface - consumed);
        let surface_new = g.snap_z(interface + existing + rise);
        if surface_new <= interface_new {
            continue;
        }
        set_interval(&mut out.cols[k], interface_new, surface_new, Material::SiliconDioxide);
        let removed = (interface - interface_new).max(0.0);
        st.silicon_removed_nm3 += removed * area;
        st.oxide_added_nm3 += (surface_new - interface_new).max(0.0) * area;
        ox_sum += surface_new - interface_new;
        rise_sum += (surface_new - interface - existing).max(0.0);
        cons_sum += removed;
        n += 1;
    }
    st.columns = n;
    if n > 0 {
        st.mean_oxide_nm = ox_sum / n as f64;
        st.mean_rise_nm = rise_sum / n as f64;
        st.mean_consumed_nm = cons_sum / n as f64;
    }
    Ok((out, st))
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;
    use atto_geom::column::{grow_film, ProcessWindow};
    use atto_geom::svo::Wafer;

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        ColumnGrid::sample(&w, &ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] }).unwrap()
    }

    #[test]
    fn dry_oxide_at_1000c_matches_the_published_curve() {
        let x10 = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 10.0, 0.0).unwrap();
        let x60 = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 60.0, 0.0).unwrap();
        let x240 = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 240.0, 0.0).unwrap();
        assert!((12.0..20.0).contains(&x10), "10 min: {x10:.1} nm");
        assert!((52.0..70.0).contains(&x60), "60 min: {x60:.1} nm");
        assert!((140.0..180.0).contains(&x240), "240 min: {x240:.1} nm");
    }

    #[test]
    fn wet_oxidation_is_far_faster_than_dry() {
        let dry = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 60.0, 0.0).unwrap();
        let wet = oxidation_nm(Ambient::Wet, Kelvin::from_celsius(1000.0), 60.0, 0.0).unwrap();
        assert!(wet > 3.0 * dry, "wet {wet:.1} vs dry {dry:.1}");
    }

    #[test]
    fn growth_is_monotone_in_time_and_temperature() {
        let mut prev = 0.0;
        for m in [0.0, 1.0, 5.0, 15.0, 60.0, 240.0] {
            let x = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), m, 0.0).unwrap();
            assert!(x >= prev, "time {m}: {x} < {prev}");
            prev = x;
        }
        let mut prev = 0.0;
        for t in [800.0, 900.0, 1000.0, 1100.0] {
            let x = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(t), 60.0, 0.0).unwrap();
            assert!(x > prev, "temperature {t}: {x} <= {prev}");
            prev = x;
        }
    }

    #[test]
    fn short_time_follows_the_linear_law() {
        let a = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(900.0), 2.0, 0.0).unwrap() / 2.0;
        let b = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(900.0), 4.0, 0.0).unwrap() / 4.0;
        assert!((a - b).abs() / a < 0.15, "linear regime: {a:.3} vs {b:.3} nm/min");
    }

    #[test]
    fn long_time_approaches_the_parabolic_law() {
        let x1 = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1100.0), 240.0, 0.0).unwrap();
        let x2 = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1100.0), 960.0, 0.0).unwrap();
        let ratio = x2 / x1;
        assert!((1.8..2.15).contains(&ratio), "parabolic ratio {ratio:.3}");
    }

    #[test]
    fn starting_from_an_existing_oxide_continues_rather_than_restarting() {
        let a = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 60.0, 0.0).unwrap();
        let b = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 30.0, a).unwrap();
        let c = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), 90.0, 0.0).unwrap();
        assert!((b - c).abs() / c < 0.02, "60+30 min ({b:.2}) should equal 90 min ({c:.2})");
    }

    #[test]
    fn silicon_consumption_and_surface_rise_sum_to_the_oxide() {
        let g = geometry_for(Nm::new(100.0), 0.0);
        assert!((g.silicon_consumed_nm - 44.0).abs() < 1e-9);
        assert!((g.surface_rise_nm - 56.0).abs() < 1e-9);
        assert!((g.silicon_consumed_nm + g.surface_rise_nm - 100.0).abs() < 1e-9);
    }

    #[test]
    fn bad_inputs_are_rejected() {
        assert!(oxidation_nm(Ambient::Dry, Kelvin::from_celsius(f64::NAN), 10.0, 0.0).is_err());
        assert!(oxidation_nm(Ambient::Dry, Kelvin::from_celsius(1000.0), -1.0, 0.0).is_err());
        assert!(oxidation_nm(Ambient::Dry, Kelvin::from_celsius(-300.0), 10.0, 0.0).is_err());
    }

    // ---- geometry (O-02, O-03) -------------------------------------------

    #[test]
    fn oxidation_moves_both_interfaces_and_grows_the_stack() {
        let g = grid();
        // 240 min gives about 160 nm, thick enough that the 16 nm lattice does
        // not dominate the 44/56 split.
        let (o, st) = oxidise_geometry(&g, Ambient::Dry, 1000.0, 240.0, 1.0).unwrap();
        assert!(st.columns > 0);
        // the surface rises and the silicon interface falls
        let h0 = g.surface_heights()[0];
        let h1 = o.surface_heights()[0];
        assert!(h1 > h0, "surface rose {h0} -> {h1}");
        assert!((h1 - h0 - st.mean_rise_nm).abs() < 2.0);
        let oxide = o.cols[0].iter().filter(|r| r.material == Material::SiliconDioxide).map(|r| r.thickness()).sum::<f64>();
        assert!((oxide - st.mean_oxide_nm).abs() < 1.0, "oxide {oxide} vs {:.1}", st.mean_oxide_nm);
        // 44% of the oxide came out of the silicon, to within one cell
        let split = st.mean_consumed_nm / st.mean_oxide_nm;
        assert!((split - 0.44).abs() < 0.08, "split {split:.3}");
    }

    #[test]
    fn oxidation_conserves_the_expansion_ratio() {
        let g = grid();
        let (_, st) = oxidise_geometry(&g, Ambient::Dry, 1000.0, 240.0, 1.0).unwrap();
        // The stack grows by the expansion: (oxide - consumed) / oxide ~ 0.56.
        let expansion = (st.oxide_added_nm3 - st.silicon_removed_nm3) / st.oxide_added_nm3;
        assert!((expansion - 0.56).abs() < 0.10, "expansion fraction {expansion:.3}");
        assert!(st.oxide_added_nm3 > st.silicon_removed_nm3, "oxide volume exceeds the silicon it consumed");
    }

    #[test]
    fn oxidising_twice_equals_oxidising_once_for_longer() {
        let g = grid();
        let (a, _) = oxidise_geometry(&g, Ambient::Dry, 1000.0, 30.0, 1.0).unwrap();
        let (b, _) = oxidise_geometry(&a, Ambient::Dry, 1000.0, 30.0, 1.0).unwrap();
        let (c, _) = oxidise_geometry(&g, Ambient::Dry, 1000.0, 60.0, 1.0).unwrap();
        let ox = |gr: &ColumnGrid| {
            gr.cols[0].iter().filter(|r| r.material == Material::SiliconDioxide).map(|r| r.thickness()).sum::<f64>()
        };
        assert!((ox(&b) - ox(&c)).abs() < 2.0, "30+30 min {:.1} vs 60 min {:.1}", ox(&b), ox(&c));
    }

    #[test]
    fn a_mask_blocks_oxidation_except_at_the_edge() {
        let g = grid();
        // cover everything except a four-column window in the middle
        let mut masked = g.clone();
        for iy in 0..masked.ny {
            for ix in 0..masked.nx {
                if !(10..14).contains(&ix) {
                    let k = masked.idx(ix, iy);
                    let top = solid_top(&masked.cols[k]).unwrap().0;
                    set_interval(&mut masked.cols[k], top, top + 96.0, Material::SiliconNitride);
                }
            }
        }
        let (o, st) = oxidise_geometry(&masked, Ambient::Dry, 1000.0, 60.0, 1.0).unwrap();
        let ox = |gr: &ColumnGrid, ix: usize| {
            gr.cols[gr.idx(ix, gr.ny / 2)].iter().filter(|r| r.material == Material::SiliconDioxide).map(|r| r.thickness()).sum::<f64>()
        };
        let open = ox(&o, 11);
        let edge = ox(&o, 9);
        let far = ox(&o, 0);
        assert!(open > 20.0, "the opening oxidises: {open}");
        assert!(edge < open && edge > 0.0, "the beak under the mask edge: {edge} vs {open}");
        assert!(far < 1.0, "far under the mask nothing grows: {far}");
        assert!(st.masked_columns > 0);
    }

    #[test]
    fn a_thicker_existing_oxide_grows_less() {
        let g = grid();
        let (pre, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(96.0), 0.0).unwrap();
        let (a, _) = oxidise_geometry(&g, Ambient::Dry, 1000.0, 60.0, 1.0).unwrap();
        let (b, _) = oxidise_geometry(&pre, Ambient::Dry, 1000.0, 60.0, 1.0).unwrap();
        let ox = |gr: &ColumnGrid| {
            gr.cols[0].iter().filter(|r| r.material == Material::SiliconDioxide).map(|r| r.thickness()).sum::<f64>()
        };
        assert!(ox(&b) - ox(&pre) < ox(&a), "a 96 nm start grows less than bare silicon");
    }

    #[test]
    fn oxidation_records_the_interface_charge() {
        let g = grid();
        let charge = InterfaceCharge { qf_cm2: 5.0e11, dit_cm2ev: 2.0e11 };
        let (_, st) = oxidise_geometry_with_charge(&g, Ambient::Dry, 1000.0, 60.0, 1.0, charge).unwrap();
        assert_eq!(st.interface, charge, "the charge is carried on the result");
        assert!(st.columns > 0);
        let default = oxidise_geometry(&g, Ambient::Dry, 1000.0, 60.0, 1.0).unwrap().1;
        assert_eq!(default.interface, InterfaceCharge::default());
    }

    #[test]
    fn no_silicon_means_no_oxidation() {
        let g = grid();
        let mut hollow = g.clone();
        for k in 0..hollow.cols.len() {
            hollow.cols[k].clear();
            hollow.active[k] = false;
        }
        let (_, st) = oxidise_geometry(&hollow, Ambient::Dry, 1000.0, 60.0, 1.0).unwrap();
        assert_eq!(st.columns, 0);
    }
}
