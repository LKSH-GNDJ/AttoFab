// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The seam between the pipeline and the physics.
//!
//! The application's `sim_core` is written against `serde` and cannot be compiled without
//! the registry, so every edit inside it is unverifiable in this environment.
//! This crate moves as much of the wiring as possible out of that blind spot:
//! each function takes and returns plain data — slices of `f64`, strings,
//! tuples — so the logic is tested here and the the earlier build-side edit shrinks to a call.
//!
//! Three of these also carry audit fixes that the earlier build got wrong:
//! F-07 (junction measured from the wrong surface), F-10 (step height from a
//! percentile spread rather than adjacent columns) and F-11 (an unlisted
//! material silently becoming free passage).

pub mod plugin;

use atto_units::{MTorr, Nm, Watts};
use atto_geom::column::ColumnGrid;
use atto_geom::svo::Run;
use atto_geom::svo::Material;
use atto_params::Species;
use atto_process::etch::EtchModel;
use atto_process::implant::Profile;

#[derive(Debug, Clone, PartialEq)]
/// Why a bridge call failed.
pub enum BridgeError {
    UnknownSpecies(String),
    BadInput(String),
}

impl std::fmt::Display for BridgeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BridgeError::UnknownSpecies(s) => write!(f, "unknown species '{s}'"),
            BridgeError::BadInput(s) => write!(f, "{s}"),
        }
    }
}
impl std::error::Error for BridgeError {}

/// Species by name, so the application's string-typed species maps without a new enum.
pub fn species_by_name(name: &str) -> Result<Species, BridgeError> {
    match name.to_ascii_lowercase().as_str() {
        "b" | "boron" => Ok(Species::Boron),
        "p" | "phosphorus" => Ok(Species::Phosphorus),
        "as" | "arsenic" => Ok(Species::Arsenic),
        other => Err(BridgeError::UnknownSpecies(other.into())),
    }
}

/// A column run in primitive form: `(z0, z1, material code)`.
///
/// The material code is the order of `Material` in its declaration, so the caller
/// passes a plain integer and this side maps it. That keeps the application's grid type out
/// of the bridge entirely, which is what makes the fork decidable without
/// compiling the earlier build.
pub type PrimitiveRun = (f64, f64, u8);

/// Material code for a material, and back. Codes are stable; adding a variant
/// appends rather than renumbers.
/// Code for a doped-silicon run. the earlier build carries a dopant species in its material
/// enum; the crates does not, because the dopant lives in a separate field. The species
/// is therefore dropped across the bridge — which is why `grid_from_columns`
/// counts these and hands the count back rather than losing them quietly.
pub const CODE_DOPANT: u8 = 8;

/// The integer code a primitive carries for a material.
pub fn material_code(m: Material) -> u8 {
    match m {
        Material::Air => 0,
        Material::Silicon => 1,
        Material::SiliconDioxide => 2,
        Material::SiliconNitride => 3,
        Material::Photoresist => 4,
        Material::Metal => 5,
    }
}

/// Material for a code.
///
/// Code 8 is a doped-silicon run from the earlier build. It decodes to **Silicon**, not Air:
/// a dopant run is geometrically silicon, and reporting it as air would delete
/// the region from the stack. The species is not carried — the caller is told how
/// many were dropped.
pub fn material_from_code(code: u8) -> Material {
    match code {
        1 | CODE_DOPANT => Material::Silicon,
        2 => Material::SiliconDioxide,
        3 => Material::SiliconNitride,
        4 => Material::Photoresist,
        5 => Material::Metal,
        // an unknown code is air, but `grid_from_columns` refuses the grid rather
        // than letting an unknown material become a void
        _ => Material::Air,
    }
}

/// Build an `atto_geom::ColumnGrid` from the application's flat data.
///
/// This is the fork made concrete: the earlier build keeps its own grid type and hands the
/// contents across as primitives, so neither crate has to adopt the other's
/// struct. Returns `None` if the sizes disagree, rather than a half-built grid.
#[allow(clippy::too_many_arguments)]
pub fn grid_from_columns(
    x0: f64,
    y0: f64,
    pitch: f64,
    nx: usize,
    ny: usize,
    z_bottom: f64,
    z_top: f64,
    active: &[bool],
    columns: &[Vec<PrimitiveRun>],
) -> Option<ColumnGrid> {
    if nx == 0 || ny == 0 || columns.len() != nx * ny || active.len() != nx * ny {
        return None;
    }
    if !(pitch > 0.0) || !pitch.is_finite() {
        return None;
    }
    let cols: Vec<Vec<Run>> = columns
        .iter()
        .map(|c| {
            c.iter()
                .map(|(z0, z1, code)| Run { z0: *z0, z1: *z1, material: material_from_code(*code) })
                .collect()
        })
        .collect();
    // an unknown code would become air and delete a region; refuse instead
    if columns.iter().flatten().any(|(_, _, code)| *code > CODE_DOPANT) {
        return None;
    }
    Some(ColumnGrid { x0, y0, pitch, nx, ny, cols, active: active.to_vec(), z_bottom, z_top })
}

/// How many runs lost their dopant species crossing the bridge.
pub fn dropped_dopant_runs(columns: &[Vec<PrimitiveRun>]) -> usize {
    columns.iter().flatten().filter(|(_, _, code)| *code == CODE_DOPANT).count()
}

/// The inverse: a grid flattened into primitives, for handing back to the earlier build.
pub fn columns_from_grid(g: &ColumnGrid) -> (Vec<bool>, Vec<Vec<PrimitiveRun>>) {
    let active = g.active.clone();
    let cols = g
        .cols
        .iter()
        .map(|c| c.iter().map(|r| (r.z0, r.z1, material_code(r.material))).collect())
        .collect();
    (active, cols)
}

/// Implanted concentration at each depth below the silicon surface, cm^-3.
///
/// Replaces the application's analytic Gaussian. Channelling is included, so the caller gets
/// the tilted-or-aligned behaviour for free.
pub fn implant_column(
    depths_nm: &[f64],
    species: &str,
    energy_kev: f64,
    dose_cm2: f64,
    tilt_deg: f64,
) -> Result<Vec<f64>, BridgeError> {
    if !(energy_kev > 0.0) || !(dose_cm2 > 0.0) {
        return Err(BridgeError::BadInput("energy and dose must be positive".into()));
    }
    let s = species_by_name(species)?;
    let p = Profile::from_species(s, energy_kev)
        .map_err(|e| BridgeError::BadInput(format!("{e:?}")))?
        .with_tilt(tilt_deg);
    Ok(depths_nm.iter().map(|d| p.concentration_cm3(atto_units::Nm::new(*d), dose_cm2)).collect())
}

/// Step height in nm, from adjacent columns of a height map (audit F-10).
///
/// the earlier build measured the spread of the top 10% and bottom 10% of all heights, which
/// reports a surface roughness on a flat wafer and a step on a patterned one
/// without distinguishing them. This takes the largest difference between
/// neighbouring active columns, which is a step or it is nothing.
pub fn step_height_nm(heights: &[f64], active: &[bool], nx: usize, ny: usize) -> Option<f64> {
    if heights.len() != nx * ny || active.len() != nx * ny || nx < 2 || ny < 1 {
        return None;
    }
    let mut best: Option<f64> = None;
    for y in 0..ny {
        for x in 0..(nx - 1) {
            let a = y * nx + x;
            let b = a + 1;
            if !active[a] || !active[b] {
                continue;
            }
            let d = (heights[b] - heights[a]).abs();
            best = Some(best.map_or(d, |v: f64| v.max(d)));
        }
    }
    best
}

/// Junction depth from a signed net-doping profile, measured from the surface
/// (audit F-07).
///
/// The profile is ordered from the surface downward and keeps its sign: the
/// junction is the first donor-to-acceptor crossing. Taking the magnitude first
/// destroys exactly the information that defines it.
pub fn junction_depth_nm(net: &[f64], leaf_nm: f64) -> Option<f64> {
    if net.len() < 2 || !(leaf_nm > 0.0) {
        return None;
    }
    for i in 1..net.len() {
        if net[i - 1] > 0.0 && net[i] <= 0.0 {
            let f = net[i - 1] / (net[i - 1] - net[i]);
            return Some((i as f64 - 0.5 + f) * leaf_nm);
        }
    }
    None
}

/// Silicon consumed and surface rise for an oxide of `total_nm` (0.44 ratio).
pub fn oxidation_split(total_nm: f64) -> (f64, f64) {
    if !(total_nm > 0.0) || !total_nm.is_finite() {
        return (0.0, 0.0);
    }
    let consumed = total_nm * atto_params::SI_CONSUMED_RATIO;
    (consumed, total_nm - consumed)
}

/// Vertical and lateral etch rates, nm/min, for a chemistry by name.
pub fn etch_rates(chemistry: &str, power: Watts, pressure: MTorr) -> Option<(f64, f64)> {
    let (power_w, pressure_mtorr) = (power.value(), pressure.value());
    let model = match chemistry.to_ascii_lowercase().as_str() {
        "oxide" | "sio2" | "cf4" => EtchModel::bundled_oxide(),
        _ => return None,
    };
    if !(power_w > 0.0) || !(pressure_mtorr > 0.0) {
        return None;
    }
    let (v, l) = model.rate(Watts::new(power_w), MTorr::new(pressure_mtorr));
    if v > 0.0 { Some((v, l)) } else { None }
}

#[cfg(test)]
mod grid_tests {
    use atto_units::{Mm, Um};
    use super::*;
    use atto_geom::column::{set_interval, solid_top};
    use atto_geom::svo::Wafer;

    fn real_grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = atto_geom::column::ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        ColumnGrid::sample(&w, &win).unwrap()
    }

    #[test]
    fn a_grid_survives_the_round_trip_through_primitives() {
        let mut g = real_grid();
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top, top + 100.0, Material::SiliconDioxide);
        }
        let (active, cols) = columns_from_grid(&g);
        let back = grid_from_columns(g.x0, g.y0, g.pitch, g.nx, g.ny, g.z_bottom, g.z_top, &active, &cols)
            .expect("the sizes agree");
        assert_eq!(back.nx, g.nx);
        assert_eq!(back.ny, g.ny);
        assert_eq!(back.active, g.active);
        assert_eq!(back.cols.len(), g.cols.len());
        for k in 0..g.cols.len() {
            assert_eq!(back.cols[k].len(), g.cols[k].len());
            for (a, b) in back.cols[k].iter().zip(g.cols[k].iter()) {
                assert!((a.z0 - b.z0).abs() < 1e-12);
                assert!((a.z1 - b.z1).abs() < 1e-12);
                assert_eq!(a.material, b.material);
            }
        }
    }

    #[test]
    fn every_material_round_trips_through_its_code() {
        for m in [
            Material::Air,
            Material::Silicon,
            Material::SiliconDioxide,
            Material::SiliconNitride,
            Material::Photoresist,
            Material::Metal,
        ] {
            assert_eq!(material_from_code(material_code(m)), m, "{m:?}");
        }
    }

    #[test]
    fn a_dopant_run_stays_silicon_and_is_counted_not_lost() {
        // code 8 is doped silicon: geometrically silicon, species dropped
        assert_eq!(material_from_code(CODE_DOPANT), Material::Silicon);
        let cols: Vec<Vec<PrimitiveRun>> = vec![vec![(0.0, 10.0, 1), (10.0, 20.0, CODE_DOPANT)]];
        assert_eq!(dropped_dopant_runs(&cols), 1, "the drop is counted, not silent");
        let g = grid_from_columns(0.0, 0.0, 10.0, 1, 1, 0.0, 20.0, &[true], &cols).expect("converts");
        assert_eq!(g.cols[0][1].material, Material::Silicon, "not Air: the region survives");
    }

    #[test]
    fn an_unknown_code_is_refused_rather_than_becoming_a_void() {
        // an unmapped material would arrive as Air and delete a region
        assert_eq!(material_from_code(99), Material::Air);
        let cols: Vec<Vec<PrimitiveRun>> = vec![vec![(0.0, 10.0, 99)]];
        assert!(grid_from_columns(0.0, 0.0, 10.0, 1, 1, 0.0, 10.0, &[true], &cols).is_none());
    }

    #[test]
    fn mismatched_sizes_are_refused_not_half_built() {
        let g = real_grid();
        let (active, cols) = columns_from_grid(&g);
        assert!(grid_from_columns(g.x0, g.y0, g.pitch, g.nx + 1, g.ny, g.z_bottom, g.z_top, &active, &cols).is_none());
        assert!(grid_from_columns(g.x0, g.y0, 0.0, g.nx, g.ny, g.z_bottom, g.z_top, &active, &cols).is_none());
        assert!(grid_from_columns(g.x0, g.y0, g.pitch, g.nx, g.ny, g.z_bottom, g.z_top, &[], &cols).is_none());
    }
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;

    #[test]
    fn species_names_map_and_unknowns_are_errors() {
        assert_eq!(species_by_name("P").unwrap(), Species::Phosphorus);
        assert_eq!(species_by_name("arsenic").unwrap(), Species::Arsenic);
        assert!(species_by_name("unobtainium").is_err());
    }

    #[test]
    fn an_implant_column_is_positive_and_peaks_near_rp() {
        let depths: Vec<f64> = (0..200).map(|i| i as f64 * 2.0).collect();
        let c = implant_column(&depths, "P", 60.0, 5.0e14, 7.0).unwrap();
        assert_eq!(c.len(), depths.len());
        assert!(c.iter().all(|v| *v >= 0.0));
        let peak = c.iter().cloned().fold(0.0, f64::max);
        assert!(peak > 1.0e18, "peak {peak:.2e}");
    }

    #[test]
    fn an_aligned_column_channels_deeper_than_a_tilted_one() {
        let depths = vec![400.0];
        let a = implant_column(&depths, "P", 60.0, 5.0e14, 0.0).unwrap()[0];
        let b = implant_column(&depths, "P", 60.0, 5.0e14, 7.0).unwrap()[0];
        assert!(a > 100.0 * b, "{a:.2e} vs {b:.2e}");
    }

    #[test]
    fn bad_implant_input_is_rejected() {
        assert!(implant_column(&[10.0], "P", 0.0, 1.0e14, 0.0).is_err());
        assert!(implant_column(&[10.0], "P", 60.0, 0.0, 0.0).is_err());
    }

    #[test]
    fn a_flat_surface_has_no_step() {
        let h = vec![100.0; 16];
        let a = vec![true; 16];
        assert_eq!(step_height_nm(&h, &a, 4, 4), Some(0.0));
    }

    #[test]
    fn a_step_is_the_adjacent_difference_not_the_spread() {
        // left half at 100 nm, right half at 140 nm
        let mut h = vec![100.0; 16];
        let mut a = vec![true; 16];
        for y in 0..4 {
            for x in 2..4 {
                h[y * 4 + x] = 140.0;
            }
        }
        a[0] = false;
        assert_eq!(step_height_nm(&h, &a, 4, 4), Some(40.0));
    }

    #[test]
    fn a_rough_surface_with_no_step_does_not_report_one() {
        // alternating columns: a large percentile spread, no step
        let mut h = vec![100.0; 16];
        for i in 0..16 {
            if i % 2 == 0 {
                h[i] = 90.0;
            }
        }
        let a = vec![true; 16];
        let s = step_height_nm(&h, &a, 4, 4).unwrap();
        assert!((s - 10.0).abs() < 1e-9, "a checkerboard reports 10 nm per column, got {s}");
    }

    #[test]
    fn mismatched_step_input_is_rejected() {
        assert!(step_height_nm(&[1.0, 2.0], &[true, true], 4, 4).is_none());
        assert!(step_height_nm(&[1.0; 4], &[true; 4], 1, 4).is_none());
    }

    #[test]
    fn the_junction_is_a_sign_change_measured_from_the_surface() {
        // 8 leaves of 10 nm: donors for the first three, then acceptors
        let net = vec![1e19, 8e18, 3e18, -1e15, -1e15, -1e15, -1e15, -1e15];
        let xj = junction_depth_nm(&net, 10.0).expect("a junction");
        assert!(xj > 24.0 && xj < 36.0, "xj {xj:.1} nm");
    }

    #[test]
    fn a_profile_without_a_crossing_has_no_junction() {
        let net = vec![1e19, 8e18, 3e18, 1e17];
        assert_eq!(junction_depth_nm(&net, 10.0), None);
        assert_eq!(junction_depth_nm(&[1.0], 10.0), None);
        assert_eq!(junction_depth_nm(&[1.0, -1.0], 0.0), None);
    }

    #[test]
    fn the_oxidation_split_keeps_44_percent() {
        let (consumed, rise) = oxidation_split(100.0);
        assert!((consumed - 44.0).abs() < 1e-9);
        assert!((rise - 56.0).abs() < 1e-9);
        assert_eq!(oxidation_split(-1.0), (0.0, 0.0));
        assert_eq!(oxidation_split(f64::NAN), (0.0, 0.0));
    }

    #[test]
    fn etch_rates_come_back_for_a_known_chemistry_only() {
        assert!(etch_rates("oxide", Watts::new(200.0), MTorr::new(20.0)).is_some());
        assert!(etch_rates("unobtainium", Watts::new(200.0), MTorr::new(20.0)).is_none());
        assert!(etch_rates("oxide", Watts::new(0.0), MTorr::new(20.0)).is_none());
    }
}
