// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Verification suite for AttoFab's own process operators and FEM.
//!
//! *Verification* asks "are we solving the equations right?" It needs no
//! experiment: it compares against exact solutions, checks conservation, and studies
//! discretisation error. *Validation* ("are these the right equations for a real
//! fab?") needs measured data and is handled separately in [`crate::signoff`].
//!
//! The SVO lattice quantises every surface to the octree cell, so operator errors are
//! sawtooth, not asymptotic: classical GCI does not apply to them. Those cases are
//! therefore judged against an explicit **cell-size bound** (error ≤ k · cell) and the
//! convergence table is reported, while GCI is used where the response really is smooth
//! (the FEM mesh study).

use crate::fem::benchmark::{strip_curvature, timoshenko};
use crate::geom::{self, ColumnGrid, EtchChemistry, EtchSpec, MaskPattern, ProcessWindow};
use crate::params::{DepositionParams, DepositionSource, DopingParams, LithoParams, StageParams};
use crate::physics::Stage;
use crate::pipeline::evaluate_stage;
use crate::validation::gci::{calculate_gci, GridResult};
use crate::validation::{CheckStatus, Finding};
use crate::wafer::{Material, WaferState};
use serde::Serialize;

#[derive(Debug, Clone)]
pub struct VvConfig {
    /// Octree cell sizes, coarse → fine (tile edge is kept at 2048 nm).
    pub leaf_sizes_nm: Vec<f64>,
    pub fem_film_cells: [usize; 3],
}

impl Default for VvConfig { fn default() -> Self { Self { leaf_sizes_nm: vec![32.0, 16.0, 8.0], fem_film_cells: [1, 2, 3] } } }

#[derive(Debug, Clone, Serialize)]
pub struct VvCase {
    pub id: String,
    pub group: &'static str,
    pub status: CheckStatus,
    /// Observed error (or violation) and the limit it is judged against.
    pub value: f64,
    pub limit: f64,
    pub unit: &'static str,
    pub message: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct VvReport { pub cases: Vec<VvCase>, pub leaf_sizes_nm: Vec<f64> }

impl VvReport {
    pub fn failed(&self) -> Vec<&VvCase> { self.cases.iter().filter(|c| c.status == CheckStatus::Fail).collect() }
    pub fn all_pass(&self) -> bool { !self.cases.is_empty() && self.cases.iter().all(|c| c.status == CheckStatus::Pass) }
    pub fn findings(&self) -> Vec<Finding> { self.cases.iter().map(|c| Finding { id: c.id.clone(), status: c.status, message: c.message.clone() }).collect() }
}

fn win() -> ProcessWindow { ProcessWindow { center: [0.0, 0.0], size: [1200.0, 1200.0] } }

pub fn wafer_for_leaf(leaf_nm: f64) -> WaferState {
    let levels = (2048.0 / leaf_nm).log2().round().max(1.0) as u32;
    WaferState::new_wafer_tile(300.0, 500.0, leaf_nm, levels, (0.0, 0.0)).expect("valid wafer")
}

fn case(id: impl Into<String>, group: &'static str, value: f64, limit: f64, unit: &'static str, message: String) -> VvCase {
    let status = if value.is_finite() && value <= limit { CheckStatus::Pass } else { CheckStatus::Fail };
    VvCase { id: id.into(), group, status, value, limit, unit, message }
}

fn solid_volume(g: &ColumnGrid) -> f64 { geom::material_volume_nm3(g, |m| m != Material::Air) }
fn vol(g: &ColumnGrid, m: Material) -> f64 { geom::material_volume_nm3(g, |x| x == m) }

fn run(w: &WaferState, p: &StageParams) -> crate::pipeline::StageEval { evaluate_stage(w, &win(), p, 1, 1.0).expect("stage evaluates") }

fn dep(thickness: f64, step: f64) -> StageParams { StageParams::Deposition(DepositionParams { material: Material::SiliconDioxide, thickness_nm: thickness, step_coverage: step, bow_radius_m: 40.0, flux3d: false, sticking: 1.0, source: DepositionSource::Cosine, flux_iterations: 1, flux_tolerance: 0.0 }) }

/// Isolated dark line: retained-resist width from the exact Gaussian-blurred aerial image.
fn analytic_line_cd(width: f64, sigma: f64, dose_ratio: f64) -> f64 {
    let erf = |x: f64| { // Abramowitz-Stegun 7.1.26, |err| < 1.5e-7
        let t = 1.0 / (1.0 + 0.3275911 * x.abs());
        let y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * (-x * x).exp();
        if x >= 0.0 { y } else { -y }
    };
    let bright = |x: f64| 0.5 * (erf((x + width / 2.0) / (std::f64::consts::SQRT_2 * sigma)) - erf((x - width / 2.0) / (std::f64::consts::SQRT_2 * sigma)));
    let thr = 1.0 - 1.0 / dose_ratio; // retained where dark-line intensity 1 - bright < 1/dose  <=>  bright > 1 - 1/dose
    let (mut lo, mut hi) = (0.0, width * 2.0);
    for _ in 0..80 { let mid = 0.5 * (lo + hi); if bright(mid) > thr { lo = mid } else { hi = mid } }
    2.0 * lo
}

/// Isolated dark line after development: the same 1-D problem with the developer included.
///
/// The retained width is where the clear-field intensity ratio at which the developer just
/// reaches the substrate (found by bisection on the Dill/Mack depth) crosses the Gaussian line
/// image. A pure threshold at 1/dose ignores that the developer needs a finite dose to clear
/// the film within the develop time, so it is the wrong reference for a developed profile. A
/// column counts as cleared once less than half a cell of resist is left, because the grid
/// stores heights in whole cells.
fn analytic_developed_line_cd(width: f64, sigma: f64, dose_ratio: f64, thickness_nm: f64, time_min: f64, leaf_nm: f64) -> f64 {
    let mack = atto_process::resist::DevelopModel::default();
    let dill = atto_process::dill::DillModel::default();
    let clears = |i: f64| atto_process::dill::develop_depth(i * dose_ratio, atto_units::Nm::new(thickness_nm), &mack, &dill, atto_units::Minutes::new(time_min), 16).0 >= thickness_nm - 0.5 * leaf_nm;
    // lowest intensity that clears the film (retained where intensity is below it)
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    if !clears(hi) { return f64::NAN; }
    for _ in 0..60 { let mid = 0.5 * (lo + hi); if clears(mid) { hi = mid } else { lo = mid } }
    let i_c = hi;
    let erf = |x: f64| {
        let t = 1.0 / (1.0 + 0.3275911 * x.abs());
        let y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * (-x * x).exp();
        if x >= 0.0 { y } else { -y }
    };
    let dark = |x: f64| 1.0 - 0.5 * (erf((x + width / 2.0) / (std::f64::consts::SQRT_2 * sigma)) - erf((x - width / 2.0) / (std::f64::consts::SQRT_2 * sigma)));
    let (mut a, mut b) = (0.0, width * 2.0);
    for _ in 0..80 { let mid = 0.5 * (a + b); if dark(mid) < i_c { a = mid } else { b = mid } }
    2.0 * a
}

pub fn run_suite(cfg: &VvConfig) -> VvReport {
    let mut cases: Vec<VvCase> = Vec::new();
    // (id, per-leaf (pitch, error)) for the convergence table.
    for &leaf in &cfg.leaf_sizes_nm {
        let w0 = wafer_for_leaf(leaf);
        let g0 = ColumnGrid::sample(&w0, &win()).expect("grid");
        let pitch = g0.pitch;
        let area = g0.active_count() as f64 * g0.cell_area();
        let tag = |s: &str| format!("{s}@{leaf:.0}nm");

        // ---- CONSERVATION ------------------------------------------------------------
        let d = run(&w0, &dep(100.0, 0.0));
        let g1 = ColumnGrid::sample(&d.wafer, &win()).unwrap();
        let realised = d.outcome.metric("Realised mean thickness").unwrap();
        let v_film = vol(&g1, Material::SiliconDioxide);
        cases.push(case(tag("CONS-DEP-svo"), "conservation", (v_film - realised * area).abs() / (realised * area), 1e-9, "rel", format!("film volume in the SVO ({v_film:.4e} nm³) equals realised thickness × window area: the octree write-back conserves volume exactly")));
        cases.push(case(tag("CONS-DEP-quantisation"), "conservation", (realised - 100.0).abs(), pitch * 0.5 + 1e-9, "nm", format!("realised {realised:.1} nm vs requested 100 nm; the lattice allows at most half a cell ({:.1} nm)", pitch * 0.5)));

        let litho = run(&d.wafer, &StageParams::default_for(Stage::Lithography));
        let etched = run(&litho.wafer, &StageParams::default_for(Stage::Etching));
        let (ga, gb) = (ColumnGrid::sample(&litho.wafer, &win()).unwrap(), ColumnGrid::sample(&etched.wafer, &win()).unwrap());
        let removed_solid = solid_volume(&ga) - solid_volume(&gb);
        let spec = if let StageParams::Etching(e) = StageParams::default_for(Stage::Etching) { EtchSpec { rates: e.rates, anisotropy: e.anisotropy, time_min: e.time_min } } else { unreachable!() };
        let (_, st) = geom::etch(&ga, &spec, 1.0).unwrap();
        let changed = ga.cols.iter().zip(&gb.cols).filter(|(a, b)| a != b).count() as f64;
        let bound = 0.5 * pitch * g0.cell_area() * changed.max(1.0);
        let acct = st.vertical_removed_nm3 + st.lateral_removed_nm3;
        cases.push(case(tag("CONS-ETCH"), "conservation", (removed_solid - acct).abs(), bound, "nm³", format!("solid removed from the SVO {removed_solid:.4e} vs accounted by the operator {acct:.4e}; difference bounded by lattice snapping over {changed:.0} changed columns")));
        cases.push(case(tag("CONS-ETCH-nonneg"), "conservation", if removed_solid >= 0.0 { 0.0 } else { 1.0 }, 0.0, "", "etching never adds material".into()));

        let dop = run(&etched.wafer, &StageParams::default_for(Stage::Doping));
        let gc = ColumnGrid::sample(&dop.wafer, &win()).unwrap();
        cases.push(case(tag("CONS-DOPING-tags"), "conservation", (solid_volume(&gc) - solid_volume(&gb)).abs() / solid_volume(&gb), 1e-12, "rel", "implant re-labels material: total solid volume is unchanged (no growth is implied)".into()));

        let clean = run(&dop.wafer, &StageParams::default_for(Stage::FinalClean));
        let gd = ColumnGrid::sample(&clean.wafer, &win()).unwrap();
        cases.push(case(tag("CONS-STRIP"), "conservation", vol(&gd, Material::Photoresist), 0.0, "nm³", "a 100% strip leaves no resist".into()));
        let other = |g: &ColumnGrid| solid_volume(g) - vol(g, Material::Photoresist);
        cases.push(case(tag("CONS-STRIP-only-resist"), "conservation", (other(&gd) - other(&gc)).abs() / other(&gc), 1e-12, "rel", "stripping removes resist and nothing else".into()));

        // ---- EXACT SOLUTIONS -----------------------------------------------------------
        // Blanket etch depth = rate × time.
        let blanket = run(&w0, &dep(100.0, 0.0));
        let gb0 = ColumnGrid::sample(&blanket.wafer, &win()).unwrap();
        let (rate, t_min) = (60.0, 0.5);
        let blanket_spec = EtchSpec { rates: EtchChemistry::OxideRie.default_rates(), anisotropy: 1.0, time_min: t_min };
        let (eb, _) = geom::etch(&gb0, &blanket_spec, 1.0).unwrap();
        let (h0, h1) = (gb0.surface_heights()[gb0.idx(gb0.nx / 2, gb0.ny / 2)], eb.surface_heights()[eb.idx(eb.nx / 2, eb.ny / 2)]);
        cases.push(case(tag("EXACT-ETCH-DEPTH"), "exact-solution", ((h0 - h1) - rate * t_min).abs(), pitch * 0.5 + 1e-9, "nm", format!("blanket oxide etch removed {:.1} nm; exact = rate × time = {:.1} nm", h0 - h1, rate * t_min)));

        // Lateral undercut = (1 − A) · rate · T at the exposed top.
        let (a_iso, t_iso) = (0.5, 1.0);
        let stack = ColumnGrid::sample(&litho.wafer, &win()).unwrap();
        let iso = EtchSpec { rates: EtchChemistry::OxideRie.default_rates(), anisotropy: a_iso, time_min: t_iso };
        let (_, sti) = geom::etch(&stack, &iso, 1.0).unwrap();
        let exact_u = (1.0 - a_iso) * rate * t_iso;
        cases.push(case(tag("EXACT-UNDERCUT"), "exact-solution", (sti.max_lateral_nm - exact_u).abs(), pitch * 1.5 + 1e-9, "nm", format!("realised lateral reach {:.1} nm vs exact (1−A)·rate·T = {:.1} nm", sti.max_lateral_nm, exact_u)));

        // Isolated resist line CD from the analytic aerial image.
        let mut lp = if let StageParams::Lithography(l) = StageParams::default_for(Stage::Lithography) { l } else { unreachable!() };
        lp.mask = MaskPattern::LineSpace { line_nm: 203.0, pitch_nm: 400.0, count: 1, length_nm: 1000.0, vertical: true, clear_field: true };
        // A develop time well past the point where the clear field opens, so the contour sits on the
        // developer's saturated plateau and not on its steep edge, where one cell of rounding moves it by tens of nm.
        lp.develop_time_min = 0.7;
        let lp: LithoParams = lp;
        let sigma = lp.k1 * lp.wavelength_nm / lp.na / 2.355;
        let lres = run(&blanket.wafer, &StageParams::Lithography(lp.clone()));
        let cd = lres.outcome.metric("Resist CD (mean line width)").unwrap_or(f64::NAN);
        let sigma_bake = atto_process::dill::bake_sigma_from_temperature(atto_units::Kelvin::new(lp.peb_c + 273.15), atto_units::Minutes::new(1.0)).0;
        let exact_cd = analytic_developed_line_cd(203.0, (sigma * sigma + sigma_bake * sigma_bake).sqrt(), lp.dose_ratio, lp.resist_thickness_nm, lp.develop_time_min, leaf);
        cases.push(case(tag("EXACT-RESIST-CD"), "exact-solution", (cd - exact_cd).abs(), 1.5 * pitch + 1e-9, "nm", format!("resist CD {cd:.1} nm vs exact Gaussian-image, baked, developed CD {exact_cd:.1} nm")));

        // Implant: thickness of the ≥1e17 cm⁻³ region vs the analytic Gaussian contour.
        let dp = DopingParams { species: crate::params::Species::Boron, energy_kev: 60.0, dose_cm2: 1.0e14, anneal_temp_c: 900.0, anneal_min: 10.0, background_cm3: 1.0e15, rp_override_nm: None, tilt_deg: 7.0 };
        let dres = run(&w0, &StageParams::Doping(dp));
        let (rp, sg, peak) = (dres.outcome.metric("Projected range Rp").unwrap(), dres.outcome.metric("Straggle after anneal").unwrap(), dres.outcome.metric("Peak concentration").unwrap());
        let half = sg * (2.0 * (peak / 1.0e17).ln()).sqrt();
        let (d1, d2) = ((rp - half).max(0.0), rp + half);
        let col = dres.wafer.column_runs(0.0, 0.0);
        let measured: f64 = col.iter().filter(|r| matches!(r.material, Material::Dopant(_))).map(|r| r.z1 - r.z0).sum();
        cases.push(case(tag("EXACT-DOPANT-CONTOUR"), "exact-solution", (measured - (d2 - d1)).abs(), pitch * 2.0 + 1e-9, "nm", format!("1e17 contour thickness {measured:.1} nm vs analytic {:.1} nm", d2 - d1)));
    }

    // ---- CONVERGENCE TABLE (error must not grow as cells shrink) ----------------------------
    for key in ["EXACT-ETCH-DEPTH", "EXACT-UNDERCUT", "EXACT-RESIST-CD", "EXACT-DOPANT-CONTOUR", "CONS-DEP-quantisation"] {
        let errs: Vec<(f64, f64)> = cfg.leaf_sizes_nm.iter().filter_map(|l| cases.iter().find(|c| c.id == format!("{key}@{l:.0}nm")).map(|c| (*l, c.value))).collect();
        if errs.len() >= 3 {
            let worst_ratio = errs.iter().map(|(l, e)| e / l).fold(0.0, f64::max);
            let gci = calculate_gci(GridResult { spacing: errs[0].0, value: errs[0].1 }, GridResult { spacing: errs[1].0, value: errs[1].1 }, GridResult { spacing: errs[2].0, value: errs[2].1 }, 1.25);
            let note = match gci { Ok(g) if g.monotonic_convergence => format!("monotone (observed order {:.2})", g.observed_order), _ => "non-asymptotic (lattice-quantised): judged by the cell bound, not GCI".into() };
            cases.push(case(format!("CONV-{key}"), "convergence", worst_ratio, 2.0, "error/cell", format!("errors by cell {errs:?}: worst error is {worst_ratio:.2} cells; {note}")));
        }
    }

    // ---- FEM -------------------------------------------------------------------------------
    let (hs, hf, l, sigma) = (100.0, 10.0, 1200.0, 200.0e6);
    let exact = timoshenko(hs, hf, sigma);
    let ks: Vec<f64> = cfg.fem_film_cells.iter().map(|n| strip_curvature(hs, hf, l, *n, sigma)).collect();
    let errs: Vec<f64> = ks.iter().map(|k| (k - exact).abs() / exact).collect();
    cases.push(case("FEM-BEAM", "fem", errs[2], 0.05, "rel", format!("film-on-substrate curvature: FEM {:.4e} vs exact two-layer beam {:.4e}; error by mesh {errs:.3?}", ks[2], exact)));
    cases.push(case("FEM-REFINEMENT", "fem", if errs[2] <= errs[0] { 0.0 } else { 1.0 }, 0.0, "", "refining the mesh does not worsen the error".into()));
    let spacing = |i: usize| hf / cfg.fem_film_cells[i] as f64;
    if let Ok(g) = calculate_gci(GridResult { spacing: spacing(0), value: ks[0] }, GridResult { spacing: spacing(1), value: ks[1] }, GridResult { spacing: spacing(2), value: ks[2] }, 1.25) {
        cases.push(case("FEM-GCI", "fem", g.fine_gci_percent, 6.0, "%", format!("fine-mesh GCI {:.2}% (observed order {:.2})", g.fine_gci_percent, g.observed_order)));
    }

    VvReport { cases, leaf_sizes_nm: cfg.leaf_sizes_nm.clone() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn suite_passes_on_the_coarse_and_default_cells() {
        let r = run_suite(&VvConfig { leaf_sizes_nm: vec![32.0, 16.0], fem_film_cells: [1, 2, 3] });
        for c in &r.cases { eprintln!("{:34} {:?} value {:.4} limit {:.4} {}", c.id, c.status, c.value, c.limit, c.unit); }
        assert!(r.failed().is_empty(), "failed: {:#?}", r.failed());
        assert!(r.cases.len() >= 20);
    }

    #[test]
    fn analytic_cd_is_monotone_in_dose_and_limits_are_sane() {
        let a = analytic_line_cd(200.0, 38.0, 1.5);
        let b = analytic_line_cd(200.0, 38.0, 3.0);
        assert!(b < a, "positive resist: more dose dissolves more, so a dark line gets narrower: {a} vs {b}");
        assert!((analytic_line_cd(200.0, 1.0, 2.0) - 200.0).abs() < 1.0, "no blur: CD = drawn width at 50% threshold");
    }

}
