// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Stage evaluation (PRD §2, §3, §7 P3).
//!
//! [`evaluate_stage`] is the *only* place a fabrication stage is computed. It
//! takes the wafer state before the stage, the process window, the stage
//! parameters, a seed and a **progress** value in `[0, 1]`:
//!
//! * `progress = 1.0` is the committed result,
//! * `progress < 1.0` is an animation frame, produced by the same operators at
//!   a fractional parameter (film thickness × f, etch time × f, …).
//!
//! Because both paths are one code path, the last animation frame is the
//! committed geometry by construction, and every frame is deterministically
//! reproducible from `(wafer_before, window, params, seed, progress)`.

use crate::geom::{self, ColumnGrid, GeomError, ProcessWindow, EtchSpec};
use crate::measure::{self, MetrologyReport, Status};
use crate::params::{DopingParams, Issue, Severity, StageParams};
use crate::physics::{deposition, etching, lithography, prep, Stage};
use crate::wafer::{DopantProfile, FilmStress, Material, Run, WaferState};
use serde::{Deserialize, Serialize};

/// Animation keyframes per stage (frames `0..=FRAMES_PER_STAGE`).
pub const FRAMES_PER_STAGE: usize = 24;

pub fn frame_progress(frame: usize) -> f64 { (frame.min(FRAMES_PER_STAGE)) as f64 / FRAMES_PER_STAGE as f64 }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Fidelity {
    /// Computed directly from the simulated geometry/analytic model.
    Simulated,
    /// Physically motivated reference model; not calibrated to a tool.
    Approximated,
    /// Calibrated against measured tool data.
    Calibrated,
}

impl Fidelity {
    pub fn label(self) -> &'static str {
        match self { Fidelity::Simulated => "Simulated", Fidelity::Approximated => "Approximated", Fidelity::Calibrated => "Calibrated" }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StepKind {
    /// Changes the SVO geometry (material is added or removed).
    Geometry,
    /// Re-labels existing material (dopant concentration buckets). Never drawn as growth.
    MaterialTag,
    /// Closed-form calculation / measurement; no geometry change.
    Analytic,
    /// Recorded in the process log only; not modelled.
    Recorded,
}

impl StepKind {
    pub fn label(self) -> &'static str {
        match self { StepKind::Geometry => "geometry", StepKind::MaterialTag => "material tags", StepKind::Analytic => "analytic", StepKind::Recorded => "recorded only" }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub struct SubStep { pub name: &'static str, pub kind: StepKind, pub t0: f64, pub t1: f64 }

const fn ss(name: &'static str, kind: StepKind, t0: f64, t1: f64) -> SubStep { SubStep { name, kind, t0, t1 } }

pub fn substeps(stage: Stage) -> Vec<SubStep> {
    use StepKind::*;
    match stage {
        Stage::SubstratePrep => vec![ss("Wafer definition", Recorded, 0.0, 0.15), ss("Initial substrate", Recorded, 0.15, 0.3), ss("Surface preparation (RCA / HF)", Analytic, 0.3, 0.6), ss("Clean · rinse · dry", Analytic, 0.6, 0.85), ss("Substrate inspection", Analytic, 0.85, 1.0)],
        Stage::Deposition => vec![ss("Deposition recipe", Recorded, 0.0, 0.1), ss("Film growth", Geometry, 0.1, 0.8), ss("Thickness / uniformity preview", Analytic, 0.8, 0.9), ss("Stress estimate", Analytic, 0.9, 0.95), ss("Film inspection", Analytic, 0.95, 1.0)],
        Stage::Lithography => vec![ss("Dehydration bake", Recorded, 0.0, 0.07), ss("Adhesion promotion", Recorded, 0.07, 0.15), ss("Photoresist coat", Geometry, 0.15, 0.35), ss("Soft bake", Recorded, 0.35, 0.45), ss("Mask alignment", Analytic, 0.45, 0.5), ss("Exposure", Analytic, 0.5, 0.65), ss("Post-exposure bake", Recorded, 0.65, 0.75), ss("Develop", Geometry, 0.75, 0.95), ss("Resist inspection", Analytic, 0.95, 1.0)],
        Stage::Etching => vec![ss("Etch chemistry selection", Recorded, 0.0, 0.1), ss("Selectivity & anisotropy setup", Analytic, 0.1, 0.25), ss("Etch execution", Geometry, 0.25, 0.8), ss("Etch inspection", Analytic, 0.8, 0.85), ss("Resist strip (optional)", Geometry, 0.85, 1.0)],
        Stage::Doping => vec![ss("Implant parameters", Recorded, 0.0, 0.1), ss("Ion implantation", MaterialTag, 0.1, 0.5), ss("Thermal activation / anneal", MaterialTag, 0.55, 0.9), ss("Junction / profile inspection", Analytic, 0.9, 1.0)],
        Stage::FinalClean => vec![ss("Resist removal", Geometry, 0.0, 0.45), ss("Chemical clean", Analytic, 0.45, 0.65), ss("Rinse", Recorded, 0.65, 0.8), ss("Dry", Recorded, 0.8, 0.9), ss("Residue / contamination inspection", Analytic, 0.9, 1.0)],
        Stage::Cmp => vec![ss("Pad conditioning", Recorded, 0.0, 0.15), ss("Polish", Geometry, 0.15, 0.75), ss("Over-polish / stop layer", Geometry, 0.75, 0.9), ss("Post-CMP inspection", Analytic, 0.9, 1.0)],
        Stage::Metrology => ["Thickness", "Topography", "Critical dimensions", "Etch depth", "Dopant profile", "Stress", "Optical / material", "Final report"].iter().enumerate()
            .map(|(i, n)| ss(n, Analytic, i as f64 / 8.0, (i + 1) as f64 / 8.0)).collect(),
    }
}

pub fn phase_at(stage: Stage, progress: f64) -> SubStep {
    let steps = substeps(stage);
    steps.iter().find(|s| progress >= s.t0 && progress < s.t1).copied().unwrap_or(*steps.last().unwrap())
}

fn g(p: f64, a: f64, b: f64) -> f64 { ((p - a) / (b - a)).clamp(0.0, 1.0) }

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Metric { pub name: String, pub value: f64, pub unit: &'static str }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ExposureOverlay { pub x0: f64, pub y0: f64, pub pitch: f64, pub nx: usize, pub ny: usize, pub intensity: Vec<f32>, pub removed: Vec<bool> }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Site { pub x: f64, pub y: f64, pub value: Option<f64>, pub out_of_spec: bool }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum Overlay {
    None,
    Surface { clean_fraction: f64, particles_cm2: f64 },
    Thickness { target_nm: f64, mean_nm: f64, min_nm: f64, max_nm: f64 },
    Exposure(ExposureOverlay),
    Sites(Vec<Site>),
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct StageOutcome {
    pub stage: Stage,
    pub fidelity: Fidelity,
    pub substeps: Vec<SubStep>,
    pub metrics: Vec<Metric>,
    pub issues: Vec<Issue>,
    pub geometry_changed: bool,
    pub duration_s: f64,
    pub report: Option<MetrologyReport>,
}

impl StageOutcome {
    fn new(stage: Stage, params: &StageParams) -> Self {
        Self { stage, fidelity: if stage == Stage::Metrology { Fidelity::Simulated } else { Fidelity::Approximated }, substeps: substeps(stage), metrics: vec![],
               issues: params.validate(), geometry_changed: false, duration_s: params.duration_s(), report: None }
    }
    fn m(&mut self, name: &str, value: f64, unit: &'static str) { self.metrics.push(Metric { name: name.into(), value, unit }); }
    fn warn(&mut self, msg: impl Into<String>) { self.issues.push(Issue::warn(msg)); }
    pub fn metric(&self, name: &str) -> Option<f64> { self.metrics.iter().find(|m| m.name == name).map(|m| m.value) }
    pub fn has_errors(&self) -> bool { self.issues.iter().any(|i| i.severity == Severity::Error) }
}

pub struct StageEval { pub wafer: WaferState, pub outcome: StageOutcome, pub overlay: Overlay }

#[derive(Debug, Clone, PartialEq, thiserror::Error, Serialize, Deserialize)]
pub enum PipelineError {
    #[error("invalid parameter: {0}")]
    Invalid(String),
    #[error("missing or unusable input: {0}")]
    MissingInput(String),
    #[error("geometry: {0}")]
    Geometry(String),
    #[error("physics model: {0}")]
    Physics(String),
}

impl From<GeomError> for PipelineError {
    fn from(e: GeomError) -> Self { match e { GeomError::InvalidParameter(m) => PipelineError::Invalid(m), other => PipelineError::Geometry(other.to_string()) } }
}

/// Evaluate one stage on `before` at `progress` (see module docs).
pub fn evaluate_stage(before: &WaferState, window: &ProcessWindow, params: &StageParams, seed: u64, progress: f64) -> Result<StageEval, PipelineError> {
    if let Some(e) = params.validate().into_iter().find(|i| i.severity == Severity::Error) { return Err(PipelineError::Invalid(e.message)); }
    if !progress.is_finite() { return Err(PipelineError::Invalid("progress must be finite".into())); }
    window.validate()?;
    let p = progress.clamp(0.0, 1.0);
    let stage = params.stage();
    let mut wafer = before.clone();
    let mut out = StageOutcome::new(stage, params);
    let overlay = match params {
        StageParams::SubstratePrep(x) => run_prep(&mut wafer, x, p, &mut out)?,
        StageParams::Deposition(x) => run_deposition(&mut wafer, window, x, p, &mut out)?,
        StageParams::Lithography(x) => run_litho(&mut wafer, window, x, p, &mut out)?,
        StageParams::Etching(x) => run_etch(&mut wafer, window, x, p, &mut out)?,
        StageParams::Doping(x) => run_doping(&mut wafer, window, x, p, &mut out)?,
        StageParams::FinalClean(x) => run_clean(&mut wafer, window, x, p, &mut out)?,
        StageParams::Metrology(x) => run_metrology(&mut wafer, window, x, seed, p, &mut out)?,
        StageParams::Cmp(x) => run_cmp(&mut wafer, window, x, p, &mut out)?,
    };
    out.geometry_changed = wafer.fingerprint() != before.fingerprint();
    Ok(StageEval { wafer, outcome: out, overlay })
}

fn run_prep(w: &mut WaferState, x: &crate::params::PrepParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let r = prep::run(&prep::PrepInputs { gamma_sg: x.gamma_sg, gamma_sl: x.gamma_sl, gamma_lg: x.gamma_lg, ra_nm: x.ra_nm, rca_dilution: x.rca_dilution })
        .map_err(|e| PipelineError::Invalid(e.to_string()))?;
    if !r.wettable { out.warn(format!("surface is not wettable (contact angle {:.1}°)", r.theta_deg)); }
    if !r.validation_passed { out.warn("substrate failed the pinhole / roughness gate"); }
    let before = w.meta.surface.particle_density_cm2;
    let after = before * (1.0 - x.particle_removal_eff);
    out.m("Contact angle", r.theta_deg, "°");
    out.m("Particle density before", before, "cm⁻²");
    out.m("Particle density after", after, "cm⁻²");
    out.m("Native oxide after", if x.hf_dip { 0.0 } else { w.meta.surface.native_oxide_nm }, "nm");
    let frac = g(p, 0.3, 0.85);
    if p >= 1.0 {
        w.meta.surface.particle_density_cm2 = after;
        w.meta.surface.contact_angle_deg = Some(r.theta_deg);
        if x.hf_dip { w.meta.surface.native_oxide_nm = 0.0; }
    }
    Ok(Overlay::Surface { clean_fraction: frac, particles_cm2: before * (1.0 - x.particle_removal_eff * frac) })
}


/// Measure the openings in the window and ask the fill model what would happen.
///
/// A metal deposition over a hole either fills it or seals it with a void. The fill model
/// knows which, given a width and a depth; this finds those from the grid: the deepest
/// columns are the holes, and the run of them through the narrowest direction is the width.
fn via_fill_check(grid: &ColumnGrid, metal_thickness_nm: f64, step_coverage: f64) -> Option<(f64, f64, atto_process::via::ViaFill)> {
    if grid.cols.is_empty() || !(grid.pitch > 0.0) {
        return None;
    }
    let h = grid.surface_heights();
    let active: Vec<f64> = (0..h.len()).filter(|k| grid.active[*k]).map(|k| h[k]).collect();
    if active.len() < 4 {
        return None;
    }
    let top = active.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let bottom = active.iter().cloned().fold(f64::INFINITY, f64::min);
    let depth = top - bottom;
    if !(depth > 0.0) {
        return None;
    }
    // width: the shorter run of columns that are near the deepest level
    let near = |k: usize| grid.active[k] && h[k] <= bottom + 0.1 * depth;
    let mut best = usize::MAX;
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            let k = j * grid.nx + i;
            if !near(k) { continue; }
            let mut wx = 1usize;
            let mut a = i;
            while a + 1 < grid.nx && near(j * grid.nx + a + 1) { wx += 1; a += 1; }
            let mut a = i;
            while a > 0 && near(j * grid.nx + a - 1) { wx += 1; a -= 1; }
            let mut wy = 1usize;
            let mut b = j;
            while b + 1 < grid.ny && near((b + 1) * grid.nx + i) { wy += 1; b += 1; }
            let mut b = j;
            while b > 0 && near((b - 1) * grid.nx + i) { wy += 1; b -= 1; }
            best = best.min(wx.min(wy));
        }
    }
    if best == usize::MAX {
        return None;
    }
    let width = best as f64 * grid.pitch;
    if !(width > 0.0) {
        return None;
    }
    let spec = atto_process::via::ViaFillSpec {
        diameter_nm: width,
        depth_nm: depth,
        mode: atto_process::via::FillMode::Conformal,
        sidewall_fraction: step_coverage.clamp(0.0, 1.0),
        rate_nm_per_min: 100.0,
        time_min: (metal_thickness_nm / 100.0).max(1e-9),
    };
    let fill = atto_process::via::fill(&spec).ok()?;
    Some((depth / width, width, fill))
}


/// Per-column film thickness from the 3-D transport solve, or `None` if it cannot run.
///
/// The scalar `step_coverage` states one number for the whole grid. This asks the solver what
/// each column actually receives, including shadowing and re-emission from nearby walls. The
/// result is normalised so an open surface keeps the requested thickness: the solver reports
/// flux relative to an open plane, which is exactly 1.0.
fn flux_scaled_thickness(grid: &ColumnGrid, base_nm: f64, sticking: f64, source: crate::params::DepositionSource) -> Option<Vec<f64>> {
    use atto_process::field3d::{self, Angular, Config, Source3d, VoxelGrid};
    if grid.nx == 0 || grid.ny == 0 || !(grid.pitch > 0.0) || !(base_nm > 0.0) {
        return None;
    }
    let h = grid.surface_heights();
    let hi = h.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let lo = h.iter().cloned().fold(f64::INFINITY, f64::min);
    let nz = ((((hi - lo).max(grid.pitch)) / grid.pitch).ceil() as usize + 4).clamp(4, 128);
    let heights: Vec<atto_units::Nm> = h.iter().map(|v| atto_units::Nm::new(*v)).collect();
    let vg = VoxelGrid::from_heights(grid.nx, grid.ny, nz, atto_units::Nm::new(grid.pitch), &heights).ok()?;
    let (angular, tilt_deg, spread_deg) = match source {
        crate::params::DepositionSource::Cosine => (Angular::Cosine, 0.0, 0.0),
        crate::params::DepositionSource::Beam { tilt_deg, spread_deg } => (Angular::Beam, tilt_deg, spread_deg),
    };
    let cfg = Config {
        rays: 20_000,
        seed: 0x5EED,
        sticking: sticking.clamp(0.0, 1.0),
        source: Source3d { angular, tilt_deg, azimuth_deg: 0.0, spread_deg },
    };
    let sol = field3d::solve(&vg, &cfg).ok()?;
    let mut per = vec![0.0f64; grid.cols.len()];
    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            // no flux means the solver did not cover this column; treating it as zero would be a
        // silent "nothing arrives here", so the whole solve is abandoned instead
        let Some(f) = sol.top_flux(ix, iy) else { return None; };
            per[k] = base_nm * f.clamp(0.0, 4.0);
        }
    }
    Some(per)
}

/// Polish the top film flat (P-01).
///
/// Uniform removal at the Preston rate, applied to every active column: the pad pressure is
/// taken as constant across the window. Pattern-dependent dishing and erosion need a line
/// scan, which `CmpModel::polish` provides and this stage does not yet build — so the reported
/// number is the blanket removal, not the topography.
fn run_cmp(w: &mut WaferState, win: &ProcessWindow, x: &crate::params::CmpParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    use atto_process::cmp::{CmpMaterial, CmpModel};
    let grid = ColumnGrid::sample(w, win)?;
    let model = CmpModel::bundled();
    let t = atto_units::Minutes::new(x.time_min * g(p, 0.1, 0.9));
    let cm = match x.material { Material::Metal => CmpMaterial::Metal, _ => CmpMaterial::Dielectric };
    let removal = model
        .removal(cm, atto_units::KPa::new(x.pressure_kpa), atto_units::MetresPerSecond::new(x.velocity_mps), t)
        .0;
    if !removal.is_finite() || removal < 0.0 {
        return Err(PipelineError::MissingInput("the polish removal is not a finite length".into()));
    }
    let mut polished = grid.clone();
    let mut touched = 0usize;
    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            if !grid.active[k] { continue; }
            let Some((z_top, _mat)) = geom::solid_top(&polished.cols[k]) else { continue };
            let z_new = (z_top - removal).max(polished.z_bottom);
            if z_new < z_top {
                geom::set_interval(&mut polished.cols[k], z_new, z_top, Material::Air);
                touched += 1;
            }
        }
    }
    polished.quantize();
    out.m("Removal", removal, "nm");
    out.m("Polished columns", touched as f64, "");
    out.m("Over-polish", x.over_polish * 100.0, "%");
    out.warn("uniform removal: dishing and erosion need a line scan, which this stage does not build");
    // No overlay to draw: the polish is uniform, so there is no pattern to show. The numbers
    // are in the report. `Overlay::Thickness` would need the per-column film thickness, which
    // this stage does not compute because it does not know which run is the film.
    let _ = &polished;
    Ok(Overlay::None)
}

fn run_deposition(w: &mut WaferState, win: &ProcessWindow, x: &crate::params::DepositionParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let grid = ColumnGrid::sample(w, win)?;
    let gr = g(p, 0.1, 0.8);
    // A transport solve per column when asked for; the scalar path otherwise. The scalar path
    // is what every existing recipe uses, so its numbers do not move.
    // A solve on a fixed geometry is right for one deposition and wrong for a growing film:
    // the surface that just grew changes what it can see. Splitting the deposition into
    // `flux_iterations` increments and re-solving between them is the fixed-point form of
    // that feedback. One iteration is the step-6 behaviour, kept as the default.
    let (grown, st) = if x.flux3d {
        let iters = x.flux_iterations.clamp(1, 64) as usize;
        let slice = x.thickness_nm * gr / iters as f64;
        let mut cur = grid.clone();
        let mut prev_mean = f64::NAN;
        let mut convergence = f64::NAN;
        let mut solved = true;
        let mut used = 0usize;
        for it in 0..iters {
            match flux_scaled_thickness(&cur, slice, x.sticking, x.source) {
                Some(per) => {
                    let mean = per.iter().sum::<f64>() / per.len().max(1) as f64;
                    if prev_mean.is_finite() && mean.abs() > 0.0 {
                        convergence = ((mean - prev_mean) / mean).abs();
                    }
                    prev_mean = mean;
                    let (next, _) = geom::grow_film_per_column(&cur, x.material, &per, x.step_coverage)?;
                    cur = next;
                    used = it + 1;
                    // stop when the film stops changing: a fixed count spends work that the
                    // tolerance has already declared unnecessary
                    if x.flux_tolerance > 0.0 && convergence.is_finite() && convergence < x.flux_tolerance {
                        break;
                    }
                }
                None => { solved = false; break; }
            }
        }
        if solved {
            out.m("Film model", 1.0, "");
            out.m("Flux iterations", used as f64, "");
            // A count equal to the maximum with a change still above the floor means the
            // answer is not converged. Saying so is the difference between a number and a
            // number you may trust.
            const CONVERGENCE_FLOOR: f64 = 0.02;
            let floor = x.flux_tolerance.max(CONVERGENCE_FLOOR);
            if used == iters && convergence.is_finite() && convergence > floor {
                out.warn(format!(
                    "the flux loop used all {} iterations and the mean flux was still changing by {:.3} per iteration; the film is not converged",
                    iters, convergence));
            }
            if convergence.is_finite() { out.m("Flux convergence", convergence, ""); }
            let before = grid.surface_heights();
            let after = cur.surface_heights();
            let total: Vec<f64> = (0..grid.cols.len()).map(|k| (after[k] - before[k]).max(0.0)).collect();
            let st = geom::film_stats_from_thickness(&grid, &total)?;
            (cur, st)
        } else {
            out.warn("the 3-D transport solve could not run on this grid; the scalar step coverage was used");
            geom::grow_film(&grid, x.material, x.thickness_nm * gr, x.step_coverage)?
        }
    } else {
        geom::grow_film(&grid, x.material, x.thickness_nm * gr, x.step_coverage)?
    };
    grown.write_back(&grid, w);
    out.m("Target thickness", x.thickness_nm, "nm");
    out.m("Realised mean thickness", st.mean_nm, "nm");
    out.m("Minimum thickness", st.min_nm, "nm");
    out.m("Maximum thickness", st.max_nm, "nm");
    out.m("Non-uniformity", st.nonuniformity_pct, "%");
    out.m("Octree cell (thickness quantum)", grid.pitch, "nm");
    if st.clipped_columns > 0 { out.warn(format!("{} columns hit the top of the high-resolution tile; film was clipped", st.clipped_columns)); }
    if p >= 1.0 {
        let tol = (0.1 * x.thickness_nm).max(grid.pitch * 0.5);
        if (st.mean_nm - x.thickness_nm).abs() > tol { out.warn(format!("realised thickness {:.1} nm differs from the {:.1} nm request by more than one octree cell", st.mean_nm, x.thickness_nm)); }
        if st.mean_nm <= 0.0 { return Err(PipelineError::MissingInput("no film was deposited: thickness is below half an octree cell".into())); }
    if x.material == Material::Metal {
        if let Some((ar, width, fill)) = via_fill_check(&grid, x.thickness_nm, x.step_coverage) {
            out.m("Opening aspect ratio", ar, "");
            out.m("Via fill fraction", fill.filled_fraction * 100.0, "%");
            out.m("Predicted void height", fill.void_height_nm, "nm");
            if fill.void && p >= 1.0 {
                out.warn(format!("a {:.0} nm wide, {:.0} nm deep opening would seal with a {:.0} nm void at this step coverage; a bottom-up fill or a wider via is needed", width, ar * width, fill.void_height_nm));
            }
        }
    }
        let r = deposition::run(&deposition::DepositionInputs { t_substrate: w.substrate_thickness_nm * 1e-9, t_film: st.mean_nm * 1e-9, bow_radius: x.bow_radius_m, ..Default::default() })
            .map_err(|e| PipelineError::Physics(e.to_string()))?;
        out.m("Film stress (Stoney)", r.stress_pa / 1e6, "MPa");
        if r.delamination_risk { out.warn("film stress exceeds the delamination threshold"); }
        w.meta.film_stress.push(FilmStress { material: x.material, thickness_nm: st.mean_nm, stress_mpa: r.stress_pa / 1e6 });
    }
    Ok(Overlay::Thickness { target_nm: x.thickness_nm, mean_nm: st.mean_nm, min_nm: st.min_nm, max_nm: st.max_nm })
}

fn run_litho(w: &mut WaferState, win: &ProcessWindow, x: &crate::params::LithoParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let grid = ColumnGrid::sample(w, win)?;
    if grid.cols.iter().zip(&grid.active).any(|(c, a)| *a && geom::solid_top(c).map(|s| s.1) == Some(Material::Photoresist)) {
        out.warn("wafer already carries photoresist; the new coat is stacked on top of it");
    }
    let (g1, g2, g3) = (g(p, 0.15, 0.35), g(p, 0.5, 0.65), g(p, 0.75, 0.95));
    let (coated, cs) = geom::coat_planarizing(&grid, x.resist_thickness_nm * g1)?;
    let r_nm = lithography::run_photo(&lithography::PhotoInputs { k1: x.k1, lambda: x.wavelength_nm, na: x.na }).map_err(|e| PipelineError::Invalid(e.to_string()))?;
    let sigma = r_nm / 2.355;
    let origin = [win.center[0] + x.align_offset_nm[0], win.center[1] + x.align_offset_nm[1]];
    // The mask and the window are both centred on the crop, but nothing enforces it: a mask
    // built at a different crop, or a window resized after the mask was made, leaves them
    // apart. Then the sampler sees no polygon and the field develops clear — which reads as
    // "dose too high" and is not.
    if let crate::geom::MaskPattern::Polygons { polys, .. } = &x.mask {
        let n = polys.len();
        out.m("Mask polygons", n as f64, "");
        if n == 0 {
            out.warn("the mask has no polygon at all: every column sees clear field");
        } else {
            let mut lo = [f64::INFINITY, f64::INFINITY];
            let mut hi = [f64::NEG_INFINITY, f64::NEG_INFINITY];
            for p in polys {
                for v in p {
                    lo[0] = lo[0].min(v[0]); lo[1] = lo[1].min(v[1]);
                    hi[0] = hi[0].max(v[0]); hi[1] = hi[1].max(v[1]);
                }
            }
            let (mx0, mx1) = (lo[0] + origin[0], hi[0] + origin[0]);
            let (my0, my1) = (lo[1] + origin[1], hi[1] + origin[1]);
            let (wx0, wx1) = (win.center[0] - win.size[0] * 0.5, win.center[0] + win.size[0] * 0.5);
            let (wy0, wy1) = (win.center[1] - win.size[1] * 0.5, win.center[1] + win.size[1] * 0.5);
            if !(mx1 > wx0 && mx0 < wx1 && my1 > wy0 && my0 < wy1) {
                out.warn(format!(
                    "the mask spans x {:.0}..{:.0} nm but the window is x {:.0}..{:.0} nm: they do not overlap, so every column sees clear field. Rebuild the mask for this crop.",
                    mx0, mx1, wx0, wx1));
            }
        }
    }
    let img = geom::aerial_image(&grid, &x.mask, origin, sigma)?;
    // A positive resist develops at a rate that depends on dose, so the depth varies
    // across a feature edge and the profile gets a slope. A negative tone keeps the
    // threshold path: the model above is for positive resist only.
    let model = atto_process::resist::DevelopModel::default();
    let dill = atto_process::dill::DillModel::default();
    // bake the latent image: the compound diffuses, so a step becomes a ramp
    let sigma_bake = atto_process::dill::bake_sigma_from_temperature(atto_units::Kelvin::new(x.peb_c + 273.15), atto_units::Minutes::new(1.0)).0;
    let baked = atto_process::dill::bake_blur(&img, grid.nx, grid.ny, sigma_bake, grid.pitch);
    out.m("Post-exposure bake diffusion length", sigma_bake, "nm");
    out.m("Develop time", x.develop_time_min * g3.max(0.0), "min");
    // The developer works through the film as designed. The grid stores that film rounded to
    // whole cells (a 150 nm coat on a 64 nm grid is thicker), so the depth reached is turned into
    // a fraction of the film and applied to the stored thickness; otherwise a coarse cell leaves a
    // sliver at the bottom of every cleared column and nothing ever opens.
    let stored_nm = cs.mean_nm.max(x.resist_thickness_nm);
    let scale = stored_nm / x.resist_thickness_nm.max(1e-9);
    let depths: Vec<f64> = if x.positive_tone {
        geom::develop_depth_map_dill(&baked, x.dose_ratio, &model, &dill, x.develop_time_min * g3.max(0.0), x.resist_thickness_nm, 16)
            .into_iter().map(|d| d * scale).collect()
    } else {
        Vec::new()
    };
    let removed: Vec<bool> = if x.positive_tone {
        depths.iter().map(|d| *d >= stored_nm - 1e-6).collect()
    } else {
        geom::develop_map(&img, x.dose_ratio, false)
    };
    let (dev, ds) = if x.positive_tone {
        geom::develop_depths(&coated, &depths)?
    } else {
        geom::develop(&coated, &removed, g3)?
    };
    dev.write_back(&grid, w);
    out.m("Rayleigh resolution", r_nm, "nm");
    out.m("Optical blur σ", sigma, "nm");
    out.m("Mean resist thickness", cs.mean_nm, "nm");
    // P-04: a planarity error is a focus error, and a focus error prints. The depth of focus
    // here is the Rayleigh form with k2 = 0.5; beyond it the pattern is not expected to hold.
    let focus_error = atto_process::coupling::defocus_from_erosion(atto_units::Nm::new(x.planarity_error_nm));
    if x.planarity_error_nm.abs() > 0.0 {
        out.m("Planarity focus error", focus_error.0, "nm");
        let dof = 0.5 * x.wavelength_nm / (x.na * x.na).max(1e-9);
        if focus_error.0.abs() > dof && p >= 1.0 {
            out.warn(format!(
                "a {:.0} nm planarity error is {:.0} nm of defocus, beyond the {:.0} nm depth of focus: this field will not print to CD",
                x.planarity_error_nm, focus_error.0.abs(), dof));
        }
    }
    out.m("Opened columns", ds.opened_columns as f64, "");
    out.m("Retained columns", ds.retained_columns as f64, "");
    out.m("Partially developed columns", ds.partial_columns as f64, "");
    out.m("Mask alignment offset", (x.align_offset_nm[0].powi(2) + x.align_offset_nm[1].powi(2)).sqrt(), "nm");
    if p >= 1.0 {
        let cd = measure::line_widths(&dev, Material::Photoresist);
        if !cd.is_empty() { out.m("Resist CD (mean line width)", cd.iter().sum::<f64>() / cd.len() as f64, "nm"); }
        if ds.opened_columns == 0 { out.warn("development opened no resist — nothing downstream can be etched or implanted; raise the dose, check the mask polarity, or pick a layer with openings"); }
        if ds.retained_columns == 0 { out.warn("development removed all resist: dose is too high or the mask is fully clear"); }
    }
    Ok(if p >= 0.45 {
        Overlay::Exposure(ExposureOverlay { x0: grid.x0, y0: grid.y0, pitch: grid.pitch, nx: grid.nx, ny: grid.ny,
            intensity: img.iter().map(|v| v * g2 as f32).collect(), removed: if g3 > 0.0 { removed } else { vec![] } })
    } else { Overlay::None })
}


/// Per-column rate multiplier from the etch model: loading plus aspect ratio.
///
/// `loading_factor` accounts for the whole open area consuming reactant; `arde_multiplier`
/// accounts for how hard it is for fresh reactant to reach the bottom of a narrow opening.
/// Both are properties of the plasma, so they come from `atto_process` rather than being
/// invented here.
fn arde_scale(grid: &ColumnGrid, spec: &EtchSpec) -> Vec<f64> {
    let model = atto_process::etch::EtchModel::bundled_oxide();
    let n = grid.cols.len();
    let mut scale = vec![1.0f64; n];
    let open = grid.cols.iter().zip(&grid.active).filter(|(c, a)| **a && geom::solid_top(c).map(|(_, m)| m != Material::Photoresist && spec.rates.rate(m) > 0.0).unwrap_or(false)).count();
    let load = model.loading_factor(open as f64 / n.max(1) as f64);
    for j in 0..grid.ny {
        for i in 0..grid.nx {
            let k = j * grid.nx + i;
            if !grid.active[k] {
                continue;
            }
            // width of the opening this column sits in: the shorter run through it
            let mut wx = 1usize;
            let mut a = i;
            while a + 1 < grid.nx && grid.active[j * grid.nx + a + 1] { wx += 1; a += 1; }
            let mut a = i;
            while a > 0 && grid.active[j * grid.nx + a - 1] { wx += 1; a -= 1; }
            let mut wy = 1usize;
            let mut b = j;
            while b + 1 < grid.ny && grid.active[(b + 1) * grid.nx + i] { wy += 1; b += 1; }
            let mut b = j;
            while b > 0 && grid.active[(b - 1) * grid.nx + i] { wy += 1; b -= 1; }
            let width = (wx.min(wy) as f64) * grid.pitch;
            if !(width > 0.0) {
                continue;
            }
            // how deep this column would reach at the blanket rate
            let top_rate = geom::solid_top(&grid.cols[k]).map(|(_, m)| spec.rates.rate(m)).unwrap_or(0.0);
            let depth = (spec.time_min * top_rate).max(0.0);
            scale[k] = model.arde_multiplier(depth / width) * load;
        }
    }
    scale
}


#[cfg(test)]
mod arde_tests {
    //! FR-3.1 (depth = time x rate) and FR-3.2 (aspect ratio and loading), PRD defect D-2.
    //!
    //! Two things the first version of these tests got wrong, kept here so they are not repeated.
    //! Loading is one factor for the whole grid, so comparing a narrow slot in one grid with a wide
    //! slot in another compares two different loadings and says nothing about aspect ratio: both
    //! slots must sit in the same grid. And ARDE only acts above the critical aspect ratio (8 in
    //! the bundled pack), so a 100 nm etch into a 20 nm slot (aspect ratio 5) is not deep enough to
    //! show it: the rate must make the slot truly deep.
    use super::*;

    /// 64 x 64 columns of 500 nm silicon, 10 nm cells, with two slots running the full height of y,
    /// the narrow one starting at x = 4 cells and the wide one at x = 24 cells.
    fn grid_with_two_slots(narrow: usize, wide: usize) -> ColumnGrid {
        let (nx, ny) = (64usize, 64usize);
        let mut cols = Vec::with_capacity(nx * ny);
        for _ in 0..nx * ny {
            cols.push(vec![crate::wafer::Run { z0: -500.0, z1: 0.0, material: Material::Silicon }]);
        }
        let mut active = vec![false; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                if (i >= 4 && i < 4 + narrow) || (i >= 24 && i < 24 + wide) {
                    active[j * nx + i] = true;
                }
            }
        }
        ColumnGrid { x0: 0.0, y0: 0.0, pitch: 10.0, nx, ny, cols, active, z_bottom: -500.0, z_top: 0.0 }
    }

    /// One minute at `rate` nm/min in every material, perfectly vertical.
    fn spec(rate: f64) -> EtchSpec {
        EtchSpec {
            rates: crate::geom::EtchRates { si: rate, sio2: rate, sic: rate, gan: rate, resist: 0.0, metal: 0.0 },
            anisotropy: 1.0,
            time_min: 1.0,
        }
    }

    /// How far the surface of column `i`, in the middle row, dropped.
    fn removed_at(g: &ColumnGrid, spec: &EtchSpec, scale: Option<&[f64]>, i: usize) -> f64 {
        let (out, _) = geom::etch_scaled(g, spec, 1.0, scale).expect("etch runs");
        let k = (g.ny / 2) * g.nx + i;
        let before = geom::solid_top(&g.cols[k]).map(|(z, _)| z).unwrap_or(0.0);
        let after = geom::solid_top(&out.cols[k]).map(|(z, _)| z).unwrap_or(0.0);
        before - after
    }

    #[test]
    fn fr_3_1_depth_equals_time_times_rate_within_one_cell() {
        let g = grid_with_two_slots(2, 20);
        let removed = removed_at(&g, &spec(100.0), None, 5);
        assert!((removed - 100.0).abs() <= g.pitch, "1 min at 100 nm/min removes 100 nm within one 10 nm cell: {removed}");
    }

    #[test]
    fn fr_3_2_a_narrow_deep_slot_etches_measurably_slower_than_a_wide_one_in_the_same_wafer() {
        // 400 nm/min for one minute: a 400 nm blanket depth, so aspect ratio 20 in the 20 nm slot
        // (far past the critical 8) and 2 in the 200 nm slot (well under it)
        let g = grid_with_two_slots(2, 20);
        let s = spec(400.0);
        let scale = arde_scale(&g, &s);
        let (narrow, wide) = (removed_at(&g, &s, Some(&scale), 4), removed_at(&g, &s, Some(&scale), 34));
        assert!(wide > 0.0 && narrow > 0.0, "both slots etch: narrow {narrow}, wide {wide}");
        assert!(narrow < 0.6 * wide, "a narrow deep slot etches measurably slower: narrow {narrow} nm, wide {wide} nm");
        // and by the amount the model says: loading is shared, ARDE retains 40% at twice the critical ratio
        let load = atto_process::etch::EtchModel::bundled_oxide().loading_factor(22.0 / 64.0);
        assert!((wide - 400.0 * load).abs() <= g.pitch, "wide slot: {wide} nm, expected {}", 400.0 * load);
        assert!((narrow - 400.0 * 0.4 * load).abs() <= g.pitch, "narrow slot: {narrow} nm, expected {}", 400.0 * 0.4 * load);
    }

    #[test]
    fn a_slot_below_the_critical_aspect_ratio_is_not_penalised_for_being_narrow() {
        // 100 nm into a 20 nm slot is aspect ratio 5, under the critical 8: only loading acts, and
        // loading is the same for every column, so the narrow and wide slots etch alike
        let g = grid_with_two_slots(2, 20);
        let scale = arde_scale(&g, &spec(100.0));
        let (n, w) = (scale[(g.ny / 2) * g.nx + 4], scale[(g.ny / 2) * g.nx + 34]);
        assert!((n - w).abs() < 1e-12, "no ARDE below the critical aspect ratio: narrow {n}, wide {w}");
    }

    #[test]
    fn a_shallow_blanket_etch_is_scaled_by_loading_alone() {
        // every column open: aspect ratio ~0, so ARDE contributes nothing and the whole factor is loading
        let (nx, ny) = (64usize, 64usize);
        let g = ColumnGrid {
            x0: 0.0, y0: 0.0, pitch: 10.0, nx, ny,
            cols: (0..nx * ny).map(|_| vec![crate::wafer::Run { z0: -500.0, z1: 0.0, material: Material::Silicon }]).collect(),
            active: vec![true; nx * ny],
            z_bottom: -500.0, z_top: 0.0,
        };
        let scale = arde_scale(&g, &spec(100.0));
        let load = atto_process::etch::EtchModel::bundled_oxide().loading_factor(1.0);
        assert!((scale[(ny / 2) * nx + nx / 2] - load).abs() < 1e-12, "{} vs loading {load}", scale[(ny / 2) * nx + nx / 2]);
    }
}

fn run_etch(w: &mut WaferState, win: &ProcessWindow, x: &crate::params::EtchParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let grid = ColumnGrid::sample(w, win)?;
    let eligible = grid.cols.iter().zip(&grid.active).any(|(c, a)| *a && geom::solid_top(c).map(|s| x.rates.rate(s.1) > 0.0).unwrap_or(false));
    if !eligible { return Err(PipelineError::MissingInput("nothing in the window is etchable in this chemistry".into())); }
    // audit F-02: cross-check the tabulated rate against the model rather
    // than against itself. A disagreement is a warning, not a silent pass.
    if let Some((v_model, _)) = atto_bridge::etch_rates("oxide", atto_units::Watts::new(x.rf_power_w), atto_units::MTorr::new(x.pressure_mtorr)) {
        let v_table = x.rates.rate(Material::SiliconDioxide);
        if v_table > 0.0 && (v_model / v_table > 2.0 || v_table / v_model > 2.0) {
            out.warn("the tabulated etch rate differs from the model by more than a factor of two");
        }
    }
    let spec = EtchSpec { rates: x.rates, anisotropy: x.anisotropy, time_min: x.time_min };
    // Aspect-ratio-dependent etching: the flux reaching the bottom of a narrow, deep
    // opening is reduced, so it etches slower than a wide shallow one. Without this a via
    // bores straight down at the blanket rate, which no real plasma does.
    let scale = arde_scale(&grid, &spec);
    let (etched, st) = geom::etch_scaled(&grid, &spec, g(p, 0.25, 0.8), Some(&scale))?;
    let mut cur = etched;
    if x.strip_after { cur = geom::strip_resist(&cur, 1.0, g(p, 0.85, 1.0))?.0; }
    cur.write_back(&grid, w);
    out.m("Etch time", x.time_min, "min");
    out.m("Maximum etch depth", st.max_depth_nm, "nm");
    out.m("Mean depth in openings", st.mean_depth_open_nm, "nm");
    out.m("Mask consumed (mean)", st.mask_consumed_mean_nm, "nm");
    out.m("Mask fully consumed", st.mask_failed_columns as f64, "columns");
    out.m("Lateral undercut / recession", st.max_lateral_nm, "nm");
    if st.mask_consumed_mean_nm > 0.0 { out.m("Selectivity (measured)", st.mean_depth_open_nm / st.mask_consumed_mean_nm, ":1"); }
    if p >= 1.0 {
        // Analytic cross-check with the closed-form etch model.
        let mut tops: std::collections::HashMap<Material, usize> = Default::default();
        for (c, a) in grid.cols.iter().zip(&grid.active) { if *a { if let Some((_, m)) = geom::solid_top(c) { if m != Material::Photoresist { *tops.entry(m).or_default() += 1; } } } }
        if let Some((m, _)) = tops.into_iter().max_by_key(|(_, n)| *n) {
            let rv = x.rates.rate(m);
            if rv > 0.0 {
                if let Ok(a) = etching::run(&etching::EtchingInputs { r_lateral: (1.0 - x.anisotropy) * rv, r_vertical: rv, er_target: rv, er_mask: x.rates.resist.max(1e-6), target_depth_nm: rv * x.time_min, mask_thickness_nm: 150.0, ..Default::default() }) {
                    out.m("Selectivity (analytic model)", a.selectivity, ":1");
                    out.m("Anisotropy (analytic model)", a.anisotropy, "");
                }
            }
        }
        if st.max_depth_nm <= 0.0 { out.warn("nothing was removed: the etch time is too short for the exposed material"); }
        if st.mask_failed_columns > 0 { out.warn(format!("mask eroded through in {} columns; pattern transfer is compromised", st.mask_failed_columns)); }
    }
    Ok(Overlay::None)
}

// --- doping ---------------------------------------------------------------

struct ProfileModel { rp: f64, sigma0: f64, sigma_eff: f64, peak: f64, xj: f64, diff_cm2s: f64 }

fn interp_range(table: &[(f64, f64, f64)], e: f64) -> (f64, f64) {
    let (le, n) = (e.ln(), table.len());
    let (i, j) = if e <= table[0].0 { (0, 1) } else if e >= table[n - 1].0 { (n - 2, n - 1) } else { let k = table.iter().position(|t| t.0 >= e).unwrap(); (k - 1, k) };
    let f = (le - table[i].0.ln()) / (table[j].0.ln() - table[i].0.ln());
    let lerp = |a: f64, b: f64| (a.ln() + f * (b.ln() - a.ln())).exp();
    (lerp(table[i].1, table[j].1), lerp(table[i].2, table[j].2))
}

/// Built-in projected range (nm) for a species at an energy (log-log interpolation of the LSS-style table).
pub fn implant_range_nm(species: crate::params::Species, energy_kev: f64) -> f64 { interp_range(species.range_table(), energy_kev).0 }

fn profile_model(x: &DopingParams, g1: f64, g2: f64) -> ProfileModel {
    let (rp, sigma0) = match x.rp_override_nm { Some(r) => (r, 0.35 * r), None => interp_range(x.species.range_table(), x.energy_kev) };
    let (d0, ea) = x.species.diffusion();
    let t_k = x.anneal_temp_c + 273.15;
    let diff = d0 * (-ea / (8.617_333e-5 * t_k)).exp();
    let sigma_eff = (sigma0 * sigma0 + 2.0 * diff * x.anneal_min * 60.0 * g2 * 1.0e14).sqrt();
    let peak = x.dose_cm2 * g1 / ((2.0 * std::f64::consts::PI).sqrt() * sigma_eff * 1e-7);
    let xj = if peak > x.background_cm3 { rp + sigma_eff * (2.0 * (peak / x.background_cm3).ln()).sqrt() } else { 0.0 };
    ProfileModel { rp, sigma0, sigma_eff, peak, xj, diff_cm2s: diff }
}

fn conc(m: &ProfileModel, d_nm: f64) -> f64 { m.peak * (-(d_nm - m.rp).powi(2) / (2.0 * m.sigma_eff * m.sigma_eff)).exp() }

/// Concentration → dopant bucket: ≥1e17 cm⁻³ is bucket 1, one bucket per decade, 6 max.
pub fn dopant_bucket(n: f64) -> u8 { if n < 1e17 { 0 } else { ((n.log10().floor() - 16.0).clamp(1.0, 6.0)) as u8 } }

fn stopping(m: Material) -> f64 {
    // No wildcard: a new material must be given a stopping power rather than
    // silently becoming free passage (audit F-11).
    match m {
        Material::SiliconDioxide | Material::SiliconCarbide => 1.0,
        Material::Photoresist => 0.9,
        Material::Metal => 2.5,
        Material::GalliumNitride => 1.4,
        Material::Air | Material::Silicon | Material::Dopant(_) => 0.0,
    }
}

fn run_doping(w: &mut WaferState, win: &ProcessWindow, x: &DopingParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let grid = ColumnGrid::sample(w, win)?;
    let (g1, g2) = (g(p, 0.1, 0.5), g(p, 0.55, 0.9));
    let pm = profile_model(x, g1.max(1e-12), g2);
    let mut newg = grid.clone();
    let mut implanted = 0usize;
    let mut deepest = 0.0f64;
    let mut truncated = false;
    if g1 > 0.0 {
        for k in 0..grid.cols.len() {
            if !grid.active[k] { continue; }
            let col = &grid.cols[k];
            let (mut cover, mut found) = (0.0, None);
            for r in col.iter().rev() {
                match r.material {
                    Material::Air => continue,
                    Material::Silicon | Material::Dopant(_) => { found = Some(*r); break; }
                    m => cover += (r.z1 - r.z0) * stopping(m),
                }
            }
            let Some(si) = found else { continue };
            let cells = ((si.z1 - si.z0) / grid.pitch) as usize;
            // the implant profile comes from the crates: channelling, straggle and tilt
            let depths: Vec<f64> = (0..cells).map(|c| (c as f64 + 0.5) * grid.pitch + cover).collect();
            let concs = atto_bridge::implant_column(
                &depths,
                x.species.label(),
                x.energy_kev,
                x.dose_cm2,
                x.tilt_deg,
            )
            .map_err(|e| PipelineError::MissingInput(format!("{e}")))?;
            let mut hit = false;
            for c in 0..cells {
                let d = (c as f64 + 0.5) * grid.pitch + cover;
                let b = dopant_bucket(concs[c]);
                if b == 0 { if d > pm.rp { break; } continue; }
                let (zhi, zlo) = (si.z1 - c as f64 * grid.pitch, si.z1 - (c + 1) as f64 * grid.pitch);
                let keep = match geom::material_at_z(col, 0.5 * (zhi + zlo)) { Material::Dopant(e) => e >= b, Material::Silicon => false, _ => true };
                if keep { continue; }
                geom::set_interval(&mut newg.cols[k], zlo, zhi, Material::Dopant(b));
                hit = true; deepest = deepest.max(si.z1 - zlo);
                if c + 1 == cells { truncated = true; }
            }
            if hit { implanted += 1; }
        }
    }
    newg.write_back(&grid, w);
    let full = profile_model(x, 1.0, 1.0);
    out.m("Projected range Rp", full.rp, "nm");
    out.m("Straggle ΔRp (as implanted)", full.sigma0, "nm");
    out.m("Straggle after anneal", full.sigma_eff, "nm");
    out.m("Peak concentration", full.peak, "cm⁻³");
    out.m("Junction depth (analytic)", full.xj, "nm");
    out.m("Diffusivity at anneal T", full.diff_cm2s, "cm²/s");
    out.m("Implanted columns", implanted as f64, "");
    if p >= 1.0 {
        if implanted == 0 { return Err(PipelineError::MissingInput("no silicon is exposed to the implant: the oxide / resist stack stops every ion. If the lithography stage opened no resist, fix that stage first — this one cannot recover from it".into())); }
        if truncated || full.xj > (w.surface_z_nm - w.bounds.min.2) { out.warn("profile reaches the bottom of the high-resolution tile and was truncated"); }
        out.m("Deepest dopant cell (octree)", deepest, "nm");
        w.meta.dopant_profiles.push(DopantProfile { species: x.species.label().into(), dose_cm2: x.dose_cm2, energy_kev: x.energy_kev, rp_nm: full.rp, sigma_nm: full.sigma_eff,
            anneal_temp_c: x.anneal_temp_c, anneal_min: x.anneal_min, peak_cm3: full.peak, junction_depth_nm: full.xj });
        if full.peak < x.background_cm3 { out.warn("peak dopant concentration is below the background doping: no junction forms"); }
    }
    Ok(Overlay::None)
}

fn run_clean(w: &mut WaferState, win: &ProcessWindow, x: &crate::params::CleanParams, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let grid = ColumnGrid::sample(w, win)?;
    let before_v = geom::material_volume_nm3(&grid, |m| m == Material::Photoresist);
    let (stripped, st) = geom::strip_resist(&grid, x.strip_efficiency, g(p, 0.0, 0.45))?;
    stripped.write_back(&grid, w);
    let residue = if before_v > 0.0 { st.resist_after_nm3 / st.resist_before_nm3.max(1e-30) } else { 0.0 };
    let particles = w.meta.surface.particle_density_cm2;
    let particles_after = if x.chemical_clean { particles * (1.0 - x.particle_removal_eff) } else { particles };
    out.m("Resist before strip", before_v / (grid.active_count().max(1) as f64 * grid.cell_area()), "nm (mean)");
    out.m("Resist residue (fraction left)", geom::material_volume_nm3(&stripped, |m| m == Material::Photoresist) / before_v.max(1e-30) * if before_v > 0.0 { 1.0 } else { 0.0 }, "");
    out.m("Particle density after", particles_after, "cm⁻²");
    if before_v == 0.0 { out.warn("no photoresist in the window; nothing to strip"); }
    if p >= 1.0 {
        w.meta.surface.residue_fraction = if before_v > 0.0 { geom::material_volume_nm3(&stripped, |m| m == Material::Photoresist) / before_v } else { 0.0 };
        w.meta.surface.particle_density_cm2 = particles_after;
        if residue > 0.01 { out.warn(format!("{:.1}% of the resist survived the strip", residue * 100.0)); }
    }
    Ok(Overlay::Surface { clean_fraction: g(p, 0.45, 0.9), particles_cm2: particles * (1.0 - (particles - particles_after) / particles.max(1e-30) * g(p, 0.45, 0.9)) })
}

fn run_metrology(w: &mut WaferState, win: &ProcessWindow, x: &crate::params::MetrologyParams, seed: u64, p: f64, out: &mut StageOutcome) -> Result<Overlay, PipelineError> {
    let report = measure::measure(w, win, x, seed)?;
    for m in &report.measurements {
        if let Some(v) = m.value { out.metrics.push(Metric { name: m.name.clone(), value: v, unit: if m.unit == "MPa" { "MPa" } else { "nm" } }); }
        match m.status {
            Status::Fail => out.warn(format!("{} = {:.1} is outside {:.1} ± {:.1}", m.name, m.value.unwrap_or(f64::NAN), m.nominal, m.tolerance)),
            Status::NotMeasured => out.warn(format!("{} could not be measured on this wafer state", m.name)),
            Status::Pass => {}
        }
    }
    if let Some(o) = &report.optical { out.m("Ellipsometric Ψ", o.psi_deg, "°"); out.m("Ellipsometric Δ", o.delta_deg, "°"); }
    let n = report.wafer_map.n;
    let (lo, hi) = win.bounds_xy();
    let shown = (g(p, 0.0, 0.8) * (n * n) as f64).ceil() as usize;
    let sites = (0..n * n).take(shown).map(|i| {
        let (si, sj) = (i % n, i / n);
        Site { x: lo[0] + (si as f64 + 0.5) / n as f64 * (hi[0] - lo[0]), y: lo[1] + (sj as f64 + 0.5) / n as f64 * (hi[1] - lo[1]), value: report.wafer_map.values[i], out_of_spec: report.wafer_map.out_of_spec(i) }
    }).collect();
    out.report = Some(report);
    Ok(Overlay::Sites(sites))
}

/// Convenience: evaluate a list of stages in order at full progress,
/// threading the wafer through. Used by tests and headless tooling.
pub fn run_chain(initial: &WaferState, window: &ProcessWindow, stages: &[StageParams], seed: u64) -> Result<Vec<StageEval>, (usize, PipelineError)> {
    let mut cur = initial.clone();
    let mut out = Vec::new();
    for (i, s) in stages.iter().enumerate() {
        let ev = evaluate_stage(&cur, window, s, seed, 1.0).map_err(|e| (i, e))?;
        cur = ev.wafer.clone();
        out.push(ev);
    }
    Ok(out)
}

#[allow(dead_code)]
fn _keep(_: Run) {}

#[cfg(test)]
mod tests {
    use super::*;

    pub fn tile() -> WaferState { WaferState::new_wafer_tile(300.0, 500.0, 16.0, 7, (0.0, 0.0)).unwrap() }
    fn win() -> ProcessWindow { ProcessWindow::default() }
    fn defaults() -> Vec<StageParams> { Stage::DEFAULT_RECIPE.iter().map(|s| StageParams::default_for(*s)).collect() }

    #[test]
    fn default_recipe_runs_end_to_end_and_stays_sparse() {
        let chain = run_chain(&tile(), &win(), &defaults(), 1).map_err(|(i, e)| format!("stage {i}: {e}")).unwrap();
        assert_eq!(chain.len(), 7);
        for ev in &chain { assert!(!ev.outcome.has_errors(), "{:?}", ev.outcome.issues); }
        let last = &chain[6];
        eprintln!("final nodes: {}", last.wafer.node_count());
        assert!(last.wafer.node_count() < 120_000, "SVO stays sparse: {} nodes", last.wafer.node_count());
        let rep = last.outcome.report.as_ref().unwrap();
        for m in &rep.measurements { eprintln!("{:40} {:?} nominal {} ± {} -> {:?}", m.name, m.value, m.nominal, m.tolerance, m.status); }
        assert_eq!(rep.wafer_fingerprint, chain[5].wafer.fingerprint(), "PROC-04: metrology references the final committed wafer");
        // The wafer map judges the oxide where it exists against the 100 ± 10 nm target: a correct film must not read as out of spec.
        let wm = &rep.wafer_map;
        assert!(wm.label.contains("mean where present"), "{}", wm.label);
        let (have, bad) = (wm.values.iter().flatten().count(), (0..wm.values.len()).filter(|i| wm.out_of_spec(*i)).count());
        assert!(have > 0 && bad == 0, "{bad} of {have} measured sites flagged out of spec: {:?}", wm.values);
    }

    #[test]
    fn a_recipe_saved_before_develop_time_existed_loads_with_the_old_one_minute_behaviour() {
        let mut v = serde_json::to_value(StageParams::default_for(Stage::Lithography)).unwrap();
        let obj = v.get_mut("Lithography").and_then(|x| x.as_object_mut()).expect("the enum is externally tagged");
        assert!(obj.remove("develop_time_min").is_some(), "the field is serialised today");
        let back: StageParams = serde_json::from_value(v).unwrap();
        match back {
            StageParams::Lithography(l) => assert_eq!(l.develop_time_min, 1.0),
            other => panic!("expected lithography, got {other:?}"),
        }
    }

    #[test]
    fn a_develop_time_outside_its_range_is_refused() {
        for bad in [0.0, -1.0, 61.0, f64::NAN, f64::INFINITY] {
            let mut p = StageParams::default_for(Stage::Lithography);
            if let StageParams::Lithography(l) = &mut p { l.develop_time_min = bad; }
            assert!(p.validate().iter().any(|i| i.severity == Severity::Error), "develop time {bad} must be refused");
        }
        assert!(StageParams::default_for(Stage::Lithography).validate().iter().all(|i| i.severity != Severity::Error));
    }

    fn litho_with_develop_time(t: f64) -> StageParams {
        let mut p = StageParams::default_for(Stage::Lithography);
        if let StageParams::Lithography(l) = &mut p { l.develop_time_min = t; }
        p
    }

    fn wafer_ready_for_litho() -> WaferState {
        let mut cur = tile();
        for s in defaults().iter().take(2) { cur = evaluate_stage(&cur, &win(), s, 3, 1.0).unwrap().wafer; }
        cur
    }

    #[test]
    fn d_5_the_recipe_sets_the_develop_time_and_the_animation_only_scales_it() {
        let cur = wafer_ready_for_litho();
        let short = evaluate_stage(&cur, &win(), &litho_with_develop_time(0.4), 3, 1.0).unwrap();
        let long = evaluate_stage(&cur, &win(), &litho_with_develop_time(2.5), 3, 1.0).unwrap();
        assert_eq!(short.outcome.metric("Develop time"), Some(0.4));
        assert_eq!(long.outcome.metric("Develop time"), Some(2.5));
        let partway = evaluate_stage(&cur, &win(), &litho_with_develop_time(2.0), 3, 0.85).unwrap();
        let t = partway.outcome.metric("Develop time").unwrap();
        assert!(t > 0.0 && t < 2.0, "partway through the animation only part of the time has elapsed: {t}");
    }

    #[test]
    fn d_5_a_much_longer_develop_changes_the_developed_geometry() {
        // if this fails the develop rate is so high that even 0.01 min clears the film: raise the short time
        let cur = wafer_ready_for_litho();
        let barely = evaluate_stage(&cur, &win(), &litho_with_develop_time(0.01), 3, 1.0).unwrap();
        let fully = evaluate_stage(&cur, &win(), &litho_with_develop_time(5.0), 3, 1.0).unwrap();
        assert_ne!(barely.wafer.fingerprint(), fully.wafer.fingerprint(), "develop time must change what the developer removed");
    }

    #[test]
    fn last_animation_frame_equals_the_committed_geometry() {
        let mut cur = tile();
        for s in defaults().iter() {
            let full = evaluate_stage(&cur, &win(), s, 3, 1.0).unwrap();
            let frame = evaluate_stage(&cur, &win(), s, 3, frame_progress(FRAMES_PER_STAGE)).unwrap();
            assert_eq!(full.wafer.fingerprint(), frame.wafer.fingerprint(), "{:?}", s.stage());
            cur = full.wafer;
        }
    }

    #[test]
    fn frames_are_deterministic_and_deposition_never_grows_before_the_growth_phase() {
        let s = StageParams::default_for(Stage::Deposition);
        let a = evaluate_stage(&tile(), &win(), &s, 9, 0.5).unwrap();
        let b = evaluate_stage(&tile(), &win(), &s, 9, 0.5).unwrap();
        assert_eq!(a.wafer.fingerprint(), b.wafer.fingerprint());
        let early = evaluate_stage(&tile(), &win(), &s, 9, 0.05).unwrap();
        assert_eq!(early.wafer.fingerprint(), tile().fingerprint(), "recipe-entry phase must not change geometry");
        assert!(!early.outcome.geometry_changed);
    }

    #[test]
    fn prep_and_metrology_do_not_change_geometry() {
        let w = tile();
        for st in [Stage::SubstratePrep, Stage::Metrology] {
            let ev = evaluate_stage(&w, &win(), &StageParams::default_for(st), 1, 1.0).unwrap();
            assert!(!ev.outcome.geometry_changed || st == Stage::SubstratePrep, "{st:?}");
            assert_eq!((ev.wafer.bounds.min, ev.wafer.bounds.max), (w.bounds.min, w.bounds.max));
        }
        let m = evaluate_stage(&w, &win(), &StageParams::default_for(Stage::Metrology), 1, 1.0).unwrap();
        assert_eq!(m.wafer.fingerprint(), w.fingerprint());
    }

    #[test]
    fn proc01_reported_thickness_matches_the_committed_geometry() {
        let ev = evaluate_stage(&tile(), &win(), &StageParams::default_for(Stage::Deposition), 1, 1.0).unwrap();
        let reported = ev.outcome.metric("Realised mean thickness").unwrap();
        let rep = measure::measure(&ev.wafer, &win(), &crate::params::MetrologyParams { targets: vec![crate::params::Target { quantity: crate::params::Quantity::FilmThickness(Material::SiliconDioxide), nominal: 100.0, tolerance: 10.0 }], sites: 3, noise_nm: 0.0, wavelength_nm: 633.0, incidence_deg: 70.0 }, 0).unwrap();
        let measured = rep.measurements[0].value.unwrap();
        assert!((reported - measured).abs() < 1e-6, "reported {reported} vs measured {measured}");
        assert!((reported - 100.0).abs() <= 8.0);
    }

    #[test]
    fn proc02_mask_openings_create_the_resist_pattern() {
        let s = run_chain(&tile(), &win(), &defaults()[..3], 1).unwrap();
        let w = &s[2].wafer;
        let g = ColumnGrid::sample(w, &win()).unwrap();
        let row = g.ny / 2;
        let line = g.idx(g.nx / 2, row);
        let space = g.idx(g.nx / 2 + (160.0 / g.pitch) as usize, row);
        assert_eq!(geom::solid_top(&g.cols[line]).unwrap().1, Material::Photoresist, "resist under the opaque mask line");
        assert_eq!(geom::solid_top(&g.cols[space]).unwrap().1, Material::SiliconDioxide, "oxide exposed in the clear space");
        let cd = s[2].outcome.metric("Resist CD (mean line width)").unwrap();
        assert!((cd - 150.0).abs() <= 3.0 * g.pitch, "CD {cd}");
    }

    #[test]
    fn proc03_etch_removes_only_exposed_oxide_and_keeps_masked_oxide() {
        let s = run_chain(&tile(), &win(), &defaults()[..4], 1).unwrap();
        let w = &s[3].wafer;
        let g = ColumnGrid::sample(w, &win()).unwrap();
        let row = g.ny / 2;
        let line = g.idx(g.nx / 2, row);
        let space = g.idx(g.nx / 2 + (160.0 / g.pitch) as usize, row);
        assert!(g.cols[line].iter().any(|r| r.material == Material::SiliconDioxide));
        assert!(!g.cols[space].iter().any(|r| r.material == Material::SiliconDioxide));
        assert!(s[3].outcome.metric("Maximum etch depth").unwrap() >= 90.0);
    }

    #[test]
    fn data01_invalid_and_non_finite_inputs_never_produce_a_result() {
        let mut d = StageParams::default_for(Stage::Deposition);
        if let StageParams::Deposition(x) = &mut d { x.thickness_nm = f64::NAN; }
        assert!(matches!(evaluate_stage(&tile(), &win(), &d, 1, 1.0), Err(PipelineError::Invalid(_))));
        let etch = StageParams::default_for(Stage::Etching);
        let r = evaluate_stage(&tile(), &win(), &etch, 1, 1.0);
        assert!(matches!(r, Err(PipelineError::MissingInput(_)) | Ok(_)));
        if let Ok(ev) = r { assert!(ev.outcome.metric("Maximum etch depth").unwrap() < 10.0, "bare silicon in oxide chemistry is effectively an etch stop"); }
        assert!(evaluate_stage(&tile(), &win(), &d, 1, f64::NAN).is_err());
    }

    #[test]
    fn implant_is_blocked_by_thick_masks_and_requires_exposed_silicon() {
        let mut d = StageParams::default_for(Stage::Doping);
        if let StageParams::Doping(x) = &mut d { x.energy_kev = 10.0; }
        let blocked = run_chain(&tile(), &win(), &[StageParams::default_for(Stage::Deposition), StageParams::default_for(Stage::Deposition), d], 1);
        assert!(matches!(blocked, Err((2, PipelineError::MissingInput(_)))), "200 nm of oxide stops a 10 keV boron implant");
    }

    #[test]
    fn anneal_broadening_matches_the_numerical_diffusion_solver() {
        use crate::physics::diffusion::step_1d;
        let x = DopingParams { species: crate::params::Species::Boron, energy_kev: 30.0, dose_cm2: 1e14, anneal_temp_c: 1000.0, anneal_min: 10.0, background_cm3: 1e15, rp_override_nm: None, tilt_deg: 7.0 };
        let m0 = profile_model(&x, 1.0, 0.0);
        let m1 = profile_model(&x, 1.0, 1.0);
        let dz = 2.0; let n = 700;
        let mut c: Vec<f64> = (0..n).map(|i| conc(&m0, (i as f64 + 0.5) * dz)).collect();
        let mats = vec![Material::Silicon; n];
        let dx_cm = dz * 1e-7;
        let dt = 0.4 * dx_cm * dx_cm / m1.diff_cm2s;
        let steps = (600.0 / dt).ceil() as usize;
        let dt = 600.0 / steps as f64;
        for _ in 0..steps { c = step_1d(&c, &mats, dx_cm, m1.diff_cm2s, dt).unwrap(); }
        let analytic = conc(&m1, m1.rp);
        let numeric = c[(m1.rp / dz) as usize];
        assert!((numeric / analytic - 1.0).abs() < 0.06, "numeric {numeric:e} vs analytic {analytic:e}");
    }

    #[test]
    fn range_table_interpolation_is_monotonic_and_hits_table_points() {
        let t = crate::params::Species::Boron.range_table();
        assert!((interp_range(t, 30.0).0 - 100.0).abs() < 1e-9);
        assert!(interp_range(t, 45.0).0 > 100.0 && interp_range(t, 45.0).0 < 200.0);
        assert_eq!(dopant_bucket(5e16), 0); assert_eq!(dopant_bucket(1e17), 1); assert_eq!(dopant_bucket(3e19), 3); assert_eq!(dopant_bucket(1e25), 6);
    }
}
