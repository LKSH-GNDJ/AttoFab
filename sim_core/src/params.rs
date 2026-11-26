// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Per-stage recipe parameters, validation and nominal durations (PRD §2).

use crate::geom::{EtchChemistry, EtchRates, MaskPattern};
use crate::physics::Stage;
use crate::wafer::Material;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Severity { Warning, Error }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Issue { pub severity: Severity, pub message: String }

impl Issue {
    pub fn error(m: impl Into<String>) -> Self { Self { severity: Severity::Error, message: m.into() } }
    pub fn warn(m: impl Into<String>) -> Self { Self { severity: Severity::Warning, message: m.into() } }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Species { Boron, Phosphorus, Arsenic }

impl Species {
    pub const ALL: [Species; 3] = [Species::Boron, Species::Phosphorus, Species::Arsenic];
    pub fn label(self) -> &'static str { match self { Species::Boron => "Boron", Species::Phosphorus => "Phosphorus", Species::Arsenic => "Arsenic" } }
    /// Intrinsic diffusivity prefactor (cm²/s) and activation energy (eV).
    pub fn diffusion(self) -> (f64, f64) {
        match self { Species::Boron => (0.76, 3.46), Species::Phosphorus => (3.85, 3.66), Species::Arsenic => (22.9, 4.1) }
    }
    /// (energy keV, Rp nm, ΔRp nm) LSS-style tabulated approximations in Si.
    pub fn range_table(self) -> &'static [(f64, f64, f64)] {
        match self {
            Species::Boron => &[(10.0, 35.0, 14.0), (30.0, 100.0, 36.0), (60.0, 200.0, 55.0), (100.0, 310.0, 70.0), (200.0, 580.0, 95.0)],
            Species::Phosphorus => &[(10.0, 14.0, 6.0), (30.0, 40.0, 17.0), (60.0, 80.0, 28.0), (100.0, 134.0, 40.0), (200.0, 270.0, 62.0)],
            Species::Arsenic => &[(10.0, 9.0, 3.5), (30.0, 20.0, 8.0), (60.0, 37.0, 14.0), (100.0, 58.0, 20.0), (200.0, 112.0, 35.0)],
        }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PrepParams { pub rca_dilution: f64, pub ra_nm: f64, pub gamma_sg: f64, pub gamma_sl: f64, pub gamma_lg: f64, pub hf_dip: bool, pub particle_removal_eff: f64 }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct DepositionParams { pub material: Material, pub thickness_nm: f64, pub step_coverage: f64, pub bow_radius_m: f64, pub flux3d: bool, pub sticking: f64, pub source: DepositionSource, pub flux_iterations: u32, pub flux_tolerance: f64 }

/// Where the arriving material comes from.
///
/// `step_coverage` says how much of the requested thickness survives on a sidewall. Sticking
/// is a different quantity: what the surface does with what arrives. They were one number
/// until this step, which is why a sputter deposition could not be described at all.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum DepositionSource {
    /// No preferred direction: a CVD source. Tilt and spread do not apply.
    Cosine,
    /// A beam about a direction, for sputter and evaporation.
    Beam { tilt_deg: f64, spread_deg: f64 },
}

impl Default for DepositionSource {
    fn default() -> Self { DepositionSource::Cosine }
}

impl Default for DepositionParams {
    fn default() -> Self {
        DepositionParams { material: Material::SiliconDioxide, thickness_nm: 100.0, step_coverage: 0.3, bow_radius_m: 40.0, flux3d: false, sticking: 1.0, source: DepositionSource::Cosine, flux_iterations: 1, flux_tolerance: 0.0 }
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct LithoParams {
    pub resist_thickness_nm: f64,
    /// Planarity error the wafer arrives with, nm (P-04). Zero is an ideal flat field.
    ///
    /// The scanner focuses on a plane; a field that has been polished back sits below it, and
    /// the resulting defocus prints as a CD shift. This field is the producer's half of that
    /// coupling — a CMP stage sets it, or a recipe sets it by hand until one exists.
    #[serde(default)]
    pub planarity_error_nm: f64,
    pub mask: MaskPattern,
    pub wavelength_nm: f64,
    pub na: f64,
    pub k1: f64,
    /// Exposure dose divided by the dose-to-clear (positive resist).
    pub dose_ratio: f64,
    pub positive_tone: bool,
    /// Mask-to-wafer alignment error, nm.
    pub align_offset_nm: [f64; 2],
    pub soft_bake_c: f64,
    pub peb_c: f64,
    /// How long the developer acts on the resist, minutes. The animation progress scales it, so a
    /// half-played stage shows a half-developed film, but it no longer *is* the develop time.
    /// Recipes saved before this field existed load with one minute, which is what a finished
    /// stage always used.
    #[serde(default = "default_develop_time_min")]
    pub develop_time_min: f64,
}

fn default_develop_time_min() -> f64 { 1.0 }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct EtchParams {
    pub chemistry: EtchChemistry,
    pub rates: EtchRates,
    pub time_min: f64,
    pub anisotropy: f64,
    pub strip_after: bool,
    /// RF power and chamber pressure, which the earlier build never carried: its rates
    /// were tabulated per chemistry with no operating point.
    #[serde(default = "default_rf_power_w")]
    pub rf_power_w: f64,
    #[serde(default = "default_pressure_mtorr")]
    pub pressure_mtorr: f64,
}

fn default_rf_power_w() -> f64 { 200.0 }
fn default_pressure_mtorr() -> f64 { 20.0 }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DopingParams {
    pub species: Species,
    pub energy_kev: f64,
    pub dose_cm2: f64,
    pub anneal_temp_c: f64,
    pub anneal_min: f64,
    pub background_cm3: f64,
    pub rp_override_nm: Option<f64>,
    /// Implant tilt from the <100> axis, degrees. 7 is the standard
    /// off-axis cut; 0 channels and gives a far deeper tail.
    #[serde(default = "default_tilt_deg")]
    pub tilt_deg: f64,
}

fn default_tilt_deg() -> f64 { 7.0 }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CleanParams { pub strip_efficiency: f64, pub chemical_clean: bool, pub particle_removal_eff: f64 }

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum Quantity {
    FilmThickness(Material),
    LineWidth(Material),
    StepHeight,
    JunctionDepth,
    TopographyRms,
    FilmStress,
    ResistResidue,
    /// Peak von Mises stress in the film(s) from the tile cross-section FEM.
    PeakFilmStress,
}

impl Quantity {
    pub fn label(&self) -> String {
        match self {
            Quantity::FilmThickness(m) => format!("{} thickness", mat_name(*m)),
            Quantity::LineWidth(m) => format!("{} critical dimension", mat_name(*m)),
            Quantity::StepHeight => "Etch depth / step height".into(),
            Quantity::JunctionDepth => "Dopant depth (≥1e17 cm⁻³)".into(),
            Quantity::TopographyRms => "Surface roughness (RMS)".into(),
            Quantity::FilmStress => "Film stress".into(),
            Quantity::ResistResidue => "Resist residue".into(),
            Quantity::PeakFilmStress => "Peak film stress (FEM von Mises)".into(),
        }
    }
    pub fn unit(&self) -> &'static str { match self { Quantity::FilmStress | Quantity::PeakFilmStress => "MPa", _ => "nm" } }
}

pub fn mat_name(m: Material) -> &'static str {
    match m {
        Material::Air => "Air", Material::Silicon => "Silicon", Material::SiliconDioxide => "SiO₂", Material::SiliconCarbide => "SiC",
        Material::GalliumNitride => "GaN", Material::Photoresist => "Photoresist", Material::Metal => "Metal", Material::Dopant(_) => "Dopant",
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Target { pub quantity: Quantity, pub nominal: f64, pub tolerance: f64 }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetrologyParams { pub targets: Vec<Target>, pub sites: usize, pub noise_nm: f64, pub wavelength_nm: f64, pub incidence_deg: f64 }

/// Chemical-mechanical planarisation (P-01..P-03).
///
/// The polish removes the top film at the Preston rate `k p v t`. The recipe carries the
/// pressure, the pad velocity and the time; the model carries the rate constants.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CmpParams {
    /// What the pad is polishing.
    pub material: Material,
    /// Downward pressure, kPa.
    pub pressure_kpa: f64,
    /// Pad velocity relative to the wafer, m/s.
    pub velocity_mps: f64,
    /// Polish time, minutes.
    pub time_min: f64,
    /// Fraction of the time spent after the film clears, for P-03's stop layer.
    pub over_polish: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum StageParams {
    SubstratePrep(PrepParams),
    Deposition(DepositionParams),
    Lithography(LithoParams),
    Etching(EtchParams),
    Doping(DopingParams),
    FinalClean(CleanParams),
    Metrology(MetrologyParams),
    Cmp(CmpParams),
}

impl StageParams {
    pub fn stage(&self) -> Stage {
        match self {
            StageParams::SubstratePrep(_) => Stage::SubstratePrep, StageParams::Deposition(_) => Stage::Deposition,
            StageParams::Lithography(_) => Stage::Lithography, StageParams::Etching(_) => Stage::Etching,
            StageParams::Doping(_) => Stage::Doping, StageParams::FinalClean(_) => Stage::FinalClean, StageParams::Metrology(_) => Stage::Metrology, StageParams::Cmp(_) => Stage::Cmp,
        }
    }

    pub fn default_for(stage: Stage) -> Self {
        match stage {
            Stage::SubstratePrep => StageParams::SubstratePrep(PrepParams { rca_dilution: 5.0, ra_nm: 0.2, gamma_sg: 60.0, gamma_sl: 20.0, gamma_lg: 72.8, hf_dip: true, particle_removal_eff: 0.995 }),
            Stage::Deposition => StageParams::Deposition(DepositionParams { material: Material::SiliconDioxide, thickness_nm: 100.0, step_coverage: 0.3, bow_radius_m: 40.0, flux3d: false, sticking: 1.0, source: DepositionSource::Cosine, flux_iterations: 1, flux_tolerance: 0.0 }),
            Stage::Lithography => StageParams::Lithography(LithoParams {
                planarity_error_nm: 0.0,
                resist_thickness_nm: 150.0,
                mask: MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, length_nm: 1000.0, vertical: true, clear_field: true },
                wavelength_nm: 193.0, na: 0.85, k1: 0.4, dose_ratio: 2.0, positive_tone: true, align_offset_nm: [0.0, 0.0], soft_bake_c: 110.0, peb_c: 110.0, develop_time_min: 0.5,
            }),
            Stage::Etching => StageParams::Etching(EtchParams { chemistry: EtchChemistry::OxideRie, rates: EtchChemistry::OxideRie.default_rates(), time_min: 2.5, anisotropy: 0.95, strip_after: false, rf_power_w: 200.0, pressure_mtorr: 20.0 }),
            Stage::Doping => StageParams::Doping(DopingParams { species: Species::Boron, energy_kev: 30.0, dose_cm2: 1.0e14, anneal_temp_c: 1000.0, anneal_min: 10.0, background_cm3: 1.0e15, rp_override_nm: None, tilt_deg: 7.0 }),
            Stage::FinalClean => StageParams::FinalClean(CleanParams { strip_efficiency: 1.0, chemical_clean: true, particle_removal_eff: 0.99 }),
            Stage::Cmp => StageParams::Cmp(CmpParams { material: Material::SiliconDioxide, pressure_kpa: 20.0, velocity_mps: 0.5, time_min: 1.0, over_polish: 0.0 }),
            Stage::Metrology => StageParams::Metrology(MetrologyParams {
                targets: vec![
                    Target { quantity: Quantity::FilmThickness(Material::SiliconDioxide), nominal: 100.0, tolerance: 10.0 },
                    Target { quantity: Quantity::LineWidth(Material::SiliconDioxide), nominal: 150.0, tolerance: 40.0 },
                    Target { quantity: Quantity::StepHeight, nominal: 100.0, tolerance: 20.0 },
                    Target { quantity: Quantity::JunctionDepth, nominal: 250.0, tolerance: 100.0 },
                ],
                sites: 5, noise_nm: 0.0, wavelength_nm: 633.0, incidence_deg: 70.0,
            }),
        }
    }

    /// Parameter validation. Any `Error` issue blocks the stage from producing
    /// a committed result.
    pub fn validate(&self) -> Vec<Issue> {
        let mut v = Vec::new();
        let finite = |name: &str, x: f64, min: f64, max: f64, v: &mut Vec<Issue>| {
            if !x.is_finite() { v.push(Issue::error(format!("{name} must be a finite number"))); }
            else if x < min || x > max { v.push(Issue::error(format!("{name} = {x} is outside [{min}, {max}]"))); }
        };
        match self {
            StageParams::SubstratePrep(p) => {
                finite("RCA dilution", p.rca_dilution, 0.0, 1.0e3, &mut v); finite("roughness Ra", p.ra_nm, 0.0, 1.0e3, &mut v);
                finite("γ liquid-gas", p.gamma_lg, 1.0e-3, 1.0e3, &mut v); finite("γ solid-gas", p.gamma_sg, 0.0, 1.0e3, &mut v); finite("γ solid-liquid", p.gamma_sl, 0.0, 1.0e3, &mut v);
                finite("particle removal efficiency", p.particle_removal_eff, 0.0, 1.0, &mut v);
            }
            StageParams::Deposition(p) => {
                finite("film thickness", p.thickness_nm, 0.0, 1.0e6, &mut v); finite("step coverage", p.step_coverage, 0.0, 1.0, &mut v); finite("bow radius", p.bow_radius_m.abs(), 1.0e-3, 1.0e9, &mut v);
                if matches!(p.material, Material::Air | Material::Photoresist | Material::Dopant(_)) { v.push(Issue::error("deposition material must be a solid film material (SiO₂, SiC, GaN, metal or silicon)")); }
                if p.thickness_nm == 0.0 { v.push(Issue::error("film thickness must be > 0")); }
                finite("sticking", p.sticking, 0.0, 1.0, &mut v);
                if p.flux_iterations == 0 || p.flux_iterations > 64 {
                    v.push(Issue::error("flux iterations must be within 1..=64"));
                }
                finite("flux tolerance", p.flux_tolerance, 0.0, 1.0, &mut v);
                match p.source {
                    DepositionSource::Cosine => {}
                    DepositionSource::Beam { tilt_deg, spread_deg } => {
                        finite("beam tilt", tilt_deg.abs(), 0.0, 89.0, &mut v);
                        finite("beam spread", spread_deg, 0.0, 90.0, &mut v);
                    }
                }
            }
            StageParams::Lithography(p) => {
                finite("resist thickness", p.resist_thickness_nm, 1.0, 1.0e5, &mut v); finite("wavelength", p.wavelength_nm, 1.0, 1.0e4, &mut v);
                finite("NA", p.na, 0.01, 2.0, &mut v); finite("k1", p.k1, 0.0, 2.0, &mut v); finite("dose ratio", p.dose_ratio, 0.01, 100.0, &mut v);
                finite("alignment offset X", p.align_offset_nm[0], -1.0e5, 1.0e5, &mut v); finite("alignment offset Y", p.align_offset_nm[1], -1.0e5, 1.0e5, &mut v);
                finite("soft bake", p.soft_bake_c, 20.0, 300.0, &mut v); finite("post-exposure bake", p.peb_c, 20.0, 300.0, &mut v); finite("develop time", p.develop_time_min, 0.01, 60.0, &mut v);
                if let Err(e) = p.mask.validate() { v.push(Issue::error(e.to_string())); }
            }
            StageParams::Etching(p) => {
                finite("etch time", p.time_min, 0.0, 1.0e4, &mut v); finite("anisotropy", p.anisotropy, 0.0, 1.0, &mut v);
                for (n, r) in [("Si", p.rates.si), ("SiO₂", p.rates.sio2), ("SiC", p.rates.sic), ("GaN", p.rates.gan), ("resist", p.rates.resist), ("metal", p.rates.metal)] { finite(&format!("{n} etch rate"), r, 0.0, 1.0e5, &mut v); }
                if p.time_min == 0.0 { v.push(Issue::error("etch time must be > 0")); }
            }
            StageParams::Doping(p) => {
                finite("implant energy", p.energy_kev, 1.0, 1.0e4, &mut v); finite("dose", p.dose_cm2, 1.0e8, 1.0e18, &mut v);
                finite("anneal temperature", p.anneal_temp_c, 25.0, 1400.0, &mut v); finite("anneal time", p.anneal_min, 0.0, 1.0e4, &mut v); finite("background doping", p.background_cm3, 1.0e10, 1.0e21, &mut v);
                if let Some(r) = p.rp_override_nm { finite("projected range", r, 0.1, 1.0e5, &mut v); }
            }
            StageParams::FinalClean(p) => { finite("strip efficiency", p.strip_efficiency, 0.0, 1.0, &mut v); finite("particle removal efficiency", p.particle_removal_eff, 0.0, 1.0, &mut v); }
            StageParams::Cmp(p) => {
                finite("pressure", p.pressure_kpa, 0.1, 200.0, &mut v);
                finite("pad velocity", p.velocity_mps, 0.01, 5.0, &mut v);
                finite("polish time", p.time_min, 0.01, 60.0, &mut v);
                finite("over-polish", p.over_polish, 0.0, 1.0, &mut v);
                if matches!(p.material, Material::Air | Material::Photoresist | Material::Dopant(_)) {
                    v.push(Issue::error("a polish needs a solid film to remove, not air, resist or a dopant"));
                }
            }
            StageParams::Metrology(p) => {
                finite("sampling sites", p.sites as f64, 1.0, 15.0, &mut v); finite("noise", p.noise_nm, 0.0, 1.0e3, &mut v);
                finite("ellipsometer wavelength", p.wavelength_nm, 100.0, 2.0e3, &mut v); finite("incidence angle", p.incidence_deg, 1.0, 89.0, &mut v);
                for t in &p.targets { if !t.nominal.is_finite() || !t.tolerance.is_finite() || t.tolerance < 0.0 { v.push(Issue::error(format!("target '{}' needs a finite nominal and tolerance >= 0", t.quantity.label()))); } }
            }
        }
        v
    }

    /// Nominal simulated process time for this stage, seconds.
    pub fn duration_s(&self) -> f64 {
        match self {
            StageParams::SubstratePrep(_) => 900.0,
            StageParams::Deposition(p) => (p.thickness_nm / 20.0 * 60.0).max(30.0),
            StageParams::Lithography(_) => 600.0,
            StageParams::Etching(p) => p.time_min * 60.0 + if p.strip_after { 120.0 } else { 0.0 },
            StageParams::Doping(p) => 120.0 + p.anneal_min * 60.0,
            StageParams::FinalClean(_) => 420.0,
            StageParams::Cmp(_) => 120.0,
            StageParams::Metrology(_) => 300.0,
        }
    }
}
