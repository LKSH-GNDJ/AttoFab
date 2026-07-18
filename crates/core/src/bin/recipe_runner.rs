//! recipe_runner: a small CLI binary that reads a JSON process recipe from
//! stdin, plays it back against a `Wafer2d`, and writes the resulting
//! wafer state as JSON to stdout.
//!
//! This is the integration point the Python backend (`backend/`) shells
//! out to - deliberately a subprocess boundary rather than a native PyO3
//! extension, so the backend has no compiled-extension/toolchain-version
//! coupling to the Rust engine: any Python environment that can spawn a
//! process and read JSON can drive the real physics engine.
//!
//! Recipe format:
//! ```text
//! {
//!   "nx": 80, "ny": 100, "dx_um": 0.01, "dy_um": 0.01,
//!   "substrate": { "dopant": "Boron", "concentration_cm3": 1e15 },
//!   "steps": [
//!     { "op": "spin_photoresist", "thickness_voxels": 10 },
//!     { "op": "expose", "mask": { "range": [0, 27] } },
//!     { "op": "develop" },
//!     { "op": "oxidize", "temperature_c": 1000, "time_hours": 2, "ambient": "Wet",
//!       "mask": { "range": [0, 27], "invert": true } },
//!     { "op": "implant", "dopant": "Phosphorus", "dose_cm2": 1e15, "energy_kev": 60 },
//!     { "op": "anneal", "temperature_c": 1000, "time_minutes": 20 },
//!     { "op": "etch_anisotropic", "depth_voxels": 12 },
//!     { "op": "deposit", "material": "Metal", "thickness_voxels": 12 }
//!   ]
//! }
//! ```
//! `mask` accepts either an explicit `[true, false, ...]` array (length nx)
//! or a `{ "range": [lo, hi], "invert": bool }` shorthand.

use attofab_core::etch::{deposit, etch_anisotropic, etch_isotropic};
use attofab_core::grid2d::Wafer2d;
use attofab_core::lithography::{develop, expose};
use attofab_core::materials::{Ambient, Dopant, Material};
use serde::Deserialize;
use std::io::{self, Read, Write};
use std::process::ExitCode;

#[derive(Deserialize)]
struct RecipeFile {
    nx: usize,
    ny: usize,
    dx_um: f64,
    dy_um: f64,
    substrate: SubstrateSpec,
    steps: Vec<StepSpec>,
}

#[derive(Deserialize)]
struct SubstrateSpec {
    dopant: Dopant,
    concentration_cm3: f64,
}

#[derive(Deserialize)]
#[serde(untagged)]
enum MaskSpec {
    Explicit(Vec<bool>),
    Range {
        range: [usize; 2],
        #[serde(default)]
        invert: bool,
    },
}

fn resolve_mask(spec: &MaskSpec, nx: usize) -> Vec<bool> {
    match spec {
        MaskSpec::Explicit(v) => v.clone(),
        MaskSpec::Range { range: [lo, hi], invert } => (0..nx)
            .map(|x| {
                let inside = x >= *lo && x < *hi;
                if *invert {
                    !inside
                } else {
                    inside
                }
            })
            .collect(),
    }
}

#[derive(Deserialize)]
#[serde(tag = "op", rename_all = "snake_case")]
enum StepSpec {
    SpinPhotoresist { thickness_voxels: usize },
    Expose { mask: MaskSpec },
    Develop,
    Oxidize {
        temperature_c: f64,
        time_hours: f64,
        ambient: Ambient,
        mask: Option<MaskSpec>,
    },
    Implant {
        dopant: Dopant,
        dose_cm2: f64,
        energy_kev: f64,
        mask: Option<MaskSpec>,
    },
    Anneal {
        temperature_c: f64,
        time_minutes: f64,
        #[serde(default = "default_safety_factor")]
        safety_factor: f64,
    },
    EtchAnisotropic { depth_voxels: usize },
    EtchIsotropic { iterations: usize },
    Deposit { material: Material, thickness_voxels: usize },
}

fn default_safety_factor() -> f64 {
    0.2
}

fn run() -> Result<(), String> {
    let mut input = String::new();
    io::stdin()
        .read_to_string(&mut input)
        .map_err(|e| format!("failed to read stdin: {e}"))?;

    let recipe: RecipeFile =
        serde_json::from_str(&input).map_err(|e| format!("invalid recipe JSON: {e}"))?;

    let mut wafer = Wafer2d::new(
        recipe.nx,
        recipe.ny,
        recipe.dx_um,
        recipe.dy_um,
        recipe.substrate.dopant,
        recipe.substrate.concentration_cm3,
    );

    for (i, step) in recipe.steps.iter().enumerate() {
        match step {
            StepSpec::SpinPhotoresist { thickness_voxels } => {
                wafer.spin_photoresist(*thickness_voxels);
            }
            StepSpec::Expose { mask } => {
                let m = resolve_mask(mask, recipe.nx);
                expose(&mut wafer, &m);
            }
            StepSpec::Develop => {
                develop(&mut wafer);
            }
            StepSpec::Oxidize { temperature_c, time_hours, ambient, mask } => {
                let m = mask.as_ref().map(|s| resolve_mask(s, recipe.nx));
                wafer.oxidize(*temperature_c, *time_hours, *ambient, m.as_deref());
            }
            StepSpec::Implant { dopant, dose_cm2, energy_kev, mask } => {
                let m = mask.as_ref().map(|s| resolve_mask(s, recipe.nx));
                wafer.implant(*dopant, *dose_cm2, *energy_kev, m.as_deref());
            }
            StepSpec::Anneal { temperature_c, time_minutes, safety_factor } => {
                wafer.anneal(*temperature_c, *time_minutes, *safety_factor);
            }
            StepSpec::EtchAnisotropic { depth_voxels } => {
                etch_anisotropic(&mut wafer, *depth_voxels);
            }
            StepSpec::EtchIsotropic { iterations } => {
                etch_isotropic(&mut wafer, *iterations);
            }
            StepSpec::Deposit { material, thickness_voxels } => {
                deposit(&mut wafer, *material, *thickness_voxels);
            }
        }
        eprintln!("step {i} ok: {}", wafer.process_log.last().cloned().unwrap_or_default());
    }

    let json = serde_json::to_string(&wafer).map_err(|e| format!("failed to serialize wafer: {e}"))?;
    io::stdout()
        .write_all(json.as_bytes())
        .map_err(|e| format!("failed to write stdout: {e}"))?;
    Ok(())
}

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("recipe_runner error: {e}");
            ExitCode::FAILURE
        }
    }
}
