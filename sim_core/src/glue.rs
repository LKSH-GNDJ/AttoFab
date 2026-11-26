// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The seam between the pipeline and the physics.
//!
//! **This file compiles** against the application's real source: `tools/check_patches.sh`
//! builds `sim_core` with this module declared. The first time that was done it
//! failed, because it passed the application's `Material` to the bridge's `material_code`,
//! which takes the bridge's own type of the same name. The conversion is explicit
//! now, and the build is the evidence.
//!
//! The logic it calls lives in `atto_bridge`, which *is* tested — 16 tests,
//! covering the implant column, the step height, the junction crossing, the
//! oxidation split, the etch rates and the grid round trip. What is unverified
//! here is only the plumbing.

use crate::geom::{Column, ColumnGrid};
use crate::wafer::{Material, Run, WaferState};

/// the application's `Material` to the code `atto_bridge` uses.
///
/// The two enums have the same variant names but are different types, so the
/// bridge's own `material_code` cannot be called with the application's — that mismatch is why
/// this file did not compile the first time it was built.
pub fn code_of(m: Material) -> u8 {
    match m {
        Material::Air => 0,
        Material::Silicon => 1,
        Material::SiliconDioxide => 2,
        Material::SiliconCarbide => 3,
        Material::GalliumNitride => 4,
        Material::Photoresist => 5,
        Material::Metal => 6,
        Material::Dopant(_) => 8,
    }
}

/// And back: a dopant run becomes silicon, since the application's bucket is not a species.
pub fn material_of(code: u8) -> Material {
    match code {
        0 => Material::Air,
        2 => Material::SiliconDioxide,
        3 => Material::SiliconCarbide,
        4 => Material::GalliumNitride,
        5 => Material::Photoresist,
        6 => Material::Metal,
        _ => Material::Silicon,
    }
}

/// The material conversion, round-tripped. This is the part of the seam that was
/// wrong when the file first compiled: it passed the application's `Material` to the bridge's
/// `material_code`, which takes the bridge's own type of the same name.
#[cfg(test)]
mod conversion_tests {
    use super::*;

    #[test]
    fn every_material_has_a_code_and_comes_back() {
        let all = [
            Material::Air,
            Material::Silicon,
            Material::SiliconDioxide,
            Material::SiliconCarbide,
            Material::GalliumNitride,
            Material::Photoresist,
            Material::Metal,
        ];
        for m in all {
            let c = code_of(m);
            assert_eq!(material_of(c), m, "{m:?} did not survive the round trip");
        }
    }

    #[test]
    fn a_dopant_run_is_code_eight_and_returns_as_silicon() {
        // the application's dopant bucket carries no species, so it cannot come back as a dopant;
        // the bridge counts it as silicon and the species is recorded as dropped.
        assert_eq!(code_of(Material::Dopant(3)), 8);
        assert_eq!(material_of(8), Material::Silicon);
    }
}

/// Flatten the application's grid into the primitives `atto_bridge` accepts.
///
/// This is the fork made concrete: the earlier build keeps its own `ColumnGrid` and hands the
/// contents across as plain numbers, so neither side adopts the other's struct.
pub fn to_bridge_primitives(grid: &ColumnGrid) -> (Vec<bool>, Vec<Vec<(f64, f64, u8)>>) {
    let active = grid.active.clone();
    let cols = grid
        .cols
        .iter()
        .map(|c: &Column| {
            c.iter()
                .map(|r: &Run| (r.z0, r.z1, code_of(r.material)))
                .collect()
        })
        .collect();
    (active, cols)
}

/// Rebuild the application's grid from primitives, for the way back.
pub fn from_bridge_primitives(
    x0: f64,
    y0: f64,
    pitch: f64,
    nx: usize,
    ny: usize,
    z_bottom: f64,
    z_top: f64,
    active: &[bool],
    cols: &[Vec<(f64, f64, u8)>],
) -> Option<ColumnGrid> {
    if cols.len() != nx * ny || active.len() != nx * ny {
        return None;
    }
    let columns: Vec<Column> = cols
        .iter()
        .map(|c| {
            c.iter()
                .map(|(z0, z1, code)| Run {
                    z0: *z0,
                    z1: *z1,
                    material: material_of(*code),
                })
                .collect()
        })
        .collect();
    Some(ColumnGrid {
        x0,
        y0,
        pitch,
        nx,
        ny,
        cols: columns,
        active: active.to_vec(),
        z_bottom,
        z_top,
    })
}

/// Depth of the silicon surface in a column, below any covering film.
pub fn silicon_surface_z(col: &[Run], z_bottom: f64) -> f64 {
    col.iter()
        .rev()
        .find(|r| matches!(r.material, Material::Silicon | Material::Dopant(_)))
        .map(|r| r.z1)
        .unwrap_or(z_bottom)
}

/// Convert a wafer to the primitives the bridge needs.
pub fn wafer_grid(w: &WaferState, grid: &ColumnGrid) -> (Vec<bool>, Vec<Vec<(f64, f64, u8)>>) {
    let _ = w;
    to_bridge_primitives(grid)
}
