// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The stack tab: the column under the cursor, run by run (PRD G-12).
//!
//! The inspector shows what is actually at a point, in order, with thicknesses —
//! not a legend. A run that is present but buried is listed like any other, which
//! is the whole reason to show a stack rather than a summary.

use atto_units::{Nm};
use atto_geom::column::ColumnGrid;
use atto_geom::svo::Material;

#[derive(Debug, Clone, PartialEq)]
/// One row of the stack view: a material, its thickness and where it came from.
pub struct StackRow {
    pub material: String,
    pub z_bottom_nm: f64,
    pub z_top_nm: f64,
    pub thickness_nm: f64,
    /// False for air: a gap is part of the stack and must be shown, not skipped.
    pub solid: bool,
}

fn material_name(m: Material) -> &'static str {
    match m {
        Material::Air => "Air",
        Material::Silicon => "Silicon",
        Material::SiliconDioxide => "SiO2",
        Material::SiliconNitride => "Si3N4",
        Material::Photoresist => "Photoresist",
        Material::Metal => "Metal",
    }
}

/// The stack at a point, top-down, with air runs included.
pub fn stack_rows(grid: &ColumnGrid, x: Nm, y_nm: f64) -> Option<Vec<StackRow>> {
    let x_nm = x.value();
    if !(grid.pitch > 0.0) {
        return None;
    }
    let ix = ((x_nm - grid.x0) / grid.pitch).floor();
    let iy = ((y_nm - grid.y0) / grid.pitch).floor();
    if ix < 0.0 || iy < 0.0 || ix >= grid.nx as f64 || iy >= grid.ny as f64 {
        return None;
    }
    let k = grid.idx(ix as usize, iy as usize);
    if !grid.active[k] {
        return None;
    }
    let mut rows: Vec<StackRow> = grid.cols[k]
        .iter()
        .filter(|r| r.z1 - r.z0 > 1e-9)
        .map(|r| StackRow {
            material: material_name(r.material).to_string(),
            z_bottom_nm: r.z0,
            z_top_nm: r.z1,
            thickness_nm: r.z1 - r.z0,
            solid: r.material != Material::Air,
        })
        .collect();
    rows.reverse();
    Some(rows)
}

/// Total solid thickness at a point, nm.
pub fn solid_thickness_nm(rows: &[StackRow]) -> f64 {
    rows.iter().filter(|r| r.solid).map(|r| r.thickness_nm).sum()
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;
    use atto_geom::column::{set_interval, ProcessWindow};
    use atto_geom::svo::Wafer;

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        ColumnGrid::sample(&w, &win).unwrap()
    }

    #[test]
    fn the_stack_is_listed_top_down_with_thicknesses() {
        let mut g = grid();
        let k = g.idx(5, 5);
        let top = atto_geom::column::solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
        set_interval(&mut g.cols[k], top, top + 100.0, Material::SiliconDioxide);
        set_interval(&mut g.cols[k], top + 100.0, top + 180.0, Material::Photoresist);
        let x = g.x0 + 5.5 * g.pitch;
        let y = g.y0 + 5.5 * g.pitch;
        let rows = stack_rows(&g, Nm::new(x), y).expect("a column");
        // the column includes the air headroom above the stack, so the first
        // solid run is what a reader is looking for
        let solid: Vec<&StackRow> = rows.iter().filter(|r| r.solid).collect();
        assert_eq!(solid[0].material, "Photoresist", "topmost solid first");
        assert!((solid[0].thickness_nm - 80.0).abs() < 1e-6);
        assert_eq!(solid[1].material, "SiO2");
        assert!((solid[1].thickness_nm - 100.0).abs() < 1e-6);
        assert!(solid.len() >= 3, "the substrate is part of the stack");
        assert!(rows[0].material == "Air" && !rows[0].solid, "headroom is shown as air");
    }

    #[test]
    fn a_gap_is_shown_not_skipped() {
        let mut g = grid();
        let k = g.idx(7, 7);
        let top = atto_geom::column::solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
        set_interval(&mut g.cols[k], top, top + 40.0, Material::SiliconDioxide);
        set_interval(&mut g.cols[k], top + 40.0, top + 60.0, Material::Air);
        set_interval(&mut g.cols[k], top + 60.0, top + 100.0, Material::SiliconDioxide);
        let x = g.x0 + 7.5 * g.pitch;
        let y = g.y0 + 7.5 * g.pitch;
        let rows = stack_rows(&g, Nm::new(x), y).expect("a column");
        assert!(rows.iter().any(|r| r.material == "Air" && !r.solid), "the void is listed");
        let solid = solid_thickness_nm(&rows);
        let total: f64 = rows.iter().map(|r| r.thickness_nm).sum();
        assert!(solid < total, "the gap is excluded from solid thickness");
    }

    #[test]
    fn a_point_outside_the_grid_has_no_stack() {
        let g = grid();
        assert!(stack_rows(&g, Nm::new(-10_000.0), 0.0).is_none());
        assert!(stack_rows(&g, Nm::new(1.0e9), 0.0).is_none());
    }
}
