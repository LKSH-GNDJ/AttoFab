// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Adaptive refinement: spend the budget where the error is (PRD W-06).
//!
//! Uniform refinement costs the cube of the resolution for a linear gain. What
//! matters is placing fine leaves where the geometry is changing — an edge, a
//! corner, a surface that steps — and leaving flat regions coarse.
//!
//! This module is the *placement decision*: a per-cell error estimate and a
//! budget-constrained selection. It does not subdivide the octree, because doing
//! so without a compiler in the loop is how an adaptive scheme becomes a
//! debugging session.

use crate::column::ColumnGrid;
use crate::svo::Material;

/// Per-cell error estimate: how much the geometry changes across a cell.
///
/// The estimate is the largest difference in top-of-solid height between a cell
/// and its four neighbours, which is large exactly at an edge and zero on a flat
/// surface. Air columns take the largest neighbour difference, because the edge of
/// a feature is where a coarse grid is most wrong.
pub fn error_estimate(grid: &ColumnGrid) -> Vec<f64> {
    let n = grid.cols.len();
    let mut out = vec![0.0; n];
    if grid.nx < 2 || grid.ny < 2 {
        return out;
    }
    let top = |k: usize| -> Option<f64> {
        grid.cols[k].iter().rev().find(|r| r.material != Material::Air).map(|r| r.z1)
    };
    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            if !grid.active[k] {
                continue;
            }
            let here = top(k);
            let mut worst: f64 = 0.0;
            for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                let (jx, jy) = (ix as i32 + dx, iy as i32 + dy);
                if jx < 0 || jy < 0 || jx >= grid.nx as i32 || jy >= grid.ny as i32 {
                    continue;
                }
                let j = grid.idx(jx as usize, jy as usize);
                if !grid.active[j] {
                    continue;
                }
                let d = match (here, top(j)) {
                    (Some(a), Some(b)) => (a - b).abs(),
                    // a cell next to air is on an edge
                    _ => grid.pitch,
                };
                worst = worst.max(d);
            }
            out[k] = worst;
        }
    }
    out
}

/// The cells to refine, highest error first, at most `budget` of them.
///
/// Ties break by index so the result is deterministic: an adaptive scheme that
/// refines differently on two identical runs cannot be regression-tested.
pub fn select(errors: &[f64], budget: usize) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..errors.len()).filter(|i| errors[*i] > 0.0).collect();
    idx.sort_by(|a, b| {
        errors[*b]
            .partial_cmp(&errors[*a])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.cmp(b))
    });
    idx.truncate(budget);
    idx.sort_unstable();
    idx
}

/// A refinement plan: which cells to subdivide and what the total error they hold.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    pub cells: Vec<usize>,
    /// Sum of the error estimates of the selected cells.
    pub captured: f64,
    /// Sum over every active cell, for a coverage fraction.
    pub total: f64,
}

impl Plan {
    /// The fraction of the total error the plan captures.
    pub fn coverage(&self) -> f64 {
        if self.total > 0.0 { self.captured / self.total } else { 1.0 }
    }
}

/// Plan a refinement under a budget.
pub fn plan(grid: &ColumnGrid, budget: usize) -> Plan {
    let errors = error_estimate(grid);
    let cells = select(&errors, budget);
    let captured = cells.iter().map(|i| errors[*i]).sum();
    let total = errors.iter().filter(|e| **e > 0.0).sum();
    Plan { cells, captured, total }
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Nm, Um};
    use super::*;
    use crate::column::{set_interval, ProcessWindow};
    use crate::svo::Wafer;

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        ColumnGrid::sample(&w, &win).unwrap()
    }

    #[test]
    fn a_flat_surface_has_no_error_anywhere() {
        let g = grid();
        let e = error_estimate(&g);
        assert!(e.iter().all(|v| *v == 0.0), "a flat wafer needs no refinement");
        assert!(plan(&g, 10).cells.is_empty());
    }

    #[test]
    fn an_edge_carries_the_error() {
        let mut g = grid();
        // raise a block on the left half only: the boundary column is an edge
        for iy in 0..g.ny {
            for ix in 0..g.nx / 2 {
                let k = g.idx(ix, iy);
                let top = g.cols[k].iter().rev().find(|r| r.material != Material::Air).map(|r| r.z1).unwrap_or(g.z_bottom);
                set_interval(&mut g.cols[k], top, top + 100.0, Material::SiliconDioxide);
            }
        }
        let e = error_estimate(&g);
        let at_edge = g.idx(g.nx / 2 - 1, g.ny / 2);
        let in_flat = g.idx(1, g.ny / 2);
        assert!(e[at_edge] > 0.0, "the boundary column has error");
        assert!(e[at_edge] > e[in_flat], "and more than the flat interior");
    }

    #[test]
    fn the_budget_selects_the_worst_cells() {
        let errors = vec![0.0, 5.0, 1.0, 9.0, 3.0, 0.0];
        let picked = select(&errors, 2);
        assert_eq!(picked, vec![1, 3], "the two highest, in index order");
        assert_eq!(select(&errors, 100).len(), 4, "never more than the non-zero cells");
        assert!(select(&errors, 0).is_empty());
    }

    #[test]
    fn ties_break_deterministically() {
        let errors = vec![2.0, 2.0, 2.0, 2.0];
        let a = select(&errors, 2);
        let b = select(&errors, 2);
        assert_eq!(a, b, "the same input refines the same cells");
        assert_eq!(a, vec![0, 1], "and the lowest indices win");
    }

    #[test]
    fn a_small_budget_captures_most_of_the_error_on_a_blocky_surface() {
        let mut g = grid();
        for iy in 0..g.ny {
            for ix in 0..g.nx / 2 {
                let k = g.idx(ix, iy);
                let top = g.cols[k].iter().rev().find(|r| r.material != Material::Air).map(|r| r.z1).unwrap_or(g.z_bottom);
                set_interval(&mut g.cols[k], top, top + 100.0, Material::SiliconDioxide);
            }
        }
        let p = plan(&g, 200);
        assert!(!p.cells.is_empty());
        assert!(p.coverage() > 0.5, "a fifth of the cells should capture most of the error: {:.2}", p.coverage());
        assert!(p.coverage() <= 1.0);
    }
}
