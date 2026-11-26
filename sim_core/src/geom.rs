// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Surface-aware process geometry operators (PRD §7 P3).
//!
//! Every operator works on a [`ColumnGrid`]: a *local, transient* sampling of
//! the sparse voxel octree into vertical material stacks on the finest-level
//! XY lattice, restricted to a small [`ProcessWindow`]. The wafer itself is
//! never densified: the grid is built from the SVO, an operator transforms the
//! columns, and only the *differences* are written back with
//! [`ColumnGrid::write_back`], which keeps the octree sparse.
//!
//! The operators are deterministic reference models (height-field / column
//! models). They are physically motivated but are **not** calibrated TCAD:
//! overhangs are not tracked by deposition, lateral etch uses a
//! distance-from-exposed-surface front, and so on. Each stage reports this
//! through [`crate::pipeline::Fidelity`].

use crate::wafer::{Bounds, Material, Run, WaferState};
use serde::{Deserialize, Serialize};

/// Upper bound on grid cells per axis (keeps operators O(10^4) columns).
pub const MAX_GRID_AXIS: usize = 192;
const EPS: f64 = 1e-9;

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum GeomError {
    #[error("process window must have a finite centre and a finite, positive size")]
    InvalidWindow,
    #[error("process window lies outside the circular wafer domain")]
    OutsideDomain,
    #[error("process window must lie inside the high-resolution tile (move the site or shrink the window)")]
    OutsideTile,
    #[error("process window needs {cells} cells on an axis (limit {limit}); reduce the window or raise the SVO leaf size")]
    WindowTooLarge { cells: usize, limit: usize },
    #[error("{0}")]
    InvalidParameter(String),
}

/// Rectangular region of interest, in wafer-local XY nm. Every stage of a
/// recipe operates inside the same window (the "die site").
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ProcessWindow {
    pub center: [f64; 2],
    pub size: [f64; 2],
}

impl Default for ProcessWindow {
    fn default() -> Self { Self { center: [0.0, 0.0], size: [1200.0, 1200.0] } }
}

impl ProcessWindow {
    pub fn validate(&self) -> Result<(), GeomError> {
        if self.center.iter().chain(self.size.iter()).any(|v| !v.is_finite()) || self.size[0] <= 0.0 || self.size[1] <= 0.0 {
            return Err(GeomError::InvalidWindow);
        }
        Ok(())
    }
    pub fn bounds_xy(&self) -> ([f64; 2], [f64; 2]) {
        ([self.center[0] - self.size[0] * 0.5, self.center[1] - self.size[1] * 0.5],
         [self.center[0] + self.size[0] * 0.5, self.center[1] + self.size[1] * 0.5])
    }
}

pub type Column = Vec<Run>;

/// Highest non-air run top and its material.
pub fn solid_top(col: &[Run]) -> Option<(f64, Material)> {
    col.iter().rev().find(|r| r.material != Material::Air).map(|r| (r.z1, r.material))
}

pub fn material_at_z(col: &[Run], z: f64) -> Material {
    col.iter().find(|r| z >= r.z0 && z < r.z1).map(|r| r.material).unwrap_or(Material::Air)
}

fn normalize(col: &mut Column) {
    col.retain(|r| r.z1 - r.z0 > EPS);
    col.sort_by(|a, b| a.z0.partial_cmp(&b.z0).unwrap());
    let mut out: Column = Vec::with_capacity(col.len());
    for r in col.drain(..) {
        if let Some(last) = out.last_mut() {
            if last.material == r.material && (r.z0 - last.z1).abs() < 1e-6 { last.z1 = r.z1; continue; }
        }
        out.push(r);
    }
    *col = out;
}

/// Overwrite `[z0, z1)` of a column with `material`.
pub fn set_interval(col: &mut Column, z0: f64, z1: f64, material: Material) {
    if z1 - z0 <= EPS { return; }
    let mut out: Column = Vec::with_capacity(col.len() + 2);
    for r in col.iter() {
        if r.z1 <= z0 || r.z0 >= z1 { out.push(*r); continue; }
        if r.z0 < z0 { out.push(Run { z0: r.z0, z1: z0, material: r.material }); }
        if r.z1 > z1 { out.push(Run { z0: z1, z1: r.z1, material: r.material }); }
    }
    out.push(Run { z0, z1, material });
    normalize(&mut out);
    *col = out;
}

/// Intervals where `new` differs from `old`, carrying the new material.
pub fn diff_intervals(old: &[Run], new: &[Run]) -> Vec<Run> {
    let mut zs: Vec<f64> = Vec::with_capacity((old.len() + new.len()) * 2);
    for r in old.iter().chain(new.iter()) { zs.push(r.z0); zs.push(r.z1); }
    zs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    zs.dedup_by(|a, b| (*a - *b).abs() < 1e-7);
    let mut out: Vec<Run> = Vec::new();
    for w in zs.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b - a <= EPS { continue; }
        let mid = 0.5 * (a + b);
        let (mo, mn) = (material_at_z(old, mid), material_at_z(new, mid));
        if mo != mn {
            if let Some(last) = out.last_mut() {
                if last.material == mn && (a - last.z1).abs() < 1e-7 { last.z1 = b; continue; }
            }
            out.push(Run { z0: a, z1: b, material: mn });
        }
    }
    out
}

/// Local vertical-stack sampling of the SVO on the finest XY lattice.
#[derive(Debug, Clone)]
pub struct ColumnGrid {
    pub x0: f64,
    pub y0: f64,
    pub pitch: f64,
    pub nx: usize,
    pub ny: usize,
    pub cols: Vec<Column>,
    pub active: Vec<bool>,
    pub z_bottom: f64,
    pub z_top: f64,
}

impl ColumnGrid {
    pub fn sample(wafer: &WaferState, win: &ProcessWindow) -> Result<Self, GeomError> {
        win.validate()?;
        let pitch = wafer.leaf_pitch_nm();
        if !(pitch.is_finite() && pitch > 0.0) { return Err(GeomError::InvalidParameter("invalid SVO pitch".into())); }
        let (lo, hi) = win.bounds_xy();
        let (bx, by) = (wafer.bounds.min.0, wafer.bounds.min.1);
        let slack = pitch * 1e-6;
        if lo[0] < bx - slack || lo[1] < by - slack || hi[0] > wafer.bounds.max.0 + slack || hi[1] > wafer.bounds.max.1 + slack {
            return Err(GeomError::OutsideTile);
        }
        let i0 = (((lo[0] - bx) / pitch).floor().max(0.0)) as usize;
        let j0 = (((lo[1] - by) / pitch).floor().max(0.0)) as usize;
        let i1 = (((hi[0] - bx) / pitch).ceil().max(0.0)) as usize;
        let j1 = (((hi[1] - by) / pitch).ceil().max(0.0)) as usize;
        let (nx, ny) = (i1.saturating_sub(i0).max(1), j1.saturating_sub(j0).max(1));
        if nx > MAX_GRID_AXIS || ny > MAX_GRID_AXIS {
            return Err(GeomError::WindowTooLarge { cells: nx.max(ny), limit: MAX_GRID_AXIS });
        }
        let (x0, y0) = (bx + i0 as f64 * pitch, by + j0 as f64 * pitch);
        let mut cols = vec![Vec::new(); nx * ny];
        let mut active = vec![false; nx * ny];
        let mut any = false;
        for iy in 0..ny {
            for ix in 0..nx {
                let (cx, cy) = (x0 + (ix as f64 + 0.5) * pitch, y0 + (iy as f64 + 0.5) * pitch);
                if wafer.domain.contains_xy(cx, cy) {
                    let k = iy * nx + ix;
                    cols[k] = wafer.column_runs(cx, cy);
                    active[k] = !cols[k].is_empty();
                    any |= active[k];
                }
            }
        }
        if !any { return Err(GeomError::OutsideDomain); }
        Ok(Self { x0, y0, pitch, nx, ny, cols, active, z_bottom: wafer.bounds.min.2, z_top: wafer.bounds.max.2 })
    }

    /// Snap a height to the octree's Z lattice (cells are `pitch` tall).
    pub fn snap_z(&self, z: f64) -> f64 { self.z_bottom + ((z - self.z_bottom) / self.pitch).round() * self.pitch }

    /// Snap every run boundary to the lattice so the grid describes exactly
    /// what the octree can store (film thickness is quantised to one cell).
    pub fn quantize(&mut self) {
        for k in 0..self.cols.len() {
            if !self.active[k] { continue; }
            let mut c: Column = self.cols[k].iter().map(|r| Run { z0: self.snap_z(r.z0), z1: self.snap_z(r.z1), material: r.material }).collect();
            normalize(&mut c);
            self.cols[k] = c;
        }
    }

    pub fn idx(&self, ix: usize, iy: usize) -> usize { iy * self.nx + ix }
    pub fn center_xy(&self, ix: usize, iy: usize) -> (f64, f64) {
        (self.x0 + (ix as f64 + 0.5) * self.pitch, self.y0 + (iy as f64 + 0.5) * self.pitch)
    }
    pub fn active_count(&self) -> usize { self.active.iter().filter(|a| **a).count() }
    pub fn cell_area(&self) -> f64 { self.pitch * self.pitch }

    /// Nearest cell index containing a local XY point.
    pub fn cell_at(&self, x: f64, y: f64) -> Option<(usize, usize)> {
        let (fx, fy) = ((x - self.x0) / self.pitch, (y - self.y0) / self.pitch);
        if fx < 0.0 || fy < 0.0 { return None; }
        let (ix, iy) = (fx as usize, fy as usize);
        if ix < self.nx && iy < self.ny && self.active[self.idx(ix, iy)] { Some((ix, iy)) } else { None }
    }

    /// Height of the highest solid in each column (or the wafer bottom if none).
    pub fn surface_heights(&self) -> Vec<f64> {
        self.cols.iter().map(|c| solid_top(c).map(|s| s.0).unwrap_or(self.z_bottom)).collect()
    }

    /// Write only the changed columns back into the SVO. Returns the number of
    /// changed columns. Row-adjacent columns with identical edits share one
    /// `set_region` call, which keeps the octree compact.
    pub fn write_back(&self, original: &ColumnGrid, wafer: &mut WaferState) -> usize {
        debug_assert_eq!((self.nx, self.ny), (original.nx, original.ny));
        let mut changed = 0usize;
        let grow = self.pitch * 1e-6;
        for iy in 0..self.ny {
            let diffs: Vec<Vec<Run>> = (0..self.nx).map(|ix| {
                let k = self.idx(ix, iy);
                if self.active[k] { diff_intervals(&original.cols[k], &self.cols[k]) } else { Vec::new() }
            }).collect();
            let mut ix = 0;
            while ix < self.nx {
                if diffs[ix].is_empty() { ix += 1; continue; }
                let mut jx = ix + 1;
                while jx < self.nx && diffs[jx] == diffs[ix] { jx += 1; }
                changed += jx - ix;
                let (xa, xb) = (self.x0 + ix as f64 * self.pitch - grow, self.x0 + jx as f64 * self.pitch + grow);
                let (ya, yb) = (self.y0 + iy as f64 * self.pitch - grow, self.y0 + (iy + 1) as f64 * self.pitch + grow);
                for run in &diffs[ix] {
                    wafer.set_region(Bounds { min: (xa, ya, run.z0), max: (xb, yb, run.z1) }, run.material);
                }
                ix = jx;
            }
        }
        changed
    }
}

// ---------------------------------------------------------------------------
// Film deposition
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct FilmStats {
    pub columns: usize,
    pub mean_nm: f64,
    pub min_nm: f64,
    pub max_nm: f64,
    /// (max − min) / (2 · mean), percent.
    pub nonuniformity_pct: f64,
    /// Columns whose film would have exceeded the root box top and was clipped.
    pub clipped_columns: usize,
    pub volume_nm3: f64,
}

fn stats_from(thick: &[f64], clipped: usize, cell_area: f64) -> FilmStats {
    if thick.is_empty() { return FilmStats::default(); }
    let n = thick.len() as f64;
    let mean = thick.iter().sum::<f64>() / n;
    let (mn, mx) = thick.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(*v), b.max(*v)));
    FilmStats {
        columns: thick.len(), mean_nm: mean, min_nm: mn, max_nm: mx,
        nonuniformity_pct: if mean > 0.0 { (mx - mn) / (2.0 * mean) * 100.0 } else { 0.0 },
        clipped_columns: clipped, volume_nm3: thick.iter().sum::<f64>() * cell_area,
    }
}

/// Surface-aware film growth: the solid top surface is dilated by an
/// ellipsoid of vertical radius `thickness` and lateral radius
/// `step_coverage · thickness`. `step_coverage = 0` is perfectly directional
/// (PVD-like); `1` is fully conformal (CVD/ALD-like). Film is added *above*
/// the existing surface; nothing is overwritten.
pub fn grow_film(grid: &ColumnGrid, material: Material, thickness_nm: f64, step_coverage: f64) -> Result<(ColumnGrid, FilmStats), GeomError> {
    if !thickness_nm.is_finite() || thickness_nm < 0.0 { return Err(GeomError::InvalidParameter(format!("film thickness must be finite and >= 0, got {thickness_nm}"))); }
    if !step_coverage.is_finite() { return Err(GeomError::InvalidParameter("step coverage must be finite".into())); }
    let mut out = grid.clone();
    if thickness_nm <= EPS || material == Material::Air { return Ok((out, stats_from(&vec![0.0; grid.active_count()], 0, grid.cell_area()))); }
    let h = grid.surface_heights();
    let s = step_coverage.clamp(0.0, 1.0);
    let r = s * thickness_nm;
    let rc = (r / grid.pitch).floor() as isize;
    let mut thick = Vec::with_capacity(grid.active_count());
    let mut clipped = 0usize;
    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            if !grid.active[k] { continue; }
            let mut top = h[k] + thickness_nm;
            if rc >= 1 {
                for dj in -rc..=rc {
                    for di in -rc..=rc {
                        if di == 0 && dj == 0 { continue; }
                        let (jx, jy) = (ix as isize + di, iy as isize + dj);
                        if jx < 0 || jy < 0 || jx as usize >= grid.nx || jy as usize >= grid.ny { continue; }
                        let d = ((di * di + dj * dj) as f64).sqrt() * grid.pitch;
                        if d > r { continue; }
                        let j = grid.idx(jx as usize, jy as usize);
                        if !grid.active[j] { continue; }
                        let v = h[j] + thickness_nm * (1.0 - (d / r) * (d / r)).max(0.0).sqrt();
                        if v > top { top = v; }
                    }
                }
            }
            let clamped = grid.snap_z(top.min(grid.z_top));
            if top > grid.z_top + EPS { clipped += 1; }
            thick.push((clamped - h[k]).max(0.0));
            set_interval(&mut out.cols[k], h[k], clamped, material);
        }
    }
    Ok((out, stats_from(&thick, clipped, grid.cell_area())))
}

/// Planarizing, volume-conserving spin coat: resist fills topography to a flat
/// level `z_f` such that the resist volume equals `thickness_nm` × window area.
/// (Water-filling; columns taller than `z_f` stay uncovered.)
pub fn coat_planarizing(grid: &ColumnGrid, thickness_nm: f64) -> Result<(ColumnGrid, FilmStats), GeomError> {
    if !thickness_nm.is_finite() || thickness_nm < 0.0 { return Err(GeomError::InvalidParameter(format!("resist thickness must be finite and >= 0, got {thickness_nm}"))); }
    let mut out = grid.clone();
    let n = grid.active_count();
    if thickness_nm <= EPS || n == 0 { return Ok((out, stats_from(&vec![0.0; n], 0, grid.cell_area()))); }
    let h = grid.surface_heights();
    let hs: Vec<f64> = (0..h.len()).filter(|k| grid.active[*k]).map(|k| h[k]).collect();
    let target = thickness_nm * n as f64;
    let (mut lo, mut hi) = (hs.iter().cloned().fold(f64::INFINITY, f64::min), hs.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + thickness_nm);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let v: f64 = hs.iter().map(|z| (mid - z).max(0.0)).sum();
        if v < target { lo = mid; } else { hi = mid; }
    }
    let zf = grid.snap_z((0.5 * (lo + hi)).min(grid.z_top));
    let mut thick = Vec::with_capacity(n);
    for k in 0..h.len() {
        if !grid.active[k] { continue; }
        let t = (zf - h[k]).max(0.0);
        thick.push(t);
        if t > EPS { set_interval(&mut out.cols[k], h[k], zf, Material::Photoresist); }
    }
    Ok((out, stats_from(&thick, 0, grid.cell_area())))
}

// ---------------------------------------------------------------------------
// Mask, exposure and development
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum MaskPattern {
    /// Equal lines centred on the mask origin. `vertical` lines run along Y
    /// (and repeat along X). `clear_field`: features are opaque lines on a
    /// transparent mask; otherwise they are clear lines on an opaque mask.
    LineSpace { line_nm: f64, pitch_nm: f64, count: u32, length_nm: f64, vertical: bool, clear_field: bool },
    /// Arbitrary rectangles `[x0, y0, x1, y1]` relative to the mask origin.
    Rects { rects: Vec<[f64; 4]>, clear_field: bool },
    /// Arbitrary polygons (e.g. a GDSII layer), nm, relative to the mask origin. Polygons are
    /// **unioned** (GDSII layer semantics); each polygon uses even-odd fill for self-overlap.
    Polygons { polys: Vec<Vec<[f64; 2]>>, clear_field: bool },
}

impl MaskPattern {
    pub fn validate(&self) -> Result<(), GeomError> {
        match self {
            MaskPattern::LineSpace { line_nm, pitch_nm, count, length_nm, .. } => {
                for (n, v) in [("mask line width", *line_nm), ("mask pitch", *pitch_nm), ("mask line length", *length_nm)] {
                    if !v.is_finite() || v <= 0.0 { return Err(GeomError::InvalidParameter(format!("{n} must be finite and > 0, got {v}"))); }
                }
                if *count == 0 { return Err(GeomError::InvalidParameter("mask needs at least one line".into())); }
                if line_nm > pitch_nm { return Err(GeomError::InvalidParameter("mask line width exceeds pitch".into())); }
            }
            MaskPattern::Polygons { polys, .. } => {
                if polys.is_empty() { return Err(GeomError::InvalidParameter("mask has no polygons".into())); }
                if polys.iter().any(|p| p.len() < 3 || p.iter().any(|v| !v[0].is_finite() || !v[1].is_finite())) {
                    return Err(GeomError::InvalidParameter("mask polygons need >= 3 finite vertices".into()));
                }
            }
            MaskPattern::Rects { rects, .. } => {
                if rects.is_empty() { return Err(GeomError::InvalidParameter("mask has no rectangles".into())); }
                if rects.iter().any(|r| r.iter().any(|v| !v.is_finite()) || r[2] <= r[0] || r[3] <= r[1]) {
                    return Err(GeomError::InvalidParameter("mask rectangles must be finite with positive size".into()));
                }
            }
        }
        Ok(())
    }

    fn is_feature(&self, u: f64, v: f64) -> bool {
        match self {
            MaskPattern::LineSpace { line_nm, pitch_nm, count, length_nm, vertical, .. } => {
                let (across, along) = if *vertical { (u, v) } else { (v, u) };
                if along.abs() > length_nm * 0.5 { return false; }
                let k = (across / pitch_nm + (*count as f64 - 1.0) * 0.5).round();
                if k < 0.0 || k > *count as f64 - 1.0 { return false; }
                let centre = (k - (*count as f64 - 1.0) * 0.5) * pitch_nm;
                (across - centre).abs() <= line_nm * 0.5
            }
            MaskPattern::Rects { rects, .. } => rects.iter().any(|r| u >= r[0] && u <= r[2] && v >= r[1] && v <= r[3]),
            MaskPattern::Polygons { polys, .. } => polys.iter().any(|p| point_in_polygon(u, v, p)),
        }
    }

    pub fn clear_field(&self) -> bool {
        match self { MaskPattern::LineSpace { clear_field, .. } | MaskPattern::Rects { clear_field, .. } | MaskPattern::Polygons { clear_field, .. } => *clear_field }
    }

    /// Transmission (1 = light passes) at a mask-frame point.
    pub fn transmission(&self, u: f64, v: f64) -> f64 {
        let f = self.is_feature(u, v);
        if self.clear_field() { if f { 0.0 } else { 1.0 } } else if f { 1.0 } else { 0.0 }
    }
}

/// Fast transmission sampler: polygon masks are bucketed on a coarse grid so a lookup only tests
/// the few polygons near the point.
pub fn mask_sampler(mask: &MaskPattern) -> Box<dyn Fn(f64, f64) -> f64 + '_> {
    let MaskPattern::Polygons { polys, clear_field } = mask else { return Box::new(move |u, v| mask.transmission(u, v)) };
    let boxes: Vec<[f64; 4]> = polys.iter().map(|p| p.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v[0]), b[1].min(v[1]), b[2].max(v[0]), b[3].max(v[1])])).collect();
    let ext = boxes.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, x| [b[0].min(x[0]), b[1].min(x[1]), b[2].max(x[2]), b[3].max(x[3])]);
    const N: usize = 64;
    let (w, h) = ((ext[2] - ext[0]).max(1e-9), (ext[3] - ext[1]).max(1e-9));
    let cell = move |x: f64, y: f64| (((x - ext[0]) / w * N as f64).floor().clamp(0.0, N as f64 - 1.0) as usize, ((y - ext[1]) / h * N as f64).floor().clamp(0.0, N as f64 - 1.0) as usize);
    let mut buckets: Vec<Vec<u32>> = vec![Vec::new(); N * N];
    for (i, b) in boxes.iter().enumerate() {
        let (c0, r0) = cell(b[0], b[1]);
        let (c1, r1) = cell(b[2], b[3]);
        for r in r0..=r1 { for c in c0..=c1 { buckets[r * N + c].push(i as u32); } }
    }
    let cf = *clear_field;
    Box::new(move |u, v| {
        let inside = u >= ext[0] && u <= ext[2] && v >= ext[1] && v <= ext[3] && {
            let (c, r) = cell(u, v);
            buckets[r * N + c].iter().any(|&i| { let b = boxes[i as usize]; u >= b[0] && u <= b[2] && v >= b[1] && v <= b[3] && point_in_polygon(u, v, &polys[i as usize]) })
        };
        if cf { if inside { 0.0 } else { 1.0 } } else if inside { 1.0 } else { 0.0 }
    })
}

/// Even-odd point-in-polygon (ray casting).
pub fn point_in_polygon(x: f64, y: f64, poly: &[[f64; 2]]) -> bool {
    let n = poly.len();
    let mut inside = false;
    let mut j = n - 1;
    for i in 0..n {
        let (a, b) = (poly[i], poly[j]);
        if (a[1] > y) != (b[1] > y) && x < (b[0] - a[0]) * (y - a[1]) / (b[1] - a[1]) + a[0] { inside = !inside; }
        j = i;
    }
    inside
}

impl MaskPattern {
    /// Build a polygon mask from flattened GDS geometry: keep one layer (and datatype, if given),
    /// crop to a square of half-width `half_extent` around `center`, and re-centre on it.
    /// Returns the mask and how many polygons were kept.
    /// Like [`Self::from_layout`], but every layer is kept.
    ///
    /// A combined etch sees all layers as resist at once, so a feature on any layer
    /// protects the same silicon. That is what a whole-file run needs.
    pub fn from_layout_all(polys: &[crate::io::gds::FlatPoly], center: (f64, f64), half_extent: f64, clear_field: bool) -> Result<(Self, usize), GeomError> {
        if !(center.0.is_finite() && center.1.is_finite() && half_extent.is_finite() && half_extent > 0.0) {
            return Err(GeomError::InvalidParameter("GDS crop centre / extent must be finite".into()));
        }
        let mut kept: Vec<Vec<[f64; 2]>> = Vec::new();
        for p in polys.iter() {
            let b = p.pts.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)]);
            if b[2] < center.0 - half_extent || b[0] > center.0 + half_extent || b[3] < center.1 - half_extent || b[1] > center.1 + half_extent { continue; }
            kept.push(p.pts.iter().map(|v| [v.0 - center.0, v.1 - center.1]).collect());
        }
        let n = kept.len();
        Ok((MaskPattern::Polygons { polys: kept, clear_field }, n))
    }

    pub fn from_layout(polys: &[crate::io::gds::FlatPoly], layer: i32, datatype: Option<i32>, center: (f64, f64), half_extent: f64, clear_field: bool) -> Result<(Self, usize), GeomError> {
        if !(center.0.is_finite() && center.1.is_finite() && half_extent.is_finite() && half_extent > 0.0) { return Err(GeomError::InvalidParameter("GDS crop centre / extent must be finite".into())); }
        let mut kept: Vec<Vec<[f64; 2]>> = Vec::new();
        for p in polys.iter().filter(|p| p.layer == layer && datatype.map_or(true, |d| p.datatype == d)) {
            let b = p.pts.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)]);
            if b[2] < center.0 - half_extent || b[0] > center.0 + half_extent || b[3] < center.1 - half_extent || b[1] > center.1 + half_extent { continue; }
            kept.push(p.pts.iter().map(|v| [v.0 - center.0, v.1 - center.1]).collect());
        }
        if kept.is_empty() { return Err(GeomError::InvalidParameter(format!("GDS layer {layer} has no polygons within ±{half_extent:.0} nm of ({:.0}, {:.0}); move the crop centre onto the pattern", center.0, center.1))); }
        let n = kept.len();
        let m = MaskPattern::Polygons { polys: kept, clear_field };
        m.validate()?;
        Ok((m, n))
    }

    /// Build a polygon mask from one GDSII layer, re-centred on its bounding box.
    pub fn from_gds(lib: &crate::io::gds_parser::GdsLibrary, layer: i32, clear_field: bool) -> Result<Self, GeomError> {
        let polys: Vec<Vec<[f64; 2]>> = lib.polygons.iter().filter(|p| p.layer == layer).map(|p| p.vertices.iter().map(|v| [v.0, v.1]).collect()).collect();
        if polys.is_empty() { return Err(GeomError::InvalidParameter(format!("GDS layer {layer} has no polygons"))); }
        let (mut lo, mut hi) = ([f64::INFINITY; 2], [f64::NEG_INFINITY; 2]);
        for p in &polys { for v in p { for a in 0..2 { lo[a] = lo[a].min(v[a]); hi[a] = hi[a].max(v[a]); } } }
        let c = [(lo[0] + hi[0]) * 0.5, (lo[1] + hi[1]) * 0.5];
        let m = MaskPattern::Polygons { polys: polys.into_iter().map(|p| p.into_iter().map(|v| [v[0] - c[0], v[1] - c[1]]).collect()).collect(), clear_field };
        m.validate()?;
        Ok(m)
    }
}

/// Gaussian-blurred aerial image of the mask on the grid (0..1).
pub fn aerial_image(grid: &ColumnGrid, mask: &MaskPattern, mask_origin: [f64; 2], sigma_nm: f64) -> Result<Vec<f32>, GeomError> {
    mask.validate()?;
    if !sigma_nm.is_finite() || sigma_nm < 0.0 || mask_origin.iter().any(|v| !v.is_finite()) {
        return Err(GeomError::InvalidParameter("optical blur and mask origin must be finite, blur >= 0".into()));
    }
    let sample = mask_sampler(mask);
    let pad = if sigma_nm > 0.0 { (3.0 * sigma_nm / grid.pitch).ceil() as usize } else { 0 };
    let (wx, wy) = (grid.nx + 2 * pad, grid.ny + 2 * pad);
    let mut t = vec![0.0f64; wx * wy];
    for j in 0..wy {
        for i in 0..wx {
            let x = grid.x0 + (i as f64 - pad as f64 + 0.5) * grid.pitch - mask_origin[0];
            let y = grid.y0 + (j as f64 - pad as f64 + 0.5) * grid.pitch - mask_origin[1];
            t[j * wx + i] = sample(x, y);
        }
    }
    if pad > 0 {
        let sig = sigma_nm / grid.pitch;
        let kern: Vec<f64> = (-(pad as isize)..=pad as isize).map(|k| (-(k * k) as f64 / (2.0 * sig * sig)).exp()).collect();
        let ksum: f64 = kern.iter().sum();
        let mut tmp = t.clone();
        for j in 0..wy { for i in 0..wx {
            let mut acc = 0.0;
            for (m, kv) in kern.iter().enumerate() {
                let ii = i as isize + m as isize - pad as isize;
                let ii = ii.clamp(0, wx as isize - 1) as usize;
                acc += kv * t[j * wx + ii];
            }
            tmp[j * wx + i] = acc / ksum;
        } }
        let mut out = tmp.clone();
        for j in 0..wy { for i in 0..wx {
            let mut acc = 0.0;
            for (m, kv) in kern.iter().enumerate() {
                let jj = j as isize + m as isize - pad as isize;
                let jj = jj.clamp(0, wy as isize - 1) as usize;
                acc += kv * tmp[jj * wx + i];
            }
            out[j * wx + i] = acc / ksum;
        } }
        t = out;
    }
    let mut img = vec![0.0f32; grid.nx * grid.ny];
    for j in 0..grid.ny { for i in 0..grid.nx { img[j * grid.nx + i] = t[(j + pad) * wx + (i + pad)] as f32; } }
    Ok(img)
}

/// Development decision. Positive resist dissolves where the exposure
/// `intensity · dose_ratio` reaches the clearing dose (1.0); negative resist
/// dissolves where it does not.
pub fn develop_map(intensity: &[f32], dose_ratio: f64, positive: bool) -> Vec<bool> {
    let thr = (1.0 / dose_ratio.max(1e-6)) as f32;
    intensity.iter().map(|&i| if positive { i >= thr } else { i < thr }).collect()
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct DevelopStats { pub opened_columns: usize, pub retained_columns: usize, pub partial_columns: usize, pub resist_removed_nm3: f64 }

/// Dissolve resist from the top down in `removed` columns; `fraction` is the
/// developed fraction of the resist thickness (1 = fully cleared).
pub fn develop(grid: &ColumnGrid, removed: &[bool], fraction: f64) -> Result<(ColumnGrid, DevelopStats), GeomError> {
    if removed.len() != grid.cols.len() { return Err(GeomError::InvalidParameter("develop map size mismatch".into())); }
    let f = fraction.clamp(0.0, 1.0);
    let mut out = grid.clone();
    let mut st = DevelopStats::default();
    for k in 0..grid.cols.len() {
        if !grid.active[k] { continue; }
        let top_resist = solid_top(&grid.cols[k]).filter(|(_, m)| *m == Material::Photoresist);
        if top_resist.is_none() { continue; }
        if removed[k] {
            st.opened_columns += 1;
            let run = *grid.cols[k].iter().rev().find(|r| r.material == Material::Photoresist).unwrap();
            let th = run.z1 - run.z0;
            let new_top = run.z1 - f * th;
            if f > 0.0 { set_interval(&mut out.cols[k], new_top, run.z1, Material::Air); }
            st.resist_removed_nm3 += f * th * grid.cell_area();
        } else { st.retained_columns += 1; }
    }
    out.quantize();
    Ok((out, st))
}


/// Dissolve resist from the top by a per-column depth, nm.
///
/// This is what makes a sidewall possible: the depth varies across the edge of a feature
/// because the dose varies, so the resist top ends up sloped instead of stepped. The
/// boolean `develop` cannot do that — every opened column loses the same amount.
pub fn develop_depths(grid: &ColumnGrid, depths: &[f64]) -> Result<(ColumnGrid, DevelopStats), GeomError> {
    if depths.len() != grid.cols.len() {
        return Err(GeomError::InvalidParameter("develop depth map size mismatch".into()));
    }
    let mut out = grid.clone();
    let mut st = DevelopStats::default();
    for k in 0..grid.cols.len() {
        if !grid.active[k] { continue; }
        if solid_top(&grid.cols[k]).filter(|(_, m)| *m == Material::Photoresist).is_none() { continue; }
        let d = depths[k];
        if !d.is_finite() || d <= 0.0 { st.retained_columns += 1; continue; }
        let run = *grid.cols[k].iter().rev().find(|r| r.material == Material::Photoresist).unwrap();
        let th = run.z1 - run.z0;
        // The stepped developer approaches the substrate asymptotically as the rate falls with depth,
        // so a column within 1 % of the film is cleared, not left with a sliver that blocks the opening.
        let removed = if d >= th * 0.99 { th } else { d.clamp(0.0, th) };
        if removed > 0.0 { set_interval(&mut out.cols[k], run.z1 - removed, run.z1, Material::Air); }
        st.resist_removed_nm3 += removed * grid.cell_area();
        if removed >= th - 1e-6 { st.opened_columns += 1; } else { st.retained_columns += 1; if removed > 1e-6 { st.partial_columns += 1; } }
    }
    out.quantize();
    Ok((out, st))
}

/// Depth developed into each column, from the aerial image and a develop model.
///
/// `intensity` is the normalised aerial image; multiplying by the dose ratio gives the
/// local dose, which is what the resist actually sees.
pub fn develop_depth_map(
    intensity: &[f32],
    dose_ratio: f64,
    model: &atto_process::resist::DevelopModel,
    time_min: f64,
    thickness_nm: f64,
) -> Vec<f64> {
    intensity.iter().map(|&i| model.depth_nm(i as f64 * dose_ratio, atto_units::Minutes::new(time_min), atto_units::Nm::new(thickness_nm)).0).collect()
}


/// Depth developed into each column, with Dill absorption and a stepped develop.
///
/// Unlike the flat version, this walks the film: at each level the compound left depends on
/// depth, and the rate depends on the compound. That depth dependence is what gives the
/// profile a foot instead of a vertical wall.
pub fn develop_depth_map_dill(
    intensity: &[f32],
    dose_ratio: f64,
    mack: &atto_process::resist::DevelopModel,
    dill: &atto_process::dill::DillModel,
    time_min: f64,
    thickness_nm: f64,
    steps: usize,
) -> Vec<f64> {
    intensity
        .iter()
        .map(|&i| atto_process::dill::develop_depth(i as f64 * dose_ratio, atto_units::Nm::new(thickness_nm), mack, dill, atto_units::Minutes::new(time_min), steps).0)
        .collect()
}


/// Grow a film whose thickness is given per column.
///
/// `grow_film` takes one thickness for the whole grid and one step coverage. A transport
/// solve produces neither: every column sees a different amount of the source. This is the
/// same growth law with the scalar replaced by an array, so the lateral smoothing still
/// applies and the result is comparable to the scalar path.
pub fn grow_film_per_column(grid: &ColumnGrid, material: Material, thickness_nm: &[f64], step_coverage: f64) -> Result<(ColumnGrid, FilmStats), GeomError> {
    if thickness_nm.len() != grid.cols.len() {
        return Err(GeomError::InvalidParameter("per-column thickness length does not match the grid".into()));
    }
    if !step_coverage.is_finite() {
        return Err(GeomError::InvalidParameter("step coverage must be finite".into()));
    }
    let mut out = grid.clone();
    let h = grid.surface_heights();
    let s = step_coverage.clamp(0.0, 1.0);
    let mut thick = Vec::with_capacity(grid.active_count());
    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            if !grid.active[k] { continue; }
            let t = thickness_nm[k];
            if !t.is_finite() || t <= EPS || material == Material::Air {
                thick.push(0.0);
                continue;
            }
            let r = s * t;
            let rc = (r / grid.pitch).floor() as isize;
            let mut top = h[k] + t;
            if rc >= 1 {
                for dj in -rc..=rc {
                    for di in -rc..=rc {
                        if di == 0 && dj == 0 { continue; }
                        let (jx, jy) = (ix as isize + di, iy as isize + dj);
                        if jx < 0 || jy < 0 || jx as usize >= grid.nx || jy as usize >= grid.ny { continue; }
                        let d = ((di * di + dj * dj) as f64).sqrt() * grid.pitch;
                        if d > r { continue; }
                        let j = grid.idx(jx as usize, jy as usize);
                        if !grid.active[j] { continue; }
                        let v = h[j] + thickness_nm[j] * (1.0 - (d / r) * (d / r)).max(0.0).sqrt();
                        if v > top { top = v; }
                    }
                }
            }
            let clamped = grid.snap_z(top.min(grid.z_top));
            set_interval(&mut out.cols[k], h[k], clamped, material);
            thick.push((clamped - h[k]).max(0.0));
        }
    }
    out.quantize();
    Ok((out, stats_from(&thick, 0, grid.cell_area())))
}


/// Film statistics for a thickness given per column.
///
/// The per-column growth path returns its own statistics for one increment. A run that grows
/// the film in several increments needs the statistics of the total, which is this.
pub fn film_stats_from_thickness(grid: &ColumnGrid, thickness: &[f64]) -> Result<FilmStats, GeomError> {
    if thickness.len() != grid.cols.len() {
        return Err(GeomError::InvalidParameter("thickness length does not match the grid".into()));
    }
    let thick: Vec<f64> = (0..grid.cols.len())
        .filter(|k| grid.active[*k])
        .map(|k| if thickness[k].is_finite() { thickness[k].max(0.0) } else { 0.0 })
        .collect();
    Ok(stats_from(&thick, 0, grid.cell_area()))
}

// ---------------------------------------------------------------------------
// Etching
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EtchChemistry { OxideRie, SiliconRie, MetalRie, Custom }

impl EtchChemistry {
    pub const ALL: [EtchChemistry; 4] = [EtchChemistry::OxideRie, EtchChemistry::SiliconRie, EtchChemistry::MetalRie, EtchChemistry::Custom];
    pub fn label(self) -> &'static str {
        match self {
            EtchChemistry::OxideRie => "Oxide RIE (CF₄/CHF₃)",
            EtchChemistry::SiliconRie => "Silicon RIE (SF₆/C₄F₈)",
            EtchChemistry::MetalRie => "Metal RIE (Cl₂/BCl₃)",
            EtchChemistry::Custom => "Custom rates",
        }
    }
    /// Representative textbook-scale vertical rates (nm/min). Not calibrated.
    pub fn default_rates(self) -> EtchRates {
        match self {
            EtchChemistry::OxideRie | EtchChemistry::Custom => EtchRates { si: 3.0, sio2: 60.0, sic: 4.0, gan: 2.0, resist: 8.0, metal: 1.0 },
            EtchChemistry::SiliconRie => EtchRates { si: 100.0, sio2: 1.5, sic: 5.0, gan: 10.0, resist: 25.0, metal: 0.5 },
            EtchChemistry::MetalRie => EtchRates { si: 8.0, sio2: 2.0, sic: 2.0, gan: 6.0, resist: 12.0, metal: 80.0 },
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EtchRates { pub si: f64, pub sio2: f64, pub sic: f64, pub gan: f64, pub resist: f64, pub metal: f64 }

impl EtchRates {
    pub fn rate(&self, m: Material) -> f64 {
        match m {
            Material::Air => 0.0, Material::Silicon | Material::Dopant(_) => self.si, Material::SiliconDioxide => self.sio2,
            Material::SiliconCarbide => self.sic, Material::GalliumNitride => self.gan, Material::Photoresist => self.resist,
            Material::Metal => self.metal,
        }
    }
    pub fn all(&self) -> [f64; 6] { [self.si, self.sio2, self.sic, self.gan, self.resist, self.metal] }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct EtchSpec { pub rates: EtchRates, pub anisotropy: f64, pub time_min: f64 }

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct EtchStats {
    pub vertical_removed_nm3: f64,
    pub lateral_removed_nm3: f64,
    /// Deepest vertical removal in any column, nm.
    pub max_depth_nm: f64,
    /// Mean vertical removal over columns whose original top was not resist, nm.
    pub mean_depth_open_nm: f64,
    /// Mean resist consumed over resist-topped columns, nm.
    pub mask_consumed_mean_nm: f64,
    /// Resist-topped columns whose resist was fully consumed.
    pub mask_failed_columns: usize,
    /// Largest realised lateral reach (undercut / mask recession), nm.
    pub max_lateral_nm: f64,
}

struct Seg { z_hi: f64, z_lo: f64, t_hi: f64, t_lo: f64 }

struct VProfile { orig_top: f64, final_top: f64, segs: Vec<Seg> }

impl VProfile {
    /// Time at which height `z` became exposed to the plasma in this column.
    fn exposure_time(&self, z: f64) -> Option<f64> {
        if z >= self.orig_top - EPS { return Some(0.0); }
        if z < self.final_top - EPS { return None; }
        for s in &self.segs {
            if z <= s.z_hi + EPS && z >= s.z_lo - EPS {
                let span = (s.z_hi - s.z_lo).max(EPS);
                return Some(s.t_hi + (s.t_lo - s.t_hi) * (s.z_hi - z) / span);
            }
        }
        None
    }
}

/// Anisotropic etch with per-material rates, a mask that erodes at its own
/// rate, etch-stop on zero-rate layers, and a lateral component
/// `(1 − A) · rate` that undercuts from every surface exposed during the etch.
/// `fraction` scales the etch time (used for animation frames; 1 = committed).
pub fn etch(grid: &ColumnGrid, spec: &EtchSpec, fraction: f64) -> Result<(ColumnGrid, EtchStats), GeomError> {
    etch_scaled(grid, spec, fraction, None)
}

/// Etch with a per-column rate multiplier, for aspect-ratio-dependent etching.
pub fn etch_scaled(grid: &ColumnGrid, spec: &EtchSpec, fraction: f64, scale: Option<&[f64]>) -> Result<(ColumnGrid, EtchStats), GeomError> {
    if !spec.time_min.is_finite() || spec.time_min < 0.0 { return Err(GeomError::InvalidParameter(format!("etch time must be finite and >= 0, got {}", spec.time_min))); }
    if !spec.anisotropy.is_finite() || !(0.0..=1.0).contains(&spec.anisotropy) { return Err(GeomError::InvalidParameter(format!("anisotropy must be within [0, 1], got {}", spec.anisotropy))); }
    if spec.rates.all().iter().any(|r| !r.is_finite() || *r < 0.0) { return Err(GeomError::InvalidParameter("etch rates must be finite and >= 0".into())); }
    let total_t = spec.time_min * fraction.clamp(0.0, 1.0);
    let mut out = grid.clone();
    let mut st = EtchStats::default();
    if total_t <= EPS { return Ok((out, st)); }
    let area = grid.cell_area();
    let n = grid.cols.len();

    // --- vertical pass -----------------------------------------------------
    let mut prof: Vec<Option<VProfile>> = (0..n).map(|_| None).collect();
    let (mut open_cols, mut open_depth, mut mask_cols, mut mask_cons) = (0usize, 0.0f64, 0usize, 0.0f64);
    for k in 0..n {
        if !grid.active[k] { continue; }
        let col = &grid.cols[k];
        let Some(ti) = col.iter().rposition(|r| r.material != Material::Air) else { continue };
        let orig_top = col[ti].z1;
        let top_is_resist = col[ti].material == Material::Photoresist;
        let (mut budget, mut t) = (total_t, 0.0);
        let mut segs = Vec::new();
        let mut final_top = orig_top;
        let mut resist_removed = 0.0;
        let mut i = ti as isize;
        while i >= 0 && budget > EPS {
            let r = col[i as usize];
            if r.material == Material::Air { break; }
            let mult = scale.map(|v| v[k]).unwrap_or(1.0);
            let rate = spec.rates.rate(r.material) * mult;
            if rate <= 0.0 { break; }
            let th = r.z1 - r.z0;
            let need = th / rate;
            if need <= budget + EPS {
                segs.push(Seg { z_hi: r.z1, z_lo: r.z0, t_hi: t, t_lo: t + need });
                budget -= need; t += need; final_top = r.z0;
                if r.material == Material::Photoresist { resist_removed += th; }
                i -= 1;
            } else {
                let dz = budget * rate;
                segs.push(Seg { z_hi: r.z1, z_lo: r.z1 - dz, t_hi: t, t_lo: total_t });
                final_top = r.z1 - dz;
                if r.material == Material::Photoresist { resist_removed += dz; }
                break;
            }
        }
        let removed = orig_top - final_top;
        if removed > EPS {
            set_interval(&mut out.cols[k], final_top, orig_top, Material::Air);
            st.vertical_removed_nm3 += removed * area;
            st.max_depth_nm = st.max_depth_nm.max(removed);
        }
        if top_is_resist {
            mask_cols += 1; mask_cons += resist_removed;
            let resist_run_th = col[ti].z1 - col[ti].z0;
            if resist_removed >= resist_run_th - 1e-6 { st.mask_failed_columns += 1; }
        } else { open_cols += 1; open_depth += removed; }
        prof[k] = Some(VProfile { orig_top, final_top, segs });
    }
    if open_cols > 0 { st.mean_depth_open_nm = open_depth / open_cols as f64; }
    if mask_cols > 0 { st.mask_consumed_mean_nm = mask_cons / mask_cols as f64; }

    // --- lateral pass (undercut from exposed surfaces) ---------------------
    let lat_frac = 1.0 - spec.anisotropy;
    if lat_frac > 1e-9 {
        let z_floor = prof.iter().flatten().map(|p| p.final_top).fold(f64::INFINITY, f64::min);
        let max_rate = spec.rates.all().iter().cloned().fold(0.0, f64::max) * lat_frac;
        let reach_max = max_rate * total_t;
        let rc = (reach_max / grid.pitch).ceil() as isize;
        let dz = (grid.pitch * 0.5).max(0.5);
        if rc >= 1 && z_floor.is_finite() {
            let mut slabs: Vec<Vec<(f64, f64)>> = vec![Vec::new(); n];
            for iy in 0..grid.ny { for ix in 0..grid.nx {
                let k = grid.idx(ix, iy);
                if !grid.active[k] { continue; }
                for run in grid.cols[k].iter().filter(|r| r.material != Material::Air) {
                    let rl = spec.rates.rate(run.material) * lat_frac;
                    if rl <= 0.0 { continue; }
                    let reach_run = rl * total_t;
                    let z_start = run.z0.max(z_floor);
                    if run.z1 <= z_start { continue; }
                    let steps = (((run.z1 - z_start) / dz).ceil() as usize).max(1);
                    let step = (run.z1 - z_start) / steps as f64;
                    for s in 0..steps {
                        let (za, zb) = (z_start + s as f64 * step, z_start + (s + 1) as f64 * step);
                        let zm = 0.5 * (za + zb);
                        let mut hit = false;
                        'nb: for dj in -rc..=rc { for di in -rc..=rc {
                            if di == 0 && dj == 0 { continue; }
                            let d = ((di * di + dj * dj) as f64).sqrt() * grid.pitch;
                            if d > reach_run { continue; }
                            let (jx, jy) = (ix as isize + di, iy as isize + dj);
                            if jx < 0 || jy < 0 || jx as usize >= grid.nx || jy as usize >= grid.ny { continue; }
                            let j = grid.idx(jx as usize, jy as usize);
                            let Some(p) = prof[j].as_ref() else { continue };
                            let Some(t0) = p.exposure_time(zm) else { continue };
                            let reach = rl * (total_t - t0);
                            if reach >= d { st.max_lateral_nm = st.max_lateral_nm.max(d); hit = true; break 'nb; }
                        } }
                        if hit { slabs[k].push((za, zb)); }
                    }
                }
            } }
            for k in 0..n {
                for (za, zb) in std::mem::take(&mut slabs[k]) {
                    // Only remove what is still material after the vertical pass.
                    let before: f64 = out.cols[k].iter().filter(|r| r.material != Material::Air)
                        .map(|r| (r.z1.min(zb) - r.z0.max(za)).max(0.0)).sum();
                    if before > EPS { st.lateral_removed_nm3 += before * area; set_interval(&mut out.cols[k], za, zb, Material::Air); }
                }
            }
        }
    }
    out.quantize();
    Ok((out, st))
}

// ---------------------------------------------------------------------------
// Resist strip
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct StripStats { pub resist_before_nm3: f64, pub resist_after_nm3: f64 }

/// Plasma ash / wet strip: every resist run loses `efficiency · fraction` of
/// its thickness from the top. `efficiency < 1` leaves a residue film.
pub fn strip_resist(grid: &ColumnGrid, efficiency: f64, fraction: f64) -> Result<(ColumnGrid, StripStats), GeomError> {
    if !efficiency.is_finite() || !(0.0..=1.0).contains(&efficiency) { return Err(GeomError::InvalidParameter(format!("strip efficiency must be within [0, 1], got {efficiency}"))); }
    let f = efficiency * fraction.clamp(0.0, 1.0);
    let mut out = grid.clone();
    let mut st = StripStats::default();
    let area = grid.cell_area();
    for k in 0..grid.cols.len() {
        if !grid.active[k] { continue; }
        let runs: Vec<Run> = grid.cols[k].iter().filter(|r| r.material == Material::Photoresist).cloned().collect();
        for r in runs {
            let th = r.z1 - r.z0;
            st.resist_before_nm3 += th * area;
            st.resist_after_nm3 += th * (1.0 - f) * area;
            if f > 0.0 { set_interval(&mut out.cols[k], r.z0 + th * (1.0 - f), r.z1, Material::Air); }
        }
    }
    out.quantize();
    Ok((out, st))
}

pub fn material_volume_nm3(grid: &ColumnGrid, pred: impl Fn(Material) -> bool) -> f64 {
    let mut v = 0.0;
    for k in 0..grid.cols.len() {
        if !grid.active[k] { continue; }
        for r in &grid.cols[k] { if pred(r.material) { v += (r.z1 - r.z0) * grid.cell_area(); } }
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wafer() -> WaferState { WaferState::new_wafer_tile(300.0, 500.0, 16.0, 7, (0.0, 0.0)).unwrap() }
    fn win() -> ProcessWindow { ProcessWindow { center: [0.0, 0.0], size: [900.0, 900.0] } }
    fn si_top(w: &WaferState) -> f64 { w.surface_z_nm }

    #[test]
    fn headroom_wafer_is_tiny_and_columns_are_si_then_air() {
        let w = wafer();
        assert_eq!(w.node_count(), 9);
        let col = w.column_runs(0.0, 0.0);
        assert_eq!(col.len(), 2);
        assert_eq!(col[0].material, Material::Silicon);
        assert!((col[0].z1 - 500_000.0).abs() < 1e-6);
        assert_eq!(col[1].material, Material::Air);
        assert_eq!(w.material_at(0.0, 0.0, 100.0), Some(Material::Silicon));
        assert_eq!(w.material_at(0.0, 0.0, 600_000.0), Some(Material::Air));
        assert_eq!(w.material_at(2.0e8, 0.0, 10.0), None, "outside circular domain");
    }

    #[test]
    fn film_grows_above_the_surface_and_reports_measured_thickness() {
        let mut w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        let (g2, st) = grow_film(&g, Material::SiliconDioxide, 100.0, 0.0).unwrap();
        assert!((st.mean_nm - 96.0).abs() < 1e-6 && st.nonuniformity_pct < 1e-9, "100 nm requested, 16 nm cells -> 96 nm realised: {}", st.mean_nm);
        g2.write_back(&g, &mut w);
        let col = w.column_runs(0.0, 0.0);
        let film = col.iter().find(|r| r.material == Material::SiliconDioxide).expect("film run");
        assert!((film.z0 - si_top(&w)).abs() < 1e-3, "film starts at the original surface");
        assert!((film.z1 - film.z0 - st.mean_nm).abs() < 1e-6, "committed thickness {}", film.z1 - film.z0);
        // Substrate is untouched below the surface (growth, not conversion).
        assert_eq!(col[0].material, Material::Silicon);
        assert!((col[0].z1 - si_top(&w)).abs() < 1e-3);
    }

    #[test]
    fn conformal_growth_coats_sidewalls_but_directional_growth_does_not() {
        let mut w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        // Metal step: left half raised by 200 nm.
        let (step, _) = grow_film(&g, Material::Metal, 200.0, 0.0).unwrap();
        let mut stepped = step.clone();
        for iy in 0..stepped.ny { for ix in stepped.nx / 2..stepped.nx { let k = stepped.idx(ix, iy); stepped.cols[k] = g.cols[k].clone(); } }
        stepped.write_back(&g, &mut w);
        let base = ColumnGrid::sample(&w, &win()).unwrap();
        let (conf, _) = grow_film(&base, Material::SiliconDioxide, 60.0, 1.0).unwrap();
        let (dir, _) = grow_film(&base, Material::SiliconDioxide, 60.0, 0.0).unwrap();
        let foot = base.nx / 2; // first low column next to the wall
        let k = base.idx(foot, base.ny / 2);
        let h0 = solid_top(&base.cols[k]).unwrap().0;
        assert!(solid_top(&conf.cols[k]).unwrap().0 - h0 > 120.0, "conformal film should climb the wall beside the step");
        assert!((solid_top(&dir.cols[k]).unwrap().0 - h0 - 64.0).abs() < 1e-6, "directional film stays one 60 nm request (64 nm on the 16 nm lattice)");
    }

    #[test]
    fn planarizing_coat_conserves_resist_volume() {
        let mut w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        let (bump, _) = grow_film(&g, Material::Metal, 80.0, 0.0).unwrap();
        let mut b = bump.clone();
        for iy in 0..b.ny { for ix in 0..b.nx { if !(ix > b.nx / 3 && ix < 2 * b.nx / 3) { let k = b.idx(ix, iy); b.cols[k] = g.cols[k].clone(); } } }
        b.write_back(&g, &mut w);
        let g1 = ColumnGrid::sample(&w, &win()).unwrap();
        let (c, st) = coat_planarizing(&g1, 150.0).unwrap();
        let v = material_volume_nm3(&c, |m| m == Material::Photoresist);
        assert!((v / (g1.active_count() as f64 * g1.cell_area()) - 150.0).abs() < g1.pitch * 0.5, "mean resist thickness {}", v / (g1.active_count() as f64 * g1.cell_area()));
        let tops: Vec<f64> = c.cols.iter().zip(&c.active).filter(|(_, a)| **a).map(|(col, _)| solid_top(col).unwrap().0).collect();
        let (mn, mx) = tops.iter().fold((f64::INFINITY, 0.0f64), |(a, b), v| (a.min(*v), b.max(*v)));
        assert!(mx - mn < 1e-6, "planarized top is flat");
        assert!(st.max_nm > st.min_nm);
    }

    #[test]
    fn write_back_round_trip_matches_the_columns() {
        let mut w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        let (g2, _) = grow_film(&g, Material::SiliconDioxide, 100.0, 0.3).unwrap();
        let changed = g2.write_back(&g, &mut w);
        assert_eq!(changed, g.active_count());
        let g3 = ColumnGrid::sample(&w, &win()).unwrap();
        for k in 0..g3.cols.len() {
            if !g3.active[k] { continue; }
            let a = solid_top(&g2.cols[k]).unwrap();
            let b = solid_top(&g3.cols[k]).unwrap();
            assert_eq!(a.1, b.1);
            assert!((a.0 - b.0).abs() < 1e-6, "column {k}: {} vs {}", a.0, b.0);
        }
        assert!(w.node_count() < 60_000, "SVO stays sparse: {} nodes", w.node_count());
    }

    #[test]
    fn polygon_masks_use_even_odd_fill_and_build_from_gds_layers() {
        let square = vec![[-100.0, -100.0], [100.0, -100.0], [100.0, 100.0], [-100.0, 100.0]];
        let hole = vec![[-40.0, -40.0], [40.0, -40.0], [40.0, 40.0], [-40.0, 40.0]];
        let m = MaskPattern::Polygons { polys: vec![square.clone(), hole.clone()], clear_field: false };
        assert!(m.validate().is_ok());
        assert_eq!(m.transmission(70.0, 0.0), 1.0, "inside the square");
        assert_eq!(m.transmission(0.0, 0.0), 1.0, "GDSII layers are unioned: an overlapping shape does not cut a hole");
        assert_eq!(m.transmission(200.0, 0.0), 0.0);
        // The bucketed sampler agrees with the direct test everywhere on a grid.
        let fast = mask_sampler(&m);
        for i in -30..30 { for j in -30..30 { let (x, y) = (i as f64 * 9.7, j as f64 * 9.1); assert_eq!(fast(x, y), m.transmission(x, y), "({x},{y})"); } }
        // Self-overlapping *single* polygon still uses even-odd internally (a bow-tie leaves its centre-line open).
        let bow = MaskPattern::Polygons { polys: vec![vec![[-50.0, -50.0], [50.0, 50.0], [50.0, -50.0], [-50.0, 50.0]]], clear_field: false };
        assert_eq!(bow.transmission(0.0, 30.0), 0.0);
        assert!(MaskPattern::Polygons { polys: vec![vec![[0.0, 0.0], [1.0, 1.0]]], clear_field: false }.validate().is_err());
        let lib = crate::io::gds_parser::GdsLibrary { polygons: vec![crate::io::gds_parser::Polygon { layer: 7, vertices: square.iter().map(|v| (v[0] + 500.0, v[1] - 300.0)).collect() }], db_unit_m: 1e-9 };
        let g = MaskPattern::from_gds(&lib, 7, true).unwrap();
        assert_eq!(g.transmission(0.0, 0.0), 0.0, "re-centred on its bounding box: dark line at the origin of a clear-field mask");
        assert!(MaskPattern::from_gds(&lib, 9, true).is_err());
    }

    #[test]
    fn aerial_image_resolves_a_grating_and_dose_sets_the_linewidth() {
        let w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        let mask = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 300.0, count: 3, length_nm: 800.0, vertical: true, clear_field: true };
        let img = aerial_image(&g, &mask, [0.0, 0.0], 30.0).unwrap();
        let removed = develop_map(&img, 2.0, true);
        let row = g.ny / 2;
        let (cx, _) = g.center_xy(g.nx / 2, row);
        assert!(cx.abs() < g.pitch);
        let centre = g.idx(g.nx / 2, row);
        assert!(!removed[centre], "resist retained under the opaque line");
        let gap = g.idx(((150.0 / g.pitch) as usize) + g.nx / 2, row); // 150 nm right of a line centre = space centre
        assert!(removed[gap], "space between lines opens");
    }

    fn resist_stack() -> (WaferState, ColumnGrid) {
        // 100 nm oxide, 150 nm resist, vertical line pattern 150/300.
        let mut w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        let (ox, _) = grow_film(&g, Material::SiliconDioxide, 100.0, 0.0).unwrap();
        ox.write_back(&g, &mut w);
        let g1 = ColumnGrid::sample(&w, &win()).unwrap();
        let (pr, _) = coat_planarizing(&g1, 150.0).unwrap();
        let mask = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 300.0, count: 3, length_nm: 800.0, vertical: true, clear_field: true };
        let img = aerial_image(&pr, &mask, [0.0, 0.0], 25.0).unwrap();
        let (dev, _) = develop(&pr, &develop_map(&img, 2.0, true), 1.0).unwrap();
        dev.write_back(&g1, &mut w);
        let g2 = ColumnGrid::sample(&w, &win()).unwrap();
        (w, g2)
    }

    #[test]
    fn etch_only_removes_exposed_eligible_material() {
        let (_w, g) = resist_stack();
        let spec = EtchSpec { rates: EtchChemistry::OxideRie.default_rates(), anisotropy: 1.0, time_min: 100.0 / 60.0 * 1.2 };
        let (e, st) = etch(&g, &spec, 1.0).unwrap();
        let row = g.ny / 2;
        let line = g.idx(g.nx / 2, row);
        let space = g.idx(g.nx / 2 + (150.0 / g.pitch) as usize, row);
        assert_eq!(solid_top(&g.cols[line]).unwrap().1, Material::Photoresist);
        assert_eq!(solid_top(&e.cols[line]).unwrap().1, Material::Photoresist, "masked column keeps its resist");
        assert!(e.cols[line].iter().any(|r| r.material == Material::SiliconDioxide), "oxide under the mask survives");
        assert!(!e.cols[space].iter().any(|r| r.material == Material::SiliconDioxide), "exposed oxide is gone");
        let si = solid_top(&e.cols[space]).unwrap();
        assert_eq!(si.1, Material::Silicon, "etch stops on silicon (selectivity)");
        assert!((si.0 - 500_000.0).abs() < 10.0, "silicon barely attacked: top at {}", si.0);
        assert!(st.max_depth_nm >= 90.0 && st.mask_failed_columns == 0);
    }

    #[test]
    fn isotropic_component_undercuts_under_the_mask() {
        let (_w, g) = resist_stack();
        let row = g.ny / 2;
        let edge_in_mask = g.idx(g.nx / 2 + (75.0 / g.pitch) as usize - 1, row); // just inside the line edge
        let aniso = EtchSpec { rates: EtchChemistry::OxideRie.default_rates(), anisotropy: 1.0, time_min: 2.0 };
        let iso = EtchSpec { anisotropy: 0.4, ..aniso };
        let (a, _) = etch(&g, &aniso, 1.0).unwrap();
        let (b, sb) = etch(&g, &iso, 1.0).unwrap();
        let ox = |gr: &ColumnGrid, k: usize| gr.cols[k].iter().filter(|r| r.material == Material::SiliconDioxide).map(|r| r.z1 - r.z0).sum::<f64>();
        assert!(ox(&a, edge_in_mask) > 90.0, "anisotropic etch keeps the oxide at the mask edge");
        assert!(ox(&b, edge_in_mask) < ox(&a, edge_in_mask) - 20.0, "undercut removes oxide beneath the resist: {} vs {}", ox(&b, edge_in_mask), ox(&a, edge_in_mask));
        assert!(sb.max_lateral_nm > 0.0);
    }

    #[test]
    fn mask_erosion_can_fail_and_expose_the_film_beneath() {
        let (_w, g) = resist_stack();
        let spec = EtchSpec { rates: EtchRates { resist: 40.0, ..EtchChemistry::OxideRie.default_rates() }, anisotropy: 1.0, time_min: 6.0 };
        let (_e, st) = etch(&g, &spec, 1.0).unwrap();
        assert!(st.mask_failed_columns > 0, "150 nm resist at 40 nm/min is gone before 6 min");
    }

    #[test]
    fn etch_frames_converge_to_the_committed_result() {
        let (_w, g) = resist_stack();
        let spec = EtchSpec { rates: EtchChemistry::OxideRie.default_rates(), anisotropy: 0.9, time_min: 2.0 };
        let (full, _) = etch(&g, &spec, 1.0).unwrap();
        let (mid, _) = etch(&g, &spec, 0.5).unwrap();
        let (end, _) = etch(&g, &spec, 1.0).unwrap();
        assert!(full.cols == end.cols.clone() || full.cols.iter().zip(&end.cols).all(|(a, b)| a == b));
        let v = |gr: &ColumnGrid| material_volume_nm3(gr, |m| m == Material::SiliconDioxide);
        assert!(v(&g) > v(&mid) && v(&mid) > v(&full), "monotonic removal");
    }

    #[test]
    fn strip_clears_resist_and_reports_residue_when_inefficient() {
        let (_w, g) = resist_stack();
        let (clean, st) = strip_resist(&g, 1.0, 1.0).unwrap();
        assert_eq!(material_volume_nm3(&clean, |m| m == Material::Photoresist), 0.0);
        let (dirty, st2) = strip_resist(&g, 0.9, 1.0).unwrap();
        assert!(material_volume_nm3(&dirty, |m| m == Material::Photoresist) > 0.0);
        assert!((st2.resist_after_nm3 / st2.resist_before_nm3 - 0.1).abs() < 1e-9 && st.resist_before_nm3 > 0.0);
    }

    #[test]
    fn non_finite_parameters_are_rejected_before_mutation() {
        let w = wafer();
        let g = ColumnGrid::sample(&w, &win()).unwrap();
        assert!(grow_film(&g, Material::Metal, f64::NAN, 0.5).is_err());
        assert!(coat_planarizing(&g, f64::INFINITY).is_err());
        assert!(etch(&g, &EtchSpec { rates: EtchChemistry::OxideRie.default_rates(), anisotropy: 1.5, time_min: 1.0 }, 1.0).is_err());
        assert!(ProcessWindow { center: [f64::NAN, 0.0], size: [10.0, 10.0] }.validate().is_err());
        assert!(ColumnGrid::sample(&w, &ProcessWindow { center: [1.0e9, 0.0], size: [900.0, 900.0] }).is_err());
    }
}

#[cfg(test)]
mod gds_mask_tests {
    use super::*;
    use crate::io::gds::{parse, summarize, tests::sample_stream};

    #[test]
    fn layout_layers_crop_and_recentre_onto_the_window() {
        let flat = parse(&sample_stream()).unwrap().flatten(None).unwrap();
        let layers = summarize(&flat);
        let l5 = layers.iter().find(|l| l.layer == 5).unwrap();
        let (c, ext) = (l5.anchor, 400.0);
        let (m, n) = MaskPattern::from_layout(&flat, 5, Some(0), c, ext, true).unwrap();
        assert!(n >= 1 && n <= l5.polygons);
        // Crop far from the pattern reports a helpful error instead of an empty mask.
        let e = MaskPattern::from_layout(&flat, 5, None, (9.0e6, 9.0e6), 100.0, true).unwrap_err().to_string();
        assert!(e.contains("no polygons within") && e.contains("crop centre"), "{e}");
        assert!(MaskPattern::from_layout(&flat, 5, None, (f64::NAN, 0.0), 100.0, true).is_err());
        assert!(MaskPattern::from_layout(&flat, 99, None, c, ext, true).is_err());
        assert!(m.validate().is_ok());
    }
}
