// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Column operators: the transient sampling of the octree that every stage
//! transforms (build PRD A-06).

use atto_units::{Nm};
use crate::svo::{Bounds, Material, Run, Wafer};

const EPS: f64 = 1e-9;

#[derive(Debug, Clone, Copy, PartialEq)]
/// The region of the wafer a run covers.
pub struct ProcessWindow {
    pub center: [f64; 2],
    pub size: [f64; 2],
}

impl Default for ProcessWindow {
    fn default() -> Self { ProcessWindow { center: [0.0, 0.0], size: [1200.0, 1200.0] } }
}

impl ProcessWindow {
    /// The window corners in x and y.
    pub fn bounds_xy(&self) -> ([f64; 2], [f64; 2]) {
        ([self.center[0] - self.size[0] * 0.5, self.center[1] - self.size[1] * 0.5],
         [self.center[0] + self.size[0] * 0.5, self.center[1] + self.size[1] * 0.5])
    }
    /// The window area, nm^2.
    pub fn area_nm2(&self) -> f64 { self.size[0] * self.size[1] }
}

/// Merge adjacent runs of one material and drop empty ones.
pub fn normalize(col: &mut Vec<Run>) {
    col.retain(|r| r.z1 - r.z0 > EPS);
    col.sort_by(|a, b| a.z0.partial_cmp(&b.z0).unwrap_or(std::cmp::Ordering::Equal));
    let mut out: Vec<Run> = Vec::with_capacity(col.len());
    for r in col.drain(..) {
        if let Some(last) = out.last_mut() {
            if last.material == r.material && (r.z0 - last.z1).abs() < 1e-6 {
                last.z1 = r.z1;
                continue;
            }
        }
        out.push(r);
    }
    *col = out;
}

/// Overwrite `[z0, z1)` with `material`.
pub fn set_interval(col: &mut Vec<Run>, z0: f64, z1: f64, material: Material) {
    if z1 - z0 <= EPS { return; }
    let mut out: Vec<Run> = Vec::with_capacity(col.len() + 2);
    for r in col.iter() {
        if r.z1 <= z0 || r.z0 >= z1 {
            out.push(*r);
            continue;
        }
        if r.z0 < z0 { out.push(Run { z0: r.z0, z1: z0, material: r.material }); }
        if r.z1 > z1 { out.push(Run { z0: z1, z1: r.z1, material: r.material }); }
    }
    out.push(Run { z0, z1, material });
    normalize(&mut out);
    *col = out;
}

/// What is at a height in a column.
pub fn material_at_z(col: &[Run], z: f64) -> Material {
    col.iter().find(|r| z >= r.z0 && z < r.z1).map(|r| r.material).unwrap_or(Material::Air)
}

/// The topmost solid run and its height.
pub fn solid_top(col: &[Run]) -> Option<(f64, Material)> {
    col.iter().rev().find(|r| r.material.is_solid()).map(|r| (r.z1, r.material))
}

/// Intervals where `new` differs from `old`, carrying the new material.
pub fn diff_intervals(old: &[Run], new: &[Run]) -> Vec<Run> {
    let mut zs: Vec<f64> = Vec::with_capacity((old.len() + new.len()) * 2);
    for r in old.iter().chain(new.iter()) {
        zs.push(r.z0);
        zs.push(r.z1);
    }
    zs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    zs.dedup_by(|a, b| (*a - *b).abs() < 1e-7);
    let mut out: Vec<Run> = Vec::new();
    for w in zs.windows(2) {
        let (a, b) = (w[0], w[1]);
        if b - a <= EPS { continue; }
        let mid = 0.5 * (a + b);
        let (mo, mn) = (material_at_z(old, mid), material_at_z(new, mid));
        if mo != mn {
            if let Some(last) = out.last_mut() {
                if last.material == mn && (a - last.z1).abs() < 1e-7 {
                    last.z1 = b;
                    continue;
                }
            }
            out.push(Run { z0: a, z1: b, material: mn });
        }
    }
    out
}

#[derive(Debug, Clone)]
/// A column representation of a wafer region: one stack per cell.
pub struct ColumnGrid {
    pub x0: f64,
    pub y0: f64,
    pub pitch: f64,
    pub nx: usize,
    pub ny: usize,
    pub cols: Vec<Vec<Run>>,
    pub active: Vec<bool>,
    pub z_bottom: f64,
    pub z_top: f64,
}

impl ColumnGrid {
    /// Check the grid's invariants, returning every violation (PRD A-11).
    ///
    /// These are the states a half-applied stage leaves: a run out of order, a
    /// run with no thickness, a column that is not the same length as the index
    /// space it is addressed with.
    pub fn check(&self) -> Vec<String> {
        let mut out = Vec::new();
        if !(self.pitch > 0.0) || !self.pitch.is_finite() {
            out.push(format!("pitch is not a positive number: {}", self.pitch));
        }
        if self.cols.len() != self.nx * self.ny {
            out.push(format!("{} columns for a {}x{} grid", self.cols.len(), self.nx, self.ny));
        }
        if self.active.len() != self.nx * self.ny {
            out.push(format!("{} active flags for a {}x{} grid", self.active.len(), self.nx, self.ny));
        }
        for (k, col) in self.cols.iter().enumerate() {
            let mut prev_top: Option<f64> = None;
            for run in col {
                if !run.z0.is_finite() || !run.z1.is_finite() {
                    out.push(format!("column {k} has a run with a non-finite height"));
                    break;
                }
                if run.z1 < run.z0 {
                    out.push(format!("column {k} has a run that ends below where it starts"));
                    break;
                }
                if let Some(t) = prev_top {
                    if run.z0 < t - 1e-9 {
                        out.push(format!("column {k} has overlapping runs"));
                        break;
                    }
                }
                prev_top = Some(run.z1);
            }
        }
        out
    }

    /// Sample a wafer into a grid over a window.
    pub fn sample(w: &Wafer, win: &ProcessWindow) -> Result<Self, String> {
        let pitch = w.leaf_pitch_nm();
        if !(pitch.is_finite() && pitch > 0.0) { return Err("invalid leaf pitch".into()); }
        let (lo, hi) = win.bounds_xy();
        let b = w.bounds;
        if lo[0] < b.min.0 || lo[1] < b.min.1 || hi[0] > b.max.0 || hi[1] > b.max.1 {
            return Err(format!("window {:?}..{:?} leaves the tile {:?}..{:?}", lo, hi, b.min, b.max));
        }
        let i0 = ((lo[0] - b.min.0) / pitch).floor().max(0.0) as usize;
        let j0 = ((lo[1] - b.min.1) / pitch).floor().max(0.0) as usize;
        let i1 = ((hi[0] - b.min.0) / pitch).ceil().max(0.0) as usize;
        let j1 = ((hi[1] - b.min.1) / pitch).ceil().max(0.0) as usize;
        let (nx, ny) = (i1.saturating_sub(i0).max(1), j1.saturating_sub(j0).max(1));
        if nx > 512 || ny > 512 { return Err(format!("window needs {nx}x{ny} cells, limit 512")); }
        let (x0, y0) = (b.min.0 + i0 as f64 * pitch, b.min.1 + j0 as f64 * pitch);
        let mut cols = vec![Vec::new(); nx * ny];
        let mut active = vec![false; nx * ny];
        let mut any = false;
        for iy in 0..ny {
            for ix in 0..nx {
                let (cx, cy) = (x0 + (ix as f64 + 0.5) * pitch, y0 + (iy as f64 + 0.5) * pitch);
                if w.contains_xy(cx, cy) {
                    let k = iy * nx + ix;
                    cols[k] = w.column_runs(cx, cy);
                    active[k] = !cols[k].is_empty();
                    any |= active[k];
                }
            }
        }
        if !any { return Err("window lies outside the wafer disk".into()); }
        Ok(ColumnGrid { x0, y0, pitch, nx, ny, cols, active, z_bottom: b.min.2, z_top: b.max.2 })
    }

    /// The column index of a cell.
    pub fn idx(&self, ix: usize, iy: usize) -> usize { iy * self.nx + ix }
    /// The area of one cell, nm^2.
    pub fn cell_area(&self) -> f64 { self.pitch * self.pitch }
    pub fn active_count(&self) -> usize { self.active.iter().filter(|a| **a).count() }
    /// The nearest cell boundary to a height.
    pub fn snap_z(&self, z: f64) -> f64 { self.z_bottom + ((z - self.z_bottom) / self.pitch).round() * self.pitch }

    /// The top-of-solid height of every active cell.
    pub fn surface_heights(&self) -> Vec<f64> {
        self.cols.iter().map(|c| solid_top(c).map(|s| s.0).unwrap_or(self.z_bottom)).collect()
    }

    /// Write only the changed columns back into the octree.
    pub fn write_back(&self, original: &ColumnGrid, w: &mut Wafer) -> usize {
        let grow = self.pitch * 1e-6;
        let mut changed = 0;
        for iy in 0..self.ny {
            let diffs: Vec<Vec<Run>> = (0..self.nx)
                .map(|ix| {
                    let k = self.idx(ix, iy);
                    if self.active[k] { diff_intervals(&original.cols[k], &self.cols[k]) } else { Vec::new() }
                })
                .collect();
            let mut ix = 0;
            while ix < self.nx {
                if diffs[ix].is_empty() { ix += 1; continue; }
                let mut jx = ix + 1;
                while jx < self.nx && diffs[jx] == diffs[ix] { jx += 1; }
                changed += jx - ix;
                let xa = self.x0 + ix as f64 * self.pitch - grow;
                let xb = self.x0 + jx as f64 * self.pitch + grow;
                let ya = self.y0 + iy as f64 * self.pitch - grow;
                let yb = self.y0 + (iy + 1) as f64 * self.pitch + grow;
                for run in &diffs[ix] {
                    w.set_region(Bounds::new((xa, ya, run.z0), (xb, yb, run.z1)), run.material);
                }
                ix = jx;
            }
        }
        changed
    }

    /// A 1-D max-plus pass with a linear slope: `out[i] = max_j (v[j] - step*|i-j|)`.
    fn cone_1d(v: &[f64], step: f64, out: &mut [f64]) {
        let n = v.len();
        out.copy_from_slice(v);
        for i in 1..n {
            let cand = out[i - 1] - step;
            if cand > out[i] { out[i] = cand; }
        }
        for i in (0..n.saturating_sub(1)).rev() {
            let cand = out[i + 1] - step;
            if cand > out[i] { out[i] = cand; }
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// A deposited film's thickness and coverage.
pub struct FilmStats {
    pub columns: usize,
    pub mean_nm: f64,
    pub min_nm: f64,
    pub max_nm: f64,
    pub nonuniformity_pct: f64,
    pub volume_nm3: f64,
}

fn film_stats(thick: &[f64], area: f64) -> FilmStats {
    if thick.is_empty() { return FilmStats::default(); }
    let n = thick.len() as f64;
    let mean = thick.iter().sum::<f64>() / n;
    let (mn, mx) = thick.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(*v), b.max(*v)));
    FilmStats {
        columns: thick.len(),
        mean_nm: mean,
        min_nm: mn,
        max_nm: mx,
        nonuniformity_pct: if mean > 0.0 { (mx - mn) / (2.0 * mean) * 100.0 } else { 0.0 },
        volume_nm3: thick.iter().sum::<f64>() * area,
    }
}

/// Surface-aware film growth.
///
/// `coverage = 0` is directional (thickness everywhere above the local surface);
/// `coverage = 1` is fully conformal (the film follows the surface within a
/// lateral reach). The lateral reach is a *diamond* rather than a circle: that
/// makes the dilation separable, so the cost is O(reach) per column instead of
/// O(reach^2) (audit F-31). The shape difference is bounded and tested.
pub fn grow_film(g: &ColumnGrid, material: Material, thickness: Nm, coverage: f64) -> Result<(ColumnGrid, FilmStats), String> {
    let thickness_nm = thickness.value();
    if !thickness_nm.is_finite() || thickness_nm < 0.0 { return Err("thickness must be finite and >= 0".into()); }
    if !coverage.is_finite() { return Err("coverage must be finite".into()); }
    let mut out = g.clone();
    let s = coverage.clamp(0.0, 1.0);
    if thickness_nm <= EPS {
        return Ok((out, film_stats(&vec![0.0; g.active_count()], g.cell_area())));
    }
    let h = g.surface_heights();
    let reach = s * thickness_nm;
    let step = if reach > g.pitch { thickness_nm * g.pitch / reach } else { f64::INFINITY };

    // cone[i] = max_j (h[j] + t - step*L1(i,j)) over the grid, separable.
    let mut cone = h.clone();
    if step.is_finite() {
        let mut line = vec![0.0; g.nx.max(g.ny)];
        let mut buf = vec![0.0; g.nx.max(g.ny)];
        for iy in 0..g.ny {
            let mut v = vec![0.0; g.nx];
            for ix in 0..g.nx {
                v[ix] = if g.active[g.idx(ix, iy)] { h[g.idx(ix, iy)] + thickness_nm } else { f64::NEG_INFINITY };
            }
            let n = g.nx;
            line[..n].copy_from_slice(&v);
            ColumnGrid::cone_1d(&line[..n], step, &mut buf[..n]);
            for ix in 0..g.nx {
                cone[g.idx(ix, iy)] = buf[ix];
            }
        }
        let mut col = vec![0.0; g.ny];
        let mut cbuf = vec![0.0; g.ny];
        for ix in 0..g.nx {
            for iy in 0..g.ny { col[iy] = cone[g.idx(ix, iy)]; }
            ColumnGrid::cone_1d(&col, step, &mut cbuf);
            for iy in 0..g.ny { cone[g.idx(ix, iy)] = cbuf[iy]; }
        }
    }

    let mut thick = Vec::with_capacity(g.active_count());
    for k in 0..g.cols.len() {
        if !g.active[k] { continue; }
        let directional = h[k] + thickness_nm;
        let top_raw = if step.is_finite() {
            (1.0 - s) * directional + s * cone[k].max(directional - thickness_nm)
        } else {
            directional
        };
        let top = g.snap_z(top_raw.min(g.z_top));
        let t = (top - h[k]).max(0.0);
        thick.push(t);
        set_interval(&mut out.cols[k], h[k], top, material);
    }
    Ok((out, film_stats(&thick, g.cell_area())))
}

/// Planarising spin coat: resist fills to a flat level such that the volume
/// equals `thickness_nm` times the window area (water filling).
pub fn coat_planarizing(g: &ColumnGrid, thickness_nm: f64) -> Result<(ColumnGrid, FilmStats), String> {
    if !thickness_nm.is_finite() || thickness_nm < 0.0 { return Err("thickness must be finite and >= 0".into()); }
    let mut out = g.clone();
    let n = g.active_count();
    if thickness_nm <= EPS || n == 0 {
        return Ok((out, film_stats(&vec![0.0; n], g.cell_area())));
    }
    let h = g.surface_heights();
    let hs: Vec<f64> = (0..h.len()).filter(|k| g.active[*k]).map(|k| h[k]).collect();
    let target = thickness_nm * n as f64;
    let mut lo = hs.iter().cloned().fold(f64::INFINITY, f64::min);
    let mut hi = hs.iter().cloned().fold(f64::NEG_INFINITY, f64::max) + thickness_nm;
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        let v: f64 = hs.iter().map(|z| (mid - z).max(0.0)).sum();
        if v < target { lo = mid } else { hi = mid }
    }
    let zf = g.snap_z((0.5 * (lo + hi)).min(g.z_top));
    let mut thick = Vec::with_capacity(n);
    for k in 0..h.len() {
        if !g.active[k] { continue; }
        let t = (zf - h[k]).max(0.0);
        thick.push(t);
        if t > EPS { set_interval(&mut out.cols[k], h[k], zf, Material::Photoresist); }
    }
    Ok((out, film_stats(&thick, g.cell_area())))
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// What a development step removed.
pub struct DevelopStats {
    pub opened_columns: usize,
    pub retained_columns: usize,
    pub resist_removed_nm3: f64,
}

/// Remove the top `fraction` of the resist in every opened column.
pub fn develop(g: &ColumnGrid, removed: &[bool], fraction: f64) -> Result<(ColumnGrid, DevelopStats), String> {
    if removed.len() != g.cols.len() { return Err("develop map size mismatch".into()); }
    let f = fraction.clamp(0.0, 1.0);
    let mut out = g.clone();
    let mut st = DevelopStats::default();
    for k in 0..g.cols.len() {
        if !g.active[k] { continue; }
        let top_resist = solid_top(&g.cols[k]).filter(|(_, m)| *m == Material::Photoresist);
        if top_resist.is_none() { continue; }
        if removed[k] {
            st.opened_columns += 1;
            let run = match g.cols[k].iter().rev().find(|r| r.material == Material::Photoresist) {
                Some(r) => *r,
                None => continue,
            };
            let th = run.thickness();
            let new_top = run.z1 - f * th;
            if f > 0.0 { set_interval(&mut out.cols[k], new_top, run.z1, Material::Air); }
            st.resist_removed_nm3 += f * th * g.cell_area();
        } else {
            st.retained_columns += 1;
        }
    }
    Ok((out, st))
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Vertical and lateral rates, nm/min.
pub struct EtchRates {
    pub si: f64,
    pub sio2: f64,
    pub si3n4: f64,
    pub resist: f64,
    pub metal: f64,
}

impl EtchRates {
    /// The vertical rate for a material, nm/min.
    pub fn rate(&self, m: Material) -> f64 {
        match m {
            Material::Air => 0.0,
            Material::Silicon => self.si,
            Material::SiliconDioxide => self.sio2,
            Material::SiliconNitride => self.si3n4,
            Material::Photoresist => self.resist,
            Material::Metal => self.metal,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// An etch: its rates, its anisotropy and how long it runs.
pub struct EtchSpec {
    pub rates: EtchRates,
    pub time_min: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// What an etch removed, and where it stopped.
pub struct EtchStats {
    /// Material removed from columns whose top is not resist (the feature depth).
    pub feature_depth_mean_nm: f64,
    pub feature_depth_max_nm: f64,
    /// Resist consumed on masked columns, kept separate (audit F-23).
    pub mask_consumed_mean_nm: f64,
    pub mask_failed_columns: usize,
    pub removed_volume_nm3: f64,
    pub open_columns: usize,
}

/// Directional etch: descend the column at the material's rate until the time
/// budget is spent. Etch stops on a zero-rate layer, continues past voids
/// (audit F-24), and reports mask consumption separately from feature depth.
pub fn etch(g: &ColumnGrid, spec: &EtchSpec, fraction: f64) -> Result<(ColumnGrid, EtchStats), String> {
    etch_scaled(g, spec, fraction, None)
}

/// Etch with a per-column rate multiplier.
///
/// `scale[k]` multiplies the rate for column `k`. This is where aspect-ratio-dependent
/// etching enters: a narrow deep hole etches slower than a wide shallow one because the
/// flux reaching its bottom is reduced. The physics lives in the caller; this only applies
/// the factor, so the geometric core stays testable on its own.
pub fn etch_scaled(g: &ColumnGrid, spec: &EtchSpec, fraction: f64, scale: Option<&[f64]>) -> Result<(ColumnGrid, EtchStats), String> {
    if !spec.time_min.is_finite() || spec.time_min < 0.0 { return Err("time must be finite and >= 0".into()); }
    for r in [spec.rates.si, spec.rates.sio2, spec.rates.si3n4, spec.rates.resist, spec.rates.metal] {
        if !r.is_finite() || r < 0.0 { return Err("rates must be finite and >= 0".into()); }
    }
    let total = spec.time_min * fraction.clamp(0.0, 1.0);
    let mut out = g.clone();
    let mut st = EtchStats::default();
    if total <= EPS { return Ok((out, st)); }
    let area = g.cell_area();
    let (mut open_n, mut open_sum, mut mask_n, mut mask_sum) = (0usize, 0.0f64, 0usize, 0.0f64);
    for k in 0..g.cols.len() {
        if !g.active[k] { continue; }
        let col = &g.cols[k];
        let ti = match col.iter().rposition(|r| r.material.is_solid()) {
            Some(i) => i,
            None => continue,
        };
        let orig_top = col[ti].z1;
        let top_is_resist = col[ti].material == Material::Photoresist;
        let (mut budget, mut final_top, mut resist_removed) = (total, orig_top, 0.0f64);
        let mut i = ti as isize;
        while i >= 0 && budget > EPS {
            let r = col[i as usize];
            if !r.material.is_solid() { i -= 1; continue; }
            let mult = scale.map(|v| v[k]).unwrap_or(1.0);
            let rate = spec.rates.rate(r.material) * mult;
            if rate <= 0.0 { break; }
            let th = r.z1 - r.z0;
            let need = th / rate;
            if need <= budget + EPS {
                budget -= need;
                final_top = r.z0;
                if r.material == Material::Photoresist { resist_removed += th; }
                i -= 1;
            } else {
                let dz = budget * rate;
                final_top = r.z1 - dz;
                if r.material == Material::Photoresist { resist_removed += dz; }
                break;
            }
        }
        let removed = orig_top - final_top;
        if removed > EPS {
            set_interval(&mut out.cols[k], final_top, orig_top, Material::Air);
            st.removed_volume_nm3 += removed * area;
        }
        if top_is_resist {
            mask_n += 1;
            mask_sum += resist_removed;
            if resist_removed >= (col[ti].z1 - col[ti].z0) - 1e-6 { st.mask_failed_columns += 1; }
        } else {
            open_n += 1;
            open_sum += removed;
            st.feature_depth_max_nm = st.feature_depth_max_nm.max(removed);
        }
    }
    if open_n > 0 { st.feature_depth_mean_nm = open_sum / open_n as f64; }
    if mask_n > 0 { st.mask_consumed_mean_nm = mask_sum / mask_n as f64; }
    st.open_columns = open_n;
    Ok((out, st))
}

/// Remove the top `efficiency` fraction of every resist run.
pub fn strip_resist(g: &ColumnGrid, efficiency: f64) -> Result<(ColumnGrid, f64), String> {
    if !efficiency.is_finite() || !(0.0..=1.0).contains(&efficiency) { return Err("efficiency must be within [0,1]".into()); }
    let mut out = g.clone();
    let mut removed = 0.0;
    for k in 0..g.cols.len() {
        if !g.active[k] { continue; }
        let runs: Vec<Run> = g.cols[k].iter().copied().filter(|r| r.material == Material::Photoresist).collect();
        for r in runs {
            let th = r.z1 - r.z0;
            let f = efficiency;
            removed += th * f * g.cell_area();
            if f > 0.0 { set_interval(&mut out.cols[k], r.z1 - th * f, r.z1, Material::Air); }
        }
    }
    Ok((out, removed))
}

/// A buried void: an Air run with material above it (build PRD S-05).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Void {
    /// Column position, nm.
    pub x_nm: f64,
    pub y_nm: f64,
    /// Top of the void, nm.
    pub z_nm: f64,
    /// Void height, nm.
    pub height_nm: f64,
}

/// Find buried voids: air with solid above it, which is a cavity rather than an
/// open surface.
pub fn find_voids(grid: &ColumnGrid) -> Vec<Void> {
    let mut out = Vec::new();
    for k in 0..grid.cols.len() {
        if !grid.active[k] {
            continue;
        }
        let col = &grid.cols[k];
        let top = match solid_top(col) {
            Some((t, _)) => t,
            None => continue,
        };
        for r in col.iter() {
            if r.material != Material::Air || r.z1 > top {
                continue;
            }
            // air below the topmost solid is enclosed
            if r.z1 - r.z0 > EPS {
                out.push(Void {
                    x_nm: grid.x0 + (k % grid.nx) as f64 * grid.pitch,
                    y_nm: grid.y0 + (k / grid.nx) as f64 * grid.pitch,
                    z_nm: r.z1,
                    height_nm: r.z1 - r.z0,
                });
            }
        }
    }
    out
}

/// Width of the narrowest recess between features along a row, nm (PRD S-05).
///
/// A seam is not the absence of material: the gap column usually still has
/// substrate at its top, just far below its neighbours. It is a run of columns
/// sitting at least one cell below the row's highest feature, bounded on both
/// sides by that feature.
pub fn seam_width_nm(grid: &ColumnGrid) -> Option<f64> {
    if grid.nx < 3 || grid.ny == 0 {
        return None;
    }
    let row = grid.ny / 2;
    let tops: Vec<Option<f64>> = (0..grid.nx)
        .map(|ix| {
            let k = grid.idx(ix, row);
            if grid.active[k] { solid_top(&grid.cols[k]).map(|(t, _)| t) } else { None }
        })
        .collect();
    let highest = tops.iter().flatten().cloned().fold(f64::MIN, f64::max);
    if !highest.is_finite() {
        return None;
    }
    let low = |ix: usize| match tops[ix] {
        Some(t) => t < highest - grid.pitch * 0.5,
        None => true,
    };
    let high = |ix: usize| match tops[ix] {
        Some(t) => t >= highest - grid.pitch * 0.5,
        None => false,
    };
    let mut best: Option<f64> = None;
    let mut ix = 1;
    while ix < grid.nx - 1 {
        if low(ix) && high(ix - 1) {
            let start = ix;
            while ix < grid.nx && low(ix) {
                ix += 1;
            }
            if ix < grid.nx && high(ix) {
                let w = (ix - start) as f64 * grid.pitch;
                best = Some(best.map_or(w, |v: f64| v.min(w)));
            }
        } else {
            ix += 1;
        }
    }
    best
}

#[cfg(test)]
mod void_tests {
    use atto_units::{Mm, Um};
    use super::*;
    use crate::svo::{Material, Wafer};

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        ColumnGrid::sample(&w, &win).unwrap()
    }

    #[test]
    fn a_solid_column_has_no_void() {
        let mut g = grid();
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top, top + 100.0, Material::SiliconDioxide);
        }
        assert!(find_voids(&g).is_empty(), "nothing buried, nothing reported");
    }

    #[test]
    fn a_buried_air_run_is_a_void() {
        let mut g = grid();
        let k = g.idx(5, 5);
        let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
        // oxide, a gap, then more oxide on top: an enclosed cavity
        set_interval(&mut g.cols[k], top, top + 40.0, Material::SiliconDioxide);
        set_interval(&mut g.cols[k], top + 40.0, top + 60.0, Material::Air);
        set_interval(&mut g.cols[k], top + 60.0, top + 100.0, Material::SiliconDioxide);
        let v = find_voids(&g);
        assert_eq!(v.len(), 1, "one cavity, got {v:?}");
        assert!((v[0].height_nm - 20.0).abs() < 1e-6, "height {:.1} nm", v[0].height_nm);
    }

    #[test]
    fn an_open_trench_is_not_a_void() {
        let mut g = grid();
        // air down from the top: open, not enclosed
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top - 50.0, top, Material::Air);
        }
        assert!(find_voids(&g).is_empty(), "an open trench is not a cavity");
    }

    #[test]
    fn a_one_column_gap_between_features_is_a_seam() {
        let mut g = grid();
        for k in 0..g.cols.len() {
            let ix = k % g.nx;
            if ix == 10 {
                continue; // the seam
            }
            let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top, top + 100.0, Material::SiliconDioxide);
        }
        let seam = seam_width_nm(&g);
        assert!(seam.is_some(), "a one-column gap between solid columns is a seam");
        assert!((seam.unwrap() - g.pitch).abs() < 1e-9);
    }

    #[test]
    fn a_fully_open_row_has_no_seam() {
        let g = grid();
        assert!(seam_width_nm(&g).is_none(), "nothing solid, no seam");
    }
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;

    fn wafer() -> Wafer { Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap() }

    fn flat_grid() -> ColumnGrid {
        let w = wafer();
        let mut g = ColumnGrid::sample(&w, &ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] }).unwrap();
        g.quantize_all();
        g
    }

    impl ColumnGrid {
        fn quantize_all(&mut self) {
            for k in 0..self.cols.len() {
                let mut c: Vec<Run> = self.cols[k]
                    .iter()
                    .map(|r| Run { z0: self.snap_z(r.z0), z1: self.snap_z(r.z1), material: r.material })
                    .collect();
                normalize(&mut c);
                self.cols[k] = c;
            }
        }
    }

    #[test]
    fn grid_matches_the_octree_pitch() {
        let g = flat_grid();
        assert_eq!(g.pitch, 16.0);
        assert_eq!(g.nx, 16);
        assert!(g.active_count() > 0);
    }

    #[test]
    fn directional_film_is_uniform() {
        let g = flat_grid();
        let (out, st) = grow_film(&g, Material::SiliconDioxide, Nm::new(100.0), 0.0).unwrap();
        assert!((st.mean_nm - 96.0).abs() < 1e-9, "100 nm quantises to 6 cells of 16 nm: {st:?}");
        assert!(st.nonuniformity_pct < 1e-9, "a flat surface gives a uniform film");
        let h0 = g.surface_heights()[0];
        let h1 = out.surface_heights()[0];
        assert!((h1 - h0 - 96.0).abs() < 1e-9);
    }

    #[test]
    fn conformal_film_fills_a_trench_more_than_directional() {
        let mut g = flat_grid();
        // cut a 64 nm wide trench 96 nm deep down the middle
        for iy in 0..g.ny {
            for ix in 7..9 {
                let k = g.idx(ix, iy);
                let top = solid_top(&g.cols[k]).unwrap().0;
                set_interval(&mut g.cols[k], top - 96.0, top, Material::Air);
            }
        }
        let (dir, sd) = grow_film(&g, Material::SiliconDioxide, Nm::new(48.0), 0.0).unwrap();
        let (con, sc) = grow_film(&g, Material::SiliconDioxide, Nm::new(48.0), 1.0).unwrap();
        let trench = g.idx(8, 8);
        let h = |gr: &ColumnGrid| solid_top(&gr.cols[trench]).unwrap().0;
        assert!(h(&con) > h(&dir), "conformal covers the sidewall: {} vs {}", h(&con), h(&dir));
        assert!(sc.max_nm >= sd.max_nm);
    }

    #[test]
    fn conformal_growth_is_monotone_in_coverage() {
        let mut g = flat_grid();
        for iy in 0..g.ny {
            for ix in 0..4 {
                let k = g.idx(ix, iy);
                let top = solid_top(&g.cols[k]).unwrap().0;
                set_interval(&mut g.cols[k], top - 96.0, top, Material::Air);
            }
        }
        let mut prev = 0.0;
        for s in [0.0, 0.25, 0.5, 0.75, 1.0] {
            let (o, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(64.0), s).unwrap();
            let top = solid_top(&o.cols[g.idx(0, 8)]).unwrap().0;
            assert!(top >= prev - 1e-9, "coverage {s} gave {top}, previous {prev}");
            prev = top;
        }
    }

    #[test]
    fn separable_dilation_tracks_the_circular_one_inside_the_reach() {
        // One tall column. Compare the separable L1 cone against the exact L2
        // cone that the earlier build used, at points strictly inside the lateral reach.
        // Outside the reach the two disagree by design: the L1 cone keeps
        // decaying where the circular kernel has already stopped, so the
        // comparison is only meaningful inside it.
        let mut g = flat_grid();
        let k0 = g.idx(8, 8);
        let base = g.surface_heights()[k0];
        set_interval(&mut g.cols[k0], base, base + 96.0, Material::Metal);
        let t = 32.0;
        let coverage = 1.0;
        let reach = coverage * t;
        let (o, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(t), coverage).unwrap();
        let h = g.surface_heights();
        for (ix, iy) in [(8usize, 8usize), (9, 8), (9, 9)] {
            let k = g.idx(ix, iy);
            let mut exact = h[k] + t;
            for jy in 0..g.ny {
                for jx in 0..g.nx {
                    let j = g.idx(jx, jy);
                    if !g.active[j] { continue; }
                    let dx = (jx as f64 - ix as f64) * g.pitch;
                    let dy = (jy as f64 - iy as f64) * g.pitch;
                    let d = (dx * dx + dy * dy).sqrt();
                    if d <= reach {
                        let v = h[j] + t * (1.0 - (d / reach).powi(2)).max(0.0).sqrt();
                        if v > exact { exact = v; }
                    }
                }
            }
            let got = solid_top(&o.cols[k]).unwrap().0;
            let diff = (got - exact).abs();
            let bound = t * (1.0 - 1.0 / 2f64.sqrt()) + g.pitch;
            assert!(diff <= bound, "at ({ix},{iy}) got {got:.1}, circular {exact:.1}, diff {diff:.1} > {bound:.1}");
        }
    }

    #[test]
    fn the_lateral_kernel_reaches_no_further_than_one_reach() {
        // A film on a flat surface must not be raised by a neighbour beyond the
        // reach; this pins the boundary behaviour of the separable kernel.
        let mut g = flat_grid();
        let k0 = g.idx(8, 8);
        let base = g.surface_heights()[k0];
        set_interval(&mut g.cols[k0], base, base + 96.0, Material::Metal);
        let (o, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(32.0), 1.0).unwrap();
        let far = g.idx(0, 0);
        let got = solid_top(&o.cols[far]).unwrap().0 - g.surface_heights()[far];
        assert!((got - 32.0).abs() < 1e-9, "a column far from the step gets the plain film thickness, got {got}");
    }

    #[test]
    fn etch_stops_on_a_zero_rate_layer() {
        let g = flat_grid();
        let (o1, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(96.0), 0.0).unwrap();
        let spec = EtchSpec { rates: EtchRates { si: 0.0, sio2: 60.0, si3n4: 0.0, resist: 8.0, metal: 0.0 }, time_min: 5.0 };
        let (o2, st) = etch(&o1, &spec, 1.0).unwrap();
        assert!(st.feature_depth_max_nm >= 90.0, "the oxide is cleared: {st:?}");
        let top = solid_top(&o2.cols[0]).unwrap();
        assert_eq!(top.1, Material::Silicon, "the etch stops on silicon at zero rate");
        assert!(st.feature_depth_max_nm < 100.0, "no silicon is removed: {st:?}");
    }

    #[test]
    fn mask_consumption_is_reported_separately_from_feature_depth() {
        let mut g = flat_grid();
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).unwrap().0;
            set_interval(&mut g.cols[k], top, top + 64.0, Material::Photoresist);
        }
        let spec = EtchSpec { rates: EtchRates { si: 10.0, sio2: 60.0, si3n4: 5.0, resist: 40.0, metal: 1.0 }, time_min: 3.0 };
        let (_, st) = etch(&g, &spec, 1.0).unwrap();
        assert!(st.mask_failed_columns > 0, "64 nm of resist at 40 nm/min is gone in 3 min: {st:?}");
        assert_eq!(st.feature_depth_max_nm, 0.0, "resist consumption must not be reported as feature depth");
        assert!(st.mask_consumed_mean_nm > 0.0);
    }

    #[test]
    fn a_void_does_not_shield_the_material_below_it() {
        let mut g = flat_grid();
        let k = g.idx(8, 8);
        let top = solid_top(&g.cols[k]).unwrap().0;
        set_interval(&mut g.cols[k], top - 48.0, top - 32.0, Material::Air);
        let spec = EtchSpec { rates: EtchRates { si: 60.0, sio2: 60.0, si3n4: 60.0, resist: 0.0, metal: 60.0 }, time_min: 1.0 };
        let (o, st) = etch(&g, &spec, 1.0).unwrap();
        assert!(st.feature_depth_max_nm > 60.0, "the etch passes the void: {st:?}");
        let new_top = solid_top(&o.cols[k]).unwrap().0;
        assert!(new_top < top - 60.0);
    }

    #[test]
    fn a_healthy_grid_passes_its_own_check() {
        let g = flat_grid();
        assert!(g.check().is_empty(), "{:?}", g.check());
    }

    #[test]
    fn a_broken_grid_is_reported_not_ignored() {
        let mut g = flat_grid();
        if let Some(run) = g.cols[0].first_mut() {
            std::mem::swap(&mut run.z0, &mut run.z1);
        }
        assert!(!g.check().is_empty(), "an inverted run must be reported");

        let mut h = flat_grid();
        h.active.pop();
        assert!(h.check().iter().any(|m| m.contains("active flags")), "{:?}", h.check());
    }

    #[test]
    fn strip_removes_only_resist() {
        let mut g = flat_grid();
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).unwrap().0;
            set_interval(&mut g.cols[k], top, top + 64.0, Material::Photoresist);
        }
        let (o, removed) = strip_resist(&g, 1.0).unwrap();
        assert!(removed > 0.0);
        assert!(o.cols.iter().all(|c| !c.iter().any(|r| r.material == Material::Photoresist)));
        assert!(o.cols[0].iter().any(|r| r.material == Material::Silicon));
    }
}
