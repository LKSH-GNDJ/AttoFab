// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! 2D layout-viewer math (no UI types): world (nm, y up) ↔ screen (px, y down), zoom about the cursor,
//! fit, and polygon culling with a hard cap so a huge layout never stalls a frame.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GdsView { pub center: [f64; 2], pub nm_per_px: f64 }

impl Default for GdsView { fn default() -> Self { Self { center: [0.0, 0.0], nm_per_px: 1000.0 } } }

const MIN_SCALE: f64 = 0.05;
const MAX_SCALE: f64 = 1.0e9;

impl GdsView {
    pub fn to_screen(&self, rect_center: [f32; 2], p: [f64; 2]) -> [f32; 2] {
        [rect_center[0] + ((p[0] - self.center[0]) / self.nm_per_px) as f32, rect_center[1] - ((p[1] - self.center[1]) / self.nm_per_px) as f32]
    }
    pub fn to_world(&self, rect_center: [f32; 2], s: [f32; 2]) -> [f64; 2] {
        [self.center[0] + (s[0] - rect_center[0]) as f64 * self.nm_per_px, self.center[1] - (s[1] - rect_center[1]) as f64 * self.nm_per_px]
    }
    /// Zoom by `factor` (< 1 zooms in) keeping the world point under `cursor` fixed.
    pub fn zoom_at(&mut self, rect_center: [f32; 2], cursor: [f32; 2], factor: f64) {
        if !factor.is_finite() || factor <= 0.0 { return; }
        let before = self.to_world(rect_center, cursor);
        self.nm_per_px = (self.nm_per_px * factor).clamp(MIN_SCALE, MAX_SCALE);
        let after = self.to_world(rect_center, cursor);
        self.center[0] += before[0] - after[0];
        self.center[1] += before[1] - after[1];
    }
    pub fn pan_px(&mut self, dx: f32, dy: f32) {
        self.center[0] -= dx as f64 * self.nm_per_px;
        self.center[1] += dy as f64 * self.nm_per_px;
    }
    /// Frame `bbox = [x0, y0, x1, y1]` (nm) in a viewport of `size` px with a 6 % margin.
    pub fn fit(&mut self, size: [f32; 2], bbox: [f64; 4]) {
        let (w, h) = ((bbox[2] - bbox[0]).max(1.0), (bbox[3] - bbox[1]).max(1.0));
        self.center = [(bbox[0] + bbox[2]) * 0.5, (bbox[1] + bbox[3]) * 0.5];
        self.nm_per_px = (w / size[0].max(1.0) as f64).max(h / size[1].max(1.0) as f64).mul_add(1.06, 0.0).clamp(MIN_SCALE, MAX_SCALE);
    }
    /// World rectangle visible in a viewport of `size` px.
    pub fn visible(&self, size: [f32; 2]) -> [f64; 4] {
        let (hw, hh) = (size[0] as f64 * 0.5 * self.nm_per_px, size[1] as f64 * 0.5 * self.nm_per_px);
        [self.center[0] - hw, self.center[1] - hh, self.center[0] + hw, self.center[1] + hh]
    }
}

/// Indices (into `bboxes`, chosen from `candidates`) of polygons worth drawing: they intersect `view`, are at
/// least `min_extent_nm` across, and — if more than `cap` qualify — only the `cap` largest are kept.
pub fn pick_visible(bboxes: &[[f64; 4]], candidates: &[u32], view: [f64; 4], min_extent_nm: f64, cap: usize) -> Vec<u32> {
    let mut v: Vec<(f64, u32)> = candidates.iter().filter_map(|&i| {
        let b = bboxes[i as usize];
        if b[2] < view[0] || b[0] > view[2] || b[3] < view[1] || b[1] > view[3] { return None; }
        let ext = (b[2] - b[0]).max(b[3] - b[1]);
        if ext < min_extent_nm { None } else { Some((ext, i)) }
    }).collect();
    if v.len() > cap {
        v.select_nth_unstable_by(cap, |a, b| b.0.partial_cmp(&a.0).unwrap());
        v.truncate(cap);
    }
    v.into_iter().map(|x| x.1).collect()
}

/// Smallest-dimension 10th percentile (nm) of the polygons of `candidates` whose bbox intersects `region`.
pub fn local_feature_p10(bboxes: &[[f64; 4]], candidates: &[u32], region: [f64; 4]) -> Option<f64> {
    let mut d: Vec<f64> = candidates.iter().filter_map(|&i| {
        let b = bboxes[i as usize];
        if b[2] < region[0] || b[0] > region[2] || b[3] < region[1] || b[1] > region[3] { None } else { Some((b[2] - b[0]).min(b[3] - b[1])) }
    }).collect();
    if d.is_empty() { return None; }
    d.sort_by(|a, b| a.partial_cmp(b).unwrap());
    Some(d[(d.len() as f64 * 0.10) as usize])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn world_screen_round_trip_has_y_up() {
        let v = GdsView { center: [1000.0, -500.0], nm_per_px: 20.0 };
        let rc = [400.0, 300.0];
        let p = [1500.0, -300.0];
        let s = v.to_screen(rc, p);
        assert!(s[0] > rc[0] && s[1] < rc[1], "+x is right, +y is up (screen y decreases)");
        let q = v.to_world(rc, s);
        assert!((q[0] - p[0]).abs() < 1e-3 && (q[1] - p[1]).abs() < 1e-3);
    }

    #[test]
    fn zoom_keeps_the_point_under_the_cursor_and_clamps() {
        let mut v = GdsView::default();
        let (rc, cur) = ([300.0, 200.0], [450.0, 120.0]);
        let w0 = v.to_world(rc, cur);
        for f in [0.5, 0.1, 3.0] { v.zoom_at(rc, cur, f); let w = v.to_world(rc, cur); assert!((w[0] - w0[0]).abs() < 1e-6 * v.nm_per_px.max(1.0) && (w[1] - w0[1]).abs() < 1e-6 * v.nm_per_px.max(1.0), "factor {f}"); }
        for _ in 0..200 { v.zoom_at(rc, cur, 0.1); } assert!(v.nm_per_px >= MIN_SCALE);
        for _ in 0..200 { v.zoom_at(rc, cur, 10.0); } assert!(v.nm_per_px <= MAX_SCALE);
        let before = v; v.zoom_at(rc, cur, f64::NAN); assert_eq!(v, before, "bad factors are ignored");
    }

    #[test]
    fn fit_contains_the_box_and_pan_follows_the_drag() {
        let mut v = GdsView::default();
        let (size, bbox) = ([800.0, 500.0], [0.0, 0.0, 9.1e7, 8.5e7]);
        v.fit(size, bbox);
        let vis = v.visible(size);
        assert!(vis[0] <= bbox[0] && vis[1] <= bbox[1] && vis[2] >= bbox[2] && vis[3] >= bbox[3]);
        let rc = [400.0, 250.0];
        let p = v.to_world(rc, [500.0, 250.0]);
        v.pan_px(30.0, -10.0);
        let s = v.to_screen(rc, p);
        assert!((s[0] - 530.0).abs() < 1e-3 && (s[1] - 240.0).abs() < 1e-3, "content moves with the pointer: {s:?}");
    }

    #[test]
    fn culling_respects_view_size_and_cap() {
        let bb: Vec<[f64; 4]> = (0..1000).map(|i| { let x = i as f64 * 100.0; [x, 0.0, x + 10.0 + (i % 50) as f64, 20.0] }).collect();
        let all: Vec<u32> = (0..1000).collect();
        let vis = pick_visible(&bb, &all, [0.0, -10.0, 5000.0, 30.0], 0.0, 10_000);
        assert!(vis.len() == 51 && vis.iter().all(|&i| bb[i as usize][0] <= 5000.0), "{}", vis.len());
        assert!(pick_visible(&bb, &all, [0.0, -10.0, 5000.0, 30.0], 40.0, 10_000).iter().all(|&i| bb[i as usize][2] - bb[i as usize][0] >= 40.0));
        let capped = pick_visible(&bb, &all, [-1.0e9, -1.0e9, 1.0e9, 1.0e9], 0.0, 100);
        assert_eq!(capped.len(), 100);
        let min_kept = capped.iter().map(|&i| bb[i as usize][2] - bb[i as usize][0]).fold(f64::INFINITY, f64::min);
        assert!(min_kept >= 50.0, "over the cap the largest shapes survive: {min_kept}");
        assert_eq!(local_feature_p10(&bb, &all, [0.0, 0.0, 1050.0, 10.0]).map(|x| x >= 10.0), Some(true));
        assert!(local_feature_p10(&bb, &all, [9.0e9, 9.0e9, 9.1e9, 9.1e9]).is_none());
    }
}
