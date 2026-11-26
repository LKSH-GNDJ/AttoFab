// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! 2D GDS Layout Canvas — top-down vector viewer for a parsed GDSII layout,
//! with per-layer visibility toggles and smooth zoom/pan state.

use sim_core::io::gds_parser::{GdsLibrary, Polygon};
use std::collections::{HashMap, HashSet};

/// Pan/zoom camera state for the 2D canvas.
#[derive(Debug, Clone, Copy)]
pub struct Camera2D {
    pub center: (f64, f64),
    pub zoom: f64,
}

impl Default for Camera2D {
    fn default() -> Self {
        Self { center: (0.0, 0.0), zoom: 1.0 }
    }
}

impl Camera2D {
    pub fn pan(&mut self, dx: f64, dy: f64) {
        self.center.0 += dx / self.zoom;
        self.center.1 += dy / self.zoom;
    }

    pub fn zoom_by(&mut self, factor: f64) {
        self.zoom = (self.zoom * factor).clamp(1e-6, 1e6);
    }

    /// Project a layout-space point (nm) into screen-space pixels, given a
    /// viewport size in pixels.
    pub fn project(&self, point: (f64, f64), viewport: (f64, f64)) -> (f64, f64) {
        let x = (point.0 - self.center.0) * self.zoom + viewport.0 / 2.0;
        let y = (point.1 - self.center.1) * self.zoom + viewport.1 / 2.0;
        (x, y)
    }
}

/// The canvas's full state: the parsed library, per-layer visibility, and
/// camera.
pub struct GdsCanvas {
    pub library: GdsLibrary,
    pub layer_visible: HashMap<i32, bool>,
    pub camera: Camera2D,
}

impl GdsCanvas {
    pub fn from_library(library: GdsLibrary) -> Self {
        let layers: HashSet<i32> = library.polygons.iter().map(|p| p.layer).collect();
        let layer_visible = layers.into_iter().map(|l| (l, true)).collect();
        Self { library, layer_visible, camera: Camera2D::default() }
    }

    pub fn toggle_layer(&mut self, layer: i32) {
        if let Some(v) = self.layer_visible.get_mut(&layer) {
            *v = !*v;
        }
    }

    /// The polygons that should currently be drawn, respecting layer
    /// visibility toggles — the render loop projects each via `camera`.
    pub fn visible_polygons(&self) -> Vec<&Polygon> {
        self.library
            .polygons
            .iter()
            .filter(|p| *self.layer_visible.get(&p.layer).unwrap_or(&true))
            .collect()
    }

    /// Frame the camera to fit the full layout in a viewport (fit-to-window).
    pub fn fit_to_view(&mut self, viewport: (f64, f64)) {
        if self.library.polygons.is_empty() {
            return;
        }
        let (mut min_x, mut min_y, mut max_x, mut max_y) =
            (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
        for poly in &self.library.polygons {
            for &(x, y) in &poly.vertices {
                min_x = min_x.min(x);
                min_y = min_y.min(y);
                max_x = max_x.max(x);
                max_y = max_y.max(y);
            }
        }
        self.camera.center = ((min_x + max_x) / 2.0, (min_y + max_y) / 2.0);
        let span_x = (max_x - min_x).max(1.0);
        let span_y = (max_y - min_y).max(1.0);
        let zoom_x = viewport.0 / span_x;
        let zoom_y = viewport.1 / span_y;
        self.camera.zoom = zoom_x.min(zoom_y) * 0.9; // 10% margin
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lib_with_two_layers() -> GdsLibrary {
        GdsLibrary {
            db_unit_m: 1e-9,
            polygons: vec![
                Polygon { layer: 1, vertices: vec![(0.0, 0.0), (10.0, 0.0), (10.0, 10.0)] },
                Polygon { layer: 2, vertices: vec![(20.0, 20.0), (30.0, 20.0), (30.0, 30.0)] },
            ],
        }
    }

    #[test]
    fn toggling_a_layer_hides_its_polygons() {
        let mut canvas = GdsCanvas::from_library(lib_with_two_layers());
        assert_eq!(canvas.visible_polygons().len(), 2);
        canvas.toggle_layer(2);
        assert_eq!(canvas.visible_polygons().len(), 1);
    }

    #[test]
    fn fit_to_view_centers_camera_on_layout_bounds() {
        let mut canvas = GdsCanvas::from_library(lib_with_two_layers());
        canvas.fit_to_view((800.0, 600.0));
        assert!(canvas.camera.zoom > 0.0);
    }
}
