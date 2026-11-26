// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! GDSII workflow on top of the loader: layout viewer window (click to choose where to simulate),
//! local feature size, the simulation-scale panel and "auto-scale & use layer".

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, RichText, Sense, Stroke, Vec2};
use sim_core::params::StageParams;

use crate::app::AttoFabApp;
use crate::gds_view::{local_feature_p10, pick_visible};

const DRAW_CAP: usize = 30_000;

fn layer_color(layer: i32, datatype: i32) -> Color32 {
    // Deterministic, well-separated hues for non-selected layers.
    let h = ((layer * 47 + datatype * 23) % 360) as f32 / 360.0;
    egui::ecolor::Hsva::new(h, 0.55, 0.85, 1.0).into()
}

impl AttoFabApp {
    /// Thinnest film or resist the recipe needs to resolve vertically (nm), infinity if none.
    pub(crate) fn thinnest_film_nm(&self) -> f64 {
        let mut t = f64::INFINITY;
        for s in self.session.state.recipe.stages.iter().filter(|s| s.enabled) {
            match &s.params {
                StageParams::Deposition(d) => t = t.min(d.thickness_nm),
                StageParams::Lithography(l) => t = t.min(l.resist_thickness_nm),
                _ => {}
            }
        }
        t
    }

    /// 10th-percentile smallest dimension (nm) of the selected layer near the crop centre
    /// (falls back to the whole layer when nothing is nearby).
    pub(crate) fn local_feature_nm(&self) -> Option<f64> {
        let g = self.gds.as_ref()?;
        let l = g.layers.get(self.gds_sel)?;
        let idx = g.by_layer.get(&(l.layer, l.datatype))?;
        let r = (self.session.max_window_nm() * 2.0).max(6000.0);
        let region = [self.gds_crop[0] - r, self.gds_crop[1] - r, self.gds_crop[0] + r, self.gds_crop[1] + r];
        local_feature_p10(&g.bboxes, idx, region).or(Some(l.feature_p10_nm))
    }

    /// Move the crop centre onto the smallest feature of the selected layer (where fine cells are needed).
    pub fn jump_to_finest_feature(&mut self) -> bool {
        let Some(g) = &self.gds else { return false };
        let Some(l) = g.layers.get(self.gds_sel) else { return false };
        let Some(idx) = g.by_layer.get(&(l.layer, l.datatype)) else { return false };
        let dim = |i: &u32| { let b = g.bboxes[*i as usize]; (b[2] - b[0]).min(b[3] - b[1]) };
        let Some(&best) = idx.iter().filter(|i| dim(i) > 0.0).min_by(|a, b| dim(a).partial_cmp(&dim(b)).unwrap()) else { return false };
        let b = g.bboxes[best as usize];
        self.gds_crop = [(b[0] + b[2]) * 0.5, (b[1] + b[3]) * 0.5];
        self.gds_msg = Some((true, format!("Crop moved to the finest feature of layer {}/{} ({:.0} nm across).", l.layer, l.datatype, dim(&best))));
        true
    }


    /// Move the crop centre onto the nearest polygon of the selected layer.
    ///
    /// "Finest feature" jumps to the smallest shape anywhere on the layer, which can be
    /// metres away from what the user is looking at. This one stays close: it picks the
    /// polygon whose centre is nearest the current crop, so a small nudge lands on the
    /// pattern instead of leaving the window empty.
    pub fn jump_to_nearest_feature(&mut self) -> bool {
        let Some(g) = &self.gds else { return false };
        let Some(l) = g.layers.get(self.gds_sel) else { return false };
        let Some(idx) = g.by_layer.get(&(l.layer, l.datatype)) else { return false };
        let (cx, cy) = (self.gds_crop[0], self.gds_crop[1]);
        let dist = |i: &u32| {
            let b = g.bboxes[*i as usize];
            let dx = (b[0] + b[2]) * 0.5 - cx;
            let dy = (b[1] + b[3]) * 0.5 - cy;
            dx * dx + dy * dy
        };
        let Some(&best) = idx.iter().min_by(|a, b| dist(a).partial_cmp(&dist(b)).unwrap_or(std::cmp::Ordering::Equal)) else { return false };
        let b = g.bboxes[best as usize];
        self.gds_crop = [(b[0] + b[2]) * 0.5, (b[1] + b[3]) * 0.5];
        self.gds_msg = Some((true, format!(
            "Crop moved to the nearest polygon of layer {}/{} ({:.0} nm across).",
            l.layer, l.datatype, (b[2] - b[0]).min(b[3] - b[1]))));
        true
    }


    /// Rebuild a GDS-derived mask after the crop moved.
    ///
    /// The mask holds the polygons around the *old* crop centre. Move the crop without
    /// rebuilding and the new window is still covered by stale geometry: lithography opens
    /// no resist, the etch removes nothing, and the implant finds no silicon. One stale
    /// mask makes every later stage fail for a reason that looks like their fault.
    ///
    /// Only a mask that already came from this layout is rebuilt, and only on a Lithography
    /// stage — a recipe using a line/space mask is left alone (`None`: nothing to do).
    /// Once a rebuild is attempted the outcome is always reported, in the return value and
    /// in `gds_msg`: a rebuild that cannot happen (no geometry under the new crop, an
    /// invalid mask) must say so, because the previous mask is then still in the recipe
    /// and still describes the old crop (PRD FR-6.3, defect D-1).
    pub(crate) fn rebuild_mask_for_crop(&mut self) -> Option<String> {
        let prep: Result<(_, i32, i32, bool, f64, bool, u32), String> = {
            let g = self.gds.as_ref()?;
            let l = g.layers.get(self.gds_sel)?.clone();
            let (idx, sid) = self.session.state.recipe.stages.iter().enumerate().find_map(|(i, s)| {
                let ok = s.enabled && matches!(&s.params, StageParams::Lithography(lp) if matches!(lp.mask, sim_core::geom::MaskPattern::Polygons { .. }));
                if ok { Some((i, s.id)) } else { None }
            })?;
            let clear_field = match &self.session.state.recipe.stages[idx].params {
                StageParams::Lithography(lp) => match &lp.mask {
                    sim_core::geom::MaskPattern::Polygons { clear_field, .. } => *clear_field,
                    _ => true,
                },
                _ => true,
            };
            let half = self.session.state.recipe.window.size[0].max(self.session.state.recipe.window.size[1]) * 0.5 + 500.0;
            let region = [self.gds_crop[0] - half, self.gds_crop[1] - half, self.gds_crop[0] + half, self.gds_crop[1] + half];
            match g.layout.flatten_region(Some(&g.top), Some(region)) {
                Ok(near) => Ok((near, l.layer, l.datatype, self.gds_all_layers, half, clear_field, sid)),
                Err(e) => Err(format!("cannot read the layout around the new crop: {e}")),
            }
        };
        let (near, layer, datatype, all, half, clear_field, stage_id) = match prep {
            Ok(t) => t,
            Err(why) => return Some(self.mask_not_rebuilt(&why)),
        };
        let built = if all {
            sim_core::geom::MaskPattern::from_layout_all(&near, (self.gds_crop[0], self.gds_crop[1]), half, clear_field)
        } else {
            sim_core::geom::MaskPattern::from_layout(&near, layer, Some(datatype), (self.gds_crop[0], self.gds_crop[1]), half, clear_field)
        };
        let (m, kept) = match built {
            Ok(t) => t,
            Err(e) => return Some(self.mask_not_rebuilt(&e.to_string())),
        };
        let mut p = match self.session.state.recipe.stages.iter().find(|s| s.id == stage_id) {
            Some(s) => s.params.clone(),
            None => return Some(self.mask_not_rebuilt("the lithography stage no longer exists")),
        };
        if let StageParams::Lithography(lp) = &mut p {
            lp.mask = m;
        }
        if let Err(e) = self.session.edit_params(stage_id, p) {
            return Some(self.mask_not_rebuilt(&e.to_string()));
        }
        let stale = self.session.status_of(stage_id) == Some(sim_core::session::StageStatus::Stale);
        let note = format!(
            " Mask rebuilt: {kept} polygon(s) around the new crop.{}",
            if stale { " Lithography and every stage after it are now stale: press Run." } else { "" }
        );
        self.gds_msg = Some((true, note.trim().to_string()));
        Some(note)
    }

    /// Report a mask rebuild that did not happen. The previous mask is still in the recipe.
    fn mask_not_rebuilt(&mut self, why: &str) -> String {
        let note = format!(" Mask NOT rebuilt: {why}. The recipe still holds the mask for the previous crop, so lithography will not match this window.");
        self.gds_msg = Some((false, note.trim().to_string()));
        note
    }

    /// Choose cell size, tile depth and window from the local feature size and the thinnest film, apply them,
    /// then load the selected layer as the Lithography mask.
    pub fn auto_scale_from_gds(&mut self) {
        let Some(feature) = self.local_feature_nm() else { self.gds_msg = Some((false, "Load a GDSII file and select a layer first.".into())); return };
        let advice = sim_core::scale::advise(feature, self.thinnest_film_nm());
        if let Err(e) = self.session.set_scale(advice.cell_nm, advice.tile_levels) { self.gds_msg = Some((false, e.to_string())); return; }
        let window = (128.0 * advice.cell_nm).min(self.session.max_window_nm());
        let _ = self.session.set_window_size([window, window]);
        self.scale_edit = (self.session.state.spec.leaf_nm, self.session.state.spec.tile_levels);
        self.use_gds_layer();
        self.need_focus = true;
        let mut msg = format!("Scale: cell {:.1} nm, tile {:.2} µm, window {:.2} µm (from {:.0} nm features, {:.0} nm thinnest film).", advice.cell_nm, self.session.state.spec.tile_edge_nm() / 1e3, window / 1e3, feature, self.thinnest_film_nm());
        for n in &advice.notes { msg.push_str(&format!(" Note: {n}.")); }
        if let Some((true, m)) = &self.gds_msg { msg = format!("{msg} {m}"); }
        self.gds_msg = Some((true, msg));
    }

    pub(crate) fn scale_panel(&mut self, ui: &mut egui::Ui) {
        ui.label(RichText::new("Simulation scale").strong());
        let spec = self.session.state.spec;
        ui.horizontal(|ui| {
            ui.label("cell");
            ui.add(egui::DragValue::new(&mut self.scale_edit.0).speed(0.5).clamp_range(1.0..=2000.0).suffix(" nm"));
            ui.label("levels");
            let mut lv = self.scale_edit.1 as f64;
            if ui.add(egui::DragValue::new(&mut lv).speed(0.1).clamp_range(4.0..=12.0)).changed() { self.scale_edit.1 = lv.round() as u32; }
        });
        let tile = self.scale_edit.0 * (1u64 << self.scale_edit.1.min(30)) as f64;
        ui.label(RichText::new(format!("tile {:.2} µm · max window {:.2} µm (≤ {} cells)", tile / 1e3, (tile * 0.95).min(sim_core::geom::MAX_GRID_AXIS as f64 * self.scale_edit.0) / 1e3, sim_core::geom::MAX_GRID_AXIS)).weak().small());
        if (self.scale_edit.0, self.scale_edit.1) != (spec.leaf_nm, spec.tile_levels) && ui.button("Apply scale").on_hover_text("Rebuilds the blank tile; every committed stage becomes stale").clicked() {
            match self.session.set_scale(self.scale_edit.0, self.scale_edit.1) { Ok(()) => { self.need_focus = true; self.status = "Scale applied. Results are stale until you Run.".into(); } Err(e) => self.status = e.to_string() }
        }
        let film = self.thinnest_film_nm();
        if film.is_finite() && film / spec.leaf_nm < 2.0 { ui.colored_label(Color32::from_rgb(240, 180, 60), format!("⚠ the thinnest film ({film:.0} nm) spans {:.1} cells", film / spec.leaf_nm)); }
    }

    /// The layout viewer: all geometry of the loaded layout, the selected layer highlighted, click to set the crop centre.
    pub(crate) fn gds_viewer(&mut self, ctx: &egui::Context) {
        if !self.gds_view_open || self.gds.is_none() { return; }
        let mut open = true;
        let mut clicked_world: Option<[f64; 2]> = None;
        let mut fit = false;
        let mut zoom_crop = false;
        let mut finest = false;
        egui::Window::new("GDS layout — click to choose where to simulate").open(&mut open).default_size([760.0, 560.0]).resizable(true).show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.checkbox(&mut self.gds_show_all, "show all layers");
                if ui.button("Fit layout").clicked() { fit = true; }
                if ui.button("Zoom to crop").on_hover_text("Frame the process window").clicked() { zoom_crop = true; }
                if ui.button("Finest feature").on_hover_text("Move the crop onto the smallest shape of the selected layer").clicked() { finest = true; zoom_crop = true; }
                ui.label(RichText::new(format!("crop centre ({:.0}, {:.0}) nm", self.gds_crop[0], self.gds_crop[1])).monospace().small());
            });
            let size = ui.available_size().max(Vec2::new(200.0, 160.0));
            let (rect, resp) = ui.allocate_exact_size(size, Sense::click_and_drag());
            self.gds_canvas_rect = rect;
            let painter = ui.painter_at(rect);
            painter.rect_filled(rect, 0.0, Color32::from_rgb(12, 16, 22));
            let Some(g) = &self.gds else { return };
            let rc = [rect.center().x, rect.center().y];
            let sz = [rect.width(), rect.height()];
            if fit || self.gds_view_fit {
                let l = &g.layers[self.gds_sel.min(g.layers.len() - 1)];
                let b = if self.gds_show_all { g.layers.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |a, l| [a[0].min(l.bbox[0]), a[1].min(l.bbox[1]), a[2].max(l.bbox[2]), a[3].max(l.bbox[3])]) } else { l.bbox };
                self.gds_view.fit(sz, b);
                self.gds_view_fit = false;
            }
            if zoom_crop {
                let w = self.session.state.recipe.window.size[0].max(200.0);
                self.gds_view.center = self.gds_crop;
                self.gds_view.nm_per_px = (w * 3.0 / sz[0] as f64).max(0.05);
            }
            // Pan / zoom / click.
            if resp.dragged() { let d = resp.drag_delta(); self.gds_view.pan_px(d.x, d.y); }
            if let Some(p) = resp.hover_pos() {
                let scroll = ui.input(|i| i.smooth_scroll_delta.y);
                if scroll != 0.0 { self.gds_view.zoom_at(rc, [p.x, p.y], (-(scroll as f64) * 0.0025).exp()); }
                let w = self.gds_view.to_world(rc, [p.x, p.y]);
                painter.text(rect.left_bottom() + Vec2::new(6.0, -6.0), Align2::LEFT_BOTTOM, format!("({:.0}, {:.0}) nm", w[0], w[1]), FontId::monospace(11.0), Color32::from_gray(200));
            }
            if resp.clicked() { if let Some(p) = resp.interact_pointer_pos() { clicked_world = Some(self.gds_view.to_world(rc, [p.x, p.y])); } }

            let view = self.gds_view.visible(sz);
            let min_ext = self.gds_view.nm_per_px * 0.6;
            let draw = |idx: &[u32], c: Color32, fill: bool| {
                let to = |p: (f64, f64)| { let s = self.gds_view.to_screen(rc, [p.0, p.1]); Pos2::new(s[0], s[1]) };
                for &i in idx {
                    let poly = &g.flat[i as usize];
                    let b = g.bboxes[i as usize];
                    let rect_like = poly.pts.len() == 4 && poly.pts.iter().all(|v| (v.0 == b[0] || v.0 == b[2]) && (v.1 == b[1] || v.1 == b[3]));
                    if rect_like && fill {
                        let (a, d) = (to((b[0], b[1])), to((b[2], b[3])));
                        painter.rect_filled(Rect::from_two_pos(a, d), 0.0, c.gamma_multiply(0.55));
                        painter.rect_stroke(Rect::from_two_pos(a, d), 0.0, Stroke::new(1.0_f32, c));
                    } else {
                        let mut pts: Vec<Pos2> = poly.pts.iter().map(|v| to(*v)).collect();
                        pts.push(pts[0]);
                        painter.add(egui::Shape::line(pts, Stroke::new(1.0_f32, c)));
                    }
                }
            };
            if self.gds_show_all {
                for l in &g.layers {
                    if let Some(idx) = g.by_layer.get(&(l.layer, l.datatype)) {
                        if (l.layer, l.datatype) == (g.layers[self.gds_sel].layer, g.layers[self.gds_sel].datatype) { continue; }
                        let v = pick_visible(&g.bboxes, idx, view, min_ext, DRAW_CAP / 3);
                        draw(&v, layer_color(l.layer, l.datatype).gamma_multiply(0.45), false);
                    }
                }
            }
            if let Some(l) = g.layers.get(self.gds_sel) {
                if let Some(idx) = g.by_layer.get(&(l.layer, l.datatype)) {
                    let v = pick_visible(&g.bboxes, idx, view, min_ext * 0.5, DRAW_CAP);
                    draw(&v, Color32::from_rgb(255, 170, 60), true);
                }
            }
            // Crop box (process window) and its tile.
            let w = self.session.state.recipe.window.size;
            let tile = self.session.state.spec.tile_edge_nm();
            for (half, c, dashed) in [([tile * 0.5; 2], Color32::from_rgb(90, 160, 255), true), ([w[0] * 0.5, w[1] * 0.5], Color32::from_rgb(255, 80, 80), false)] {
                let (a, b) = (self.gds_view.to_screen(rc, [self.gds_crop[0] - half[0], self.gds_crop[1] - half[1]]), self.gds_view.to_screen(rc, [self.gds_crop[0] + half[0], self.gds_crop[1] + half[1]]));
                let r = Rect::from_two_pos(Pos2::new(a[0], a[1]), Pos2::new(b[0], b[1]));
                if dashed { painter.extend(egui::Shape::dashed_line(&[r.left_top(), r.right_top(), r.right_bottom(), r.left_bottom(), r.left_top()], Stroke::new(1.0_f32, c), 6.0, 4.0)); } else { painter.rect_stroke(r.expand(0.5), 0.0, Stroke::new(1.5_f32, c)); }
                // Always visible even when tiny.
                let ctr = self.gds_view.to_screen(rc, self.gds_crop);
                painter.circle_stroke(Pos2::new(ctr[0], ctr[1]), 6.0, Stroke::new(1.5_f32, Color32::from_rgb(255, 80, 80)));
            }
            painter.text(rect.right_top() + Vec2::new(-8.0, 6.0), Align2::RIGHT_TOP, "red: process window   blue: simulated tile   gold: selected layer", FontId::proportional(10.5), Color32::from_gray(170));
        });
        if finest { self.jump_to_finest_feature(); self.gds_view_fit = false; let w = self.session.state.recipe.window.size[0].max(200.0); self.gds_view.center = self.gds_crop; self.gds_view.nm_per_px = (w * 3.0 / self.gds_canvas_rect.width().max(1.0) as f64).max(0.05); }
        if let Some(w) = clicked_world { self.gds_crop = w; let _ = self.rebuild_mask_for_crop(); }
        self.gds_view_open = open;
    }

}
