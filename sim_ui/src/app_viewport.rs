//! Viewport: input routing (explicit tool modes), GPU paint, and 2D overlays.

use eframe::egui::{self, Align2, Color32, FontId, Pos2, Rect, Sense, Stroke, Vec2};
use sim_core::legend::{layers_for_stage, LegendLayer};
use sim_core::pipeline::{Overlay, StepKind};
use sim_core::transform::Vec3;
use sim_core::wafer::{Bounds, Material, Run};

use crate::app::{col, AttoFabApp, DisplayKind, InspectReadout, Tool};
use crate::svo_view::SvoViewportCallback;
use crate::viewport_math::*;

fn p2(v: [f32; 2]) -> Pos2 { Pos2::new(v[0], v[1]) }

impl AttoFabApp {
    pub(crate) fn work_plane(&self) -> WorkPlane {
        let slice_z = self.surface_z() + self.slice_offset_nm;
        WorkPlane::for_view(self.rig.mode, &self.transform(), self.surface_z(), slice_z, self.site())
    }

    #[cfg(test)]
    pub(crate) fn work_plane_for_test(&self) -> WorkPlane { self.drag_plane() }

    /// Plane through the wafer origin used while dragging the wafer.
    fn drag_plane(&self) -> WorkPlane { WorkPlane::for_view(self.rig.mode, &self.transform(), self.surface_z(), self.surface_z(), [0.0, 0.0]) }

    /// Section plane in wafer-local coordinates `(n, d)`: keep `dot(n, p) + d >= 0`.
    fn section_clip_local(&self) -> Option<([f64; 3], f64)> {
        if !self.rig.mode.is_section() { return None; }
        let t = self.transform();
        let n = t.dir_to_local(self.rig.forward());
        let p0 = t.world_to_local(self.work_plane().point);
        Some((n, -dot(n, p0)))
    }

    fn rect_arrays(rect: Rect) -> ([f32; 2], [f32; 2]) { ([rect.min.x, rect.min.y], [rect.width(), rect.height()]) }

    pub(crate) fn project_world(&self, rect: Rect, w: Vec3) -> Option<Pos2> {
        let n = self.rig.project(w, self.viewport_aspect)?;
        let (mn, sz) = Self::rect_arrays(rect);
        Some(p2(ndc_to_pointer(n, mn, sz)))
    }
    /// Where a wafer-local point lands, in pixels, for a viewport whose top-left is the origin and which is `w × h` px
    /// (the same projection the overlays use). Lets tests tie overlays to rendered pixels.
    pub fn project_local_px(&self, l: [f64; 3], w: f32, h: f32) -> Option<[f32; 2]> {
        self.project_local(Rect::from_min_size(Pos2::ZERO, Vec2::new(w, h)), l).map(|p| [p.x, p.y])
    }
    pub(crate) fn project_local(&self, rect: Rect, l: Vec3) -> Option<Pos2> { self.project_world(rect, self.transform().local_to_world(l)) }

    pub(crate) fn viewport(&mut self, ui: &mut egui::Ui) {
        let size = ui.available_size();
        let (rect, resp) = ui.allocate_exact_size(size, Sense::click_and_drag());
        self.viewport_aspect = (rect.width() / rect.height().max(1.0)) as f64;
        self.viewport_rect = rect;
        if self.need_focus && rect.width() > 60.0 { self.focus_site(); self.need_focus = false; }
        self.handle_input(ui.ctx(), &resp, rect);
        let painter = ui.painter_at(rect);
        painter.rect_filled(rect, 0.0, col(self.palette.background));
        if self.gpu_ok {
            let u = pack_uniform(&self.build_uniform(rect));
            SvoViewportCallback { nodes: self.display.nodes.clone(), revision: self.display.revision, uniform: u, rect }.add_to(ui, rect);
        } else {
            painter.text(rect.center(), Align2::CENTER_CENTER, "GPU renderer unavailable (wgpu not initialised). Overlays and analytics still work.", FontId::proportional(14.0), Color32::GRAY);
        }
        self.draw_overlays(&painter, rect);
        if resp.hovered() { ui.ctx().set_cursor_icon(match self.tool { Tool::MoveWafer => egui::CursorIcon::Grab, Tool::Place | Tool::Inspect | Tool::Select => egui::CursorIcon::Crosshair, Tool::Navigate => egui::CursorIcon::Default }); }
    }

    /// Nodes + shader uniform for an offscreen frame of `w × h` pixels (headless rendering and regression tests).
    pub fn offscreen_inputs(&mut self, w: f32, h: f32) -> (std::sync::Arc<Vec<sim_core::render_wgpu::GpuSvoNode>>, [[f32; 4]; UNIFORM_VEC4S]) {
        let rect = Rect::from_min_size(Pos2::ZERO, Vec2::new(w, h));
        self.viewport_aspect = (w / h) as f64;
        (self.display.nodes.clone(), pack_uniform(&self.build_uniform(rect)))
    }

    // ----- uniform -------------------------------------------------------------------------

    fn build_uniform(&self, rect: Rect) -> UniformInput {
        let t = self.transform();
        let w = &self.display.wafer;
        let o = w.tile_origin();
        let origin = [o.0, o.1, o.2];
        let rel = |l: Vec3| [(l[0] - origin[0]) as f32, (l[1] - origin[1]) as f32, (l[2] - origin[2]) as f32];
        let (f, r, u) = self.rig.basis();
        let eye_l = t.world_to_local(self.rig.eye());
        let clip = self.section_clip_local().map(|(n, d)| {
            let d_rel = d + dot(n, origin); // plane offset expressed in tile-relative coordinates
            [n[0] as f32, n[1] as f32, n[2] as f32, d_rel as f32]
        });
        let (pal, vis) = self.palette.shader_tables();
        UniformInput {
            eye: rel(eye_l), forward: rel_dir(t.dir_to_local(f)), right: rel_dir(t.dir_to_local(r)), up: rel_dir(t.dir_to_local(u)),
            viewport_px: [rect.width(), rect.height()], origin_px: [0.0, 0.0], mode: self.rig.mode, slice_z: (self.slice_offset_nm) as f32, clip,
            domain_radius: w.domain.radius_nm as f32, domain_center: [(w.domain.center_x_nm - origin[0]) as f32, (w.domain.center_y_nm - origin[1]) as f32],
            bulk_z: if w.implicit_bulk { Some([(0.0 - origin[2]) as f32, (w.surface_z_nm - origin[2]) as f32]) } else { None },
            tile_min: rel([w.bounds.min.0, w.bounds.min.1, w.bounds.min.2]), tile_max: rel([w.bounds.max.0, w.bounds.max.1, w.bounds.max.2]),
            selection: self.selection.as_ref().map(|(b, _)| (rel([b.min.0, b.min.1, b.min.2]), rel([b.max.0, b.max.1, b.max.2]))),
            palette: pal, visible: vis,
        }
    }

    // ----- input ------------------------------------------------------------------------------

    fn handle_input(&mut self, ctx: &egui::Context, resp: &egui::Response, rect: Rect) {
        let (mn, sz) = Self::rect_arrays(rect);
        let to_ndc = |p: Pos2| pointer_to_ndc([p.x, p.y], mn, sz);
        let aspect = self.viewport_aspect;
        let t = self.transform();
        let radius = self.wafer_radius();
        let plane = self.work_plane();

        // A wafer drag ends the moment the button is released anywhere, even outside the widget.
        if self.mouse.wafer_drag.is_some() && !ctx.input(|i| i.pointer.primary_down()) { self.session.end_drag(); self.mouse.wafer_drag = None; }

        self.hover = resp.hover_pos().and_then(to_ndc).and_then(|n| pick_on_wafer(&self.rig, n, aspect, &plane, &t, radius));

        if resp.hovered() {
            let (scroll, pinch) = ctx.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
            let mut k = 1.0f64;
            if scroll != 0.0 { k *= (-(scroll as f64) * 0.0025).exp(); }
            if pinch != 1.0 { k /= pinch as f64; }
            if (k - 1.0).abs() > 1e-9 {
                let anchor = self.hover.map(|h| h.world).unwrap_or(self.rig.target);
                self.rig.zoom_about(anchor, k);
            }
        }

        let press = ctx.input(|i| i.pointer.press_origin()).and_then(to_ndc);
        let cur = resp.interact_pointer_pos().and_then(to_ndc);
        let primary = resp.dragged_by(egui::PointerButton::Primary);
        let middle = resp.dragged_by(egui::PointerButton::Middle);
        let secondary = resp.dragged_by(egui::PointerButton::Secondary);
        let persp = self.rig.mode == ViewMode::Perspective;

        // Camera: orbit (3D) / pan. Never touches the wafer transform.
        let pan_active = middle || (self.tool == Tool::Navigate && !persp && primary);
        if pan_active {
            if self.mouse.pan.is_none() { self.mouse.pan = press.and_then(|n| PanDrag::begin(&self.rig, n, aspect, plane)); }
            if let (Some(p), Some(c)) = (&self.mouse.pan, cur) { if let Some(r) = p.update(c, aspect) { self.rig = r; } }
        } else { self.mouse.pan = None; }
        if secondary || (self.tool == Tool::Navigate && persp && primary) {
            let d = resp.drag_delta();
            self.rig.orbit(d.x as f64, d.y as f64);
        }

        // Move Wafer: explicit mode, drag-captured, locked while a process animation plays.
        if self.tool == Tool::MoveWafer && primary {
            if self.locked() {
                self.status = "Wafer movement is locked while a process animation is playing. Pause it first.".into();
            } else {
                if self.mouse.wafer_drag.is_none() {
                    if let Some(n) = press {
                        let dp = self.drag_plane();
                        self.mouse.wafer_drag = WaferDrag::begin(&self.rig, n, aspect, dp, t);
                        if self.mouse.wafer_drag.is_some() { self.session.begin_drag(); }
                    }
                }
                if let (Some(d), Some(c)) = (&self.mouse.wafer_drag, cur) {
                    if let Some(nt) = d.update(&self.rig, c, aspect) { self.session.drag_to(nt); }
                }
            }
        } else if self.mouse.wafer_drag.is_some() && !primary { self.session.end_drag(); self.mouse.wafer_drag = None; }

        // Place / Select / Inspect act on a click (Place also follows a drag so the preview tracks the cursor).
        let clicked = resp.clicked_by(egui::PointerButton::Primary);
        match self.tool {
            Tool::Place if (clicked || primary) && !self.locked() => match self.hover {
                Some(h) => match self.session.propose_site([h.local[0], h.local[1]]) {
                    Ok(()) => self.status = format!("Site preview at ({:.0}, {:.0}) nm. Apply it, or Esc to cancel.", h.local[0], h.local[1]),
                    Err(e) => self.status = e.to_string(),
                },
                None => if clicked { self.status = "That point is not on the wafer.".into(); },
            },
            Tool::Select if clicked => { if let Some(n) = cur { self.select_at_ndc(n); } }
            Tool::Inspect if clicked => { if let Some(n) = cur { self.inspect_at(n); } }
            _ => {}
        }
    }

    /// Select the SVO leaf under an NDC position (highlighted in the viewport, described in the inspector).
    pub fn select_at_ndc(&mut self, ndc: [f64; 2]) -> bool {
        self.selection = self.pick_leaf(ndc).map(|(b, m, _)| (b, m));
        self.status = match &self.selection {
            Some((b, m)) => format!("Selected {:?} leaf: {:.1} × {:.1} × {:.1} nm", m, b.max.0 - b.min.0, b.max.1 - b.min.1, b.max.2 - b.min.2),
            None => "No geometry under the cursor.".into(),
        };
        self.selection.is_some()
    }

    fn pick_leaf(&self, ndc: [f64; 2]) -> Option<(Bounds, Material, Vec3)> {
        let ray = self.rig.ray(ndc, self.viewport_aspect)?;
        let (o, d) = self.transform().ray_to_local(ray.0, ray.1);
        let clip = self.section_clip_local();
        let mut best: Option<(f64, Bounds, Material)> = None;
        self.display.wafer.visit_leaves(|b, m| {
            if m == Material::Air { return; }
            if let Some(l) = Palette_layer(m) { if !self.palette.is_visible(l) { return; } }
            let (mut t0, mut t1) = (f64::NEG_INFINITY, f64::INFINITY);
            for a in 0..3 {
                let (lo, hi) = ([b.min.0, b.min.1, b.min.2][a], [b.max.0, b.max.1, b.max.2][a]);
                if d[a].abs() < 1e-12 { if o[a] < lo || o[a] > hi { return; } continue; }
                let (a0, a1) = ((lo - o[a]) / d[a], (hi - o[a]) / d[a]);
                t0 = t0.max(a0.min(a1)); t1 = t1.min(a0.max(a1));
            }
            if let Some((n, dd)) = clip {
                let (num, den) = (dot(n, o) + dd, dot(n, d));
                if den.abs() < 1e-12 { if num < 0.0 { return; } } else { let tc = -num / den; if den > 0.0 { t0 = t0.max(tc); } else { t1 = t1.min(tc); } }
            }
            if t0 <= t1 && t1 > 0.0 {
                let tt = t0.max(0.0);
                if best.as_ref().map_or(true, |x| tt < x.0) { best = Some((tt, b, m)); }
            }
        });
        best.map(|(t, b, m)| (b, m, [o[0] + d[0] * t, o[1] + d[1] * t, o[2] + d[2] * t]))
    }

    fn inspect_at(&mut self, ndc: [f64; 2]) {
        let w = self.display.wafer.clone();
        let hit = self.pick_leaf(ndc).map(|(_, _, p)| p).or_else(|| self.hover.map(|h| h.local));
        let Some(p) = hit else { self.status = "Nothing to inspect there.".into(); return };
        let in_tile = w.tile_contains_xy(p[0], p[1]);
        let runs: Vec<Run> = if in_tile { w.column_runs(p[0], p[1]) } else { vec![Run { z0: 0.0, z1: w.surface_z_nm, material: Material::Silicon }, Run { z0: w.surface_z_nm, z1: w.surface_z_nm * 2.0, material: Material::Air }] };
        let bucket = runs.iter().find_map(|r| if let Material::Dopant(b) = r.material { Some(b) } else { None });
        self.status = format!("Inspected ({:.0}, {:.0}) nm", p[0], p[1]);
        self.inspect = Some(InspectReadout { x: p[0], y: p[1], runs, surface_z: w.surface_z_nm, in_tile, bucket });
    }

    // ----- overlays ---------------------------------------------------------------------------

    fn draw_overlays(&self, painter: &egui::Painter, rect: Rect) {
        let pal = &self.palette;
        let sel_c = col(pal.layer(LegendLayer::SelectionPreview));
        let z = self.surface_z();
        let w = self.session.state.recipe.window;
        let quad = |c: [f64; 2], s: [f64; 2], zz: f64| -> Option<Vec<Pos2>> {
            let (hx, hy) = (s[0] * 0.5, s[1] * 0.5);
            [[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy]].iter().map(|d| self.project_local(rect, [c[0] + d[0], c[1] + d[1], zz])).collect()
        };
        let draw_loop = |pts: &[Pos2], stroke: Stroke, dashed: bool| {
            for i in 0..pts.len() {
                let (a, b) = (pts[i], pts[(i + 1) % pts.len()]);
                if dashed { painter.extend(egui::Shape::dashed_line(&[a, b], stroke, 6.0, 4.0)); } else { painter.line_segment([a, b], stroke); }
            }
        };

        // Tile box + committed process window + pending site preview.
        if self.show_tile_box {
            let edge = self.session.state.spec.tile_edge_nm();
            let b = &self.display.wafer.bounds;
            let c = [(b.min.0 + b.max.0) * 0.5, (b.min.1 + b.max.1) * 0.5];
            if let Some(q) = quad(c, [edge, edge], z) { draw_loop(&q, Stroke::new(1.0_f32, sel_c.gamma_multiply(0.35)), true); }
        }
        if let Some(q) = quad(w.center, w.size, z) {
            draw_loop(&q, Stroke::new(1.5_f32, sel_c), false);
            painter.text(q[0] + Vec2::new(3.0, -3.0), Align2::LEFT_BOTTOM, "process window", FontId::proportional(10.0), sel_c);
        }
        if let Some(s) = self.session.pending_site {
            if let Some(q) = quad(s, w.size, z) {
                draw_loop(&q, Stroke::new(2.0_f32, col(pal.layer(LegendLayer::WarningError))), true);
                painter.text(q[0] + Vec2::new(3.0, -3.0), Align2::LEFT_BOTTOM, "pending site (uncommitted)", FontId::proportional(10.0), col(pal.layer(LegendLayer::WarningError)));
            }
        }

        // Stage overlays (exposure map, metrology sites).
        if !self.rig.mode.is_section() {
            match &self.display.overlay {
                Overlay::Exposure(e) => {
                    let open_c = col(pal.layer(LegendLayer::MaskOpenings));
                    let etched_c = col(pal.layer(LegendLayer::EtchedRegion));
                    for iy in 0..e.ny { for ix in 0..e.nx {
                        let k = iy * e.nx + ix;
                        let v = e.intensity[k];
                        let removed = e.removed.get(k).copied().unwrap_or(false);
                        if v < 0.04 && !removed { continue; }
                        let (x0, y0) = (e.x0 + ix as f64 * e.pitch, e.y0 + iy as f64 * e.pitch);
                        let pts: Option<Vec<Pos2>> = [[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]].iter().map(|d| self.project_local(rect, [x0 + d[0] * e.pitch, y0 + d[1] * e.pitch, z])).collect();
                        if let Some(pts) = pts {
                            let c = if removed { etched_c.gamma_multiply(0.35) } else { open_c.gamma_multiply((v * 0.45).clamp(0.0, 0.45)) };
                            painter.add(egui::Shape::convex_polygon(pts, c, Stroke::NONE));
                        }
                    } }
                }
                Overlay::Sites(sites) => {
                    for s in sites {
                        if let Some(p) = self.project_local(rect, [s.x, s.y, z]) {
                            let c = if s.out_of_spec { col(pal.layer(LegendLayer::WarningError)) } else { col(pal.layer(LegendLayer::Metrology)) };
                            painter.circle_stroke(p, 6.0, Stroke::new(2.0_f32, c));
                            painter.circle_filled(p, 2.0, c);
                            if let Some(v) = s.value { painter.text(p + Vec2::new(8.0, -6.0), Align2::LEFT_BOTTOM, format!("{v:.0}"), FontId::monospace(10.0), c); }
                        }
                    }
                }
                _ => {}
            }
        }

        // Inspect marker.
        if let Some(i) = &self.inspect {
            if let Some(p) = self.project_local(rect, [i.x, i.y, z]) { painter.circle_stroke(p, 7.0, Stroke::new(2.0_f32, sel_c)); painter.line_segment([p - Vec2::new(12.0, 0.0), p + Vec2::new(12.0, 0.0)], Stroke::new(1.0_f32, sel_c)); painter.line_segment([p - Vec2::new(0.0, 12.0), p + Vec2::new(0.0, 12.0)], Stroke::new(1.0_f32, sel_c)); }
        }

        // Selected leaf: an SVO leaf can be a few pixels wide, so always draw a ring with a minimum on-screen size.
        if let Some((b, m)) = &self.selection {
            let (c, e) = ([(b.min.0 + b.max.0) * 0.5, (b.min.1 + b.max.1) * 0.5, (b.min.2 + b.max.2) * 0.5], [b.max.0 - b.min.0, b.max.1 - b.min.1, b.max.2 - b.min.2]);
            if let Some(p) = self.project_local(rect, c) {
                let gold = Color32::from_rgb(255, 209, 46);
                let span = self.project_local(rect, [c[0] + e[0] * 0.5, c[1] + e[1] * 0.5, c[2] + e[2] * 0.5]).map(|q| (q - p).length()).unwrap_or(0.0);
                let r = (span * 1.4).max(11.0);
                painter.circle_stroke(p, r, Stroke::new(2.0_f32, gold));
                for d in [Vec2::X, -Vec2::X, Vec2::Y, -Vec2::Y] { painter.line_segment([p + d * (r + 2.0), p + d * (r + 8.0)], Stroke::new(2.0_f32, gold)); }
                painter.text(p + Vec2::new(r + 10.0, -r), Align2::LEFT_BOTTOM, format!("selected: {m:?}"), FontId::proportional(11.0), gold);
            }
        }

        // Wafer origin + local axes.
        let t = self.transform();
        if let Some(o) = self.project_local(rect, [0.0, 0.0, z]) {
            if rect.contains(o) {
                painter.circle_filled(o, 3.5, Color32::WHITE);
                for (d, c, l) in [([1.0, 0.0, 0.0], Color32::from_rgb(255, 100, 100), "x"), ([0.0, 1.0, 0.0], Color32::from_rgb(110, 225, 110), "y")] {
                    let tip_l = [d[0] * self.rig.distance * 0.06, d[1] * self.rig.distance * 0.06, z];
                    if let Some(tp) = self.project_local(rect, tip_l) { painter.line_segment([o, tp], Stroke::new(1.5_f32, c)); painter.text(tp, Align2::CENTER_BOTTOM, l, FontId::monospace(10.0), c); }
                }
                painter.text(o + Vec2::new(6.0, 6.0), Align2::LEFT_TOP, "wafer origin", FontId::proportional(10.0), Color32::from_gray(190));
            }
        }

        self.draw_hud(painter, rect, &t);
        self.draw_gizmo(painter, rect);
        self.draw_scale_bar(painter, rect);
        self.draw_legend_chips(painter, rect);
    }

    fn draw_hud(&self, painter: &egui::Painter, rect: Rect, t: &sim_core::transform::WaferTransform) {
        let fg = if self.dark { Color32::from_gray(225) } else { Color32::from_gray(30) };
        let f = FontId::proportional(12.0);
        let mut y = rect.min.y + 8.0;
        let mut line = |s: String, c: Color32, painter: &egui::Painter| { painter.text(Pos2::new(rect.min.x + 10.0, y), Align2::LEFT_TOP, s, f.clone(), c); y += 16.0; };
        let w = &self.display.wafer;
        line(format!("View: {}   Tool: {}", self.rig.mode.label(), self.tool.label()), fg, painter);
        line(format!("Wafer Ø{:.0} mm · origin ({:.3}, {:.3}, {:.3}) mm · yaw {:.1}°", w.domain.radius_nm * 2.0 / 1.0e6, t.position[0] / 1.0e6, t.position[1] / 1.0e6, t.position[2] / 1.0e6, t.yaw_rad.to_degrees()), fg, painter);
        let (badge, c) = match self.display.kind {
            DisplayKind::Committed => ("● COMMITTED geometry", Color32::from_rgb(90, 200, 120)),
            DisplayKind::Preview => ("◌ PREVIEW — not committed", Color32::from_rgb(255, 190, 60)),
            DisplayKind::AnimationFrame => ("▶ ANIMATION FRAME — replay, committed state unchanged", Color32::from_rgb(110, 170, 255)),
        };
        line(badge.into(), c, painter);
        if let (Some(ph), Some(st)) = (self.display.phase, self.display.stage) {
            if self.display.kind != DisplayKind::Committed || self.playback.is_playing() {
                line(format!("{} · {} ({:.0}%)", st.label(), ph.name, self.display.progress * 100.0), fg, painter);
                let note = match ph.kind {
                    StepKind::Geometry => "this phase changes SVO geometry".to_string(),
                    StepKind::MaterialTag => "this phase re-labels material (dopant tags); no growth is implied".to_string(),
                    StepKind::Analytic | StepKind::Recorded => format!("{} phase: geometry is unchanged", ph.kind.label()),
                };
                line(note, Color32::from_gray(150), painter);
            }
        }
        match &self.display.overlay {
            Overlay::Thickness { target_nm, mean_nm, min_nm, max_nm } => line(format!("film: mean {mean_nm:.1} nm (min {min_nm:.1}, max {max_nm:.1}) · target {target_nm:.0} nm"), fg, painter),
            Overlay::Surface { clean_fraction, particles_cm2 } => line(format!("surface clean {:.0}% · particles {:.2e} cm⁻²", clean_fraction * 100.0, particles_cm2), fg, painter),
            _ => {}
        }
        if self.locked() { line("Interaction lock: process animation playing (pause to move the wafer or place)".into(), Color32::from_rgb(255, 190, 60), painter); }
        // Coordinate readout.
        let txt = match self.hover {
            Some(h) => format!("world ({:.3}, {:.3}, {:.3}) mm   wafer-local ({:.1}, {:.1}) nm", h.world[0] / 1.0e6, h.world[1] / 1.0e6, h.world[2] / 1.0e6, h.local[0], h.local[1]),
            None => "cursor not over the work plane".into(),
        };
        painter.text(Pos2::new(rect.center().x, rect.max.y - 8.0), Align2::CENTER_BOTTOM, txt, FontId::monospace(11.0), fg);
    }

    fn draw_gizmo(&self, painter: &egui::Painter, rect: Rect) {
        let origin = Pos2::new(rect.min.x + 62.0, rect.max.y - 62.0);
        let (f, r, u) = self.rig.basis();
        painter.circle_filled(origin, 40.0, Color32::from_black_alpha(if self.dark { 70 } else { 20 }));
        for (e, c, l) in [([1.0, 0.0, 0.0], Color32::from_rgb(255, 100, 100), "X"), ([0.0, 1.0, 0.0], Color32::from_rgb(110, 225, 110), "Y"), ([0.0, 0.0, 1.0], Color32::from_rgb(100, 160, 255), "Z")] {
            let (sx, sy, depth) = (dot(e, r), dot(e, u), dot(e, f));
            let tip = origin + Vec2::new(sx as f32, -sy as f32) * 32.0;
            let c = c.gamma_multiply(if depth > 0.3 { 0.45 } else { 1.0 });
            painter.line_segment([origin, tip], Stroke::new(2.0_f32, c));
            painter.text(tip + Vec2::new(sx as f32, -sy as f32) * 7.0, Align2::CENTER_CENTER, l, FontId::proportional(12.0), c);
        }
    }

    fn draw_scale_bar(&self, painter: &egui::Painter, rect: Rect) {
        let npp = self.rig.nm_per_pixel(rect.height() as f64);
        let (len, label) = nice_length(npp * 130.0);
        let px = (len / npp) as f32;
        let fg = if self.dark { Color32::from_gray(230) } else { Color32::from_gray(20) };
        let (x1, y) = (rect.max.x - 20.0, rect.max.y - 34.0);
        painter.line_segment([Pos2::new(x1 - px, y), Pos2::new(x1, y)], Stroke::new(3.0_f32, fg));
        painter.line_segment([Pos2::new(x1 - px, y - 5.0), Pos2::new(x1 - px, y + 5.0)], Stroke::new(1.5_f32, fg));
        painter.line_segment([Pos2::new(x1, y - 5.0), Pos2::new(x1, y + 5.0)], Stroke::new(1.5_f32, fg));
        painter.text(Pos2::new(x1 - px * 0.5, y - 8.0), Align2::CENTER_BOTTOM, label, FontId::proportional(12.0), fg);
    }

    fn draw_legend_chips(&self, painter: &egui::Painter, rect: Rect) {
        let active: Vec<LegendLayer> = match self.display.stage.or_else(|| self.session.state.recipe.index_of(self.selected_stage).map(|i| self.session.state.recipe.stages[i].params.stage())) {
            Some(s) => layers_for_stage(s).to_vec(),
            None => vec![LegendLayer::Substrate],
        };
        let mut y = rect.min.y + 8.0;
        let fg = if self.dark { Color32::from_gray(225) } else { Color32::from_gray(30) };
        for l in active {
            let on = self.palette.is_visible(l);
            let c = col(self.palette.layer(l));
            let x = rect.max.x - 150.0;
            painter.rect_filled(Rect::from_min_size(Pos2::new(x, y + 2.0), Vec2::splat(11.0)), 2.0, if on { c } else { c.gamma_multiply(0.2) });
            painter.text(Pos2::new(x + 17.0, y), Align2::LEFT_TOP, l.label(), FontId::proportional(11.5), if on { fg } else { fg.gamma_multiply(0.4) });
            y += 15.0;
        }
    }
}

fn rel_dir(d: Vec3) -> [f32; 3] { [d[0] as f32, d[1] as f32, d[2] as f32] }

#[allow(non_snake_case)]
fn Palette_layer(m: Material) -> Option<LegendLayer> { sim_core::legend::Palette::material_layer(m) }
