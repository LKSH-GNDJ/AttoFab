//! AttoFab the earlier build application shell (PRD §5 layout).
//!
//! Layout: top toolbar · left recipe · centre 3D viewport (+ overlays) · right
//! process inspector · bottom timeline / playback / log. All state transitions go
//! through `sim_core::session::Session`; this file is presentation and input
//! routing only.

use eframe::egui::{self, Color32, RichText, Vec2};
use sim_core::legend::{Palette, Theme};
use sim_core::physics::Stage;
use sim_core::pipeline::{Overlay, SubStep, FRAMES_PER_STAGE};
use sim_core::playback::{PlayState, Playback};
use sim_core::render_wgpu::{flatten_for_gpu_relative, GpuSvoNode};
use sim_core::session::{FrameCache, Session, StageStatus, Timeline};
use sim_core::transform::WaferTransform;
use sim_core::wafer::{Bounds, Material, Run, WaferState};
use std::sync::Arc;

use crate::svo_view::SvoViewportCallback;
use crate::viewport_math::{CameraRig, PanDrag, ViewMode, WaferDrag};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tool { Navigate, Select, MoveWafer, Place, Inspect }

impl Tool {
    pub const ALL: [Tool; 5] = [Tool::Navigate, Tool::Select, Tool::MoveWafer, Tool::Place, Tool::Inspect];
    pub fn label(self) -> &'static str {
        match self { Tool::Navigate => "Navigate", Tool::Select => "Select", Tool::MoveWafer => "Move Wafer", Tool::Place => "Place / Edit", Tool::Inspect => "Inspect" }
    }
    pub fn hint(self) -> &'static str {
        match self {
            Tool::Navigate => "Left-drag: orbit (3D) or pan (side/top views). Middle-drag: pan. Wheel: zoom toward cursor.   [N]",
            Tool::Select => "Left-click selects the SVO leaf under the cursor.   [S]",
            Tool::MoveWafer => "Left-drag moves the wafer in the environment; the point you grab stays under the cursor. The SVO is never rewritten.   [M]",
            Tool::Place => "Left-click sets the process site and previews its window. Dragging never moves the wafer.   [P]",
            Tool::Inspect => "Left-click shows the material stack, heights and measurements at that point.   [I]",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum DisplayKind { Committed, Preview, AnimationFrame }

pub(crate) struct Display {
    pub wafer: Arc<WaferState>,
    pub overlay: Overlay,
    pub nodes: Arc<Vec<GpuSvoNode>>,
    pub revision: u64,
    pub ptr: usize,
    pub kind: DisplayKind,
    pub phase: Option<SubStep>,
    pub stage: Option<Stage>,
    pub progress: f64,
}

#[derive(Default)]
pub(crate) struct MouseState {
    pub wafer_drag: Option<WaferDrag>,
    pub pan: Option<PanDrag>,
}

pub(crate) struct InspectReadout {
    pub x: f64,
    pub y: f64,
    pub runs: Vec<Run>,
    pub surface_z: f64,
    pub in_tile: bool,
    pub bucket: Option<u8>,
}

pub struct AttoFabApp {
    pub(crate) session: Session,
    pub(crate) playback: Playback,
    pub(crate) timeline: Timeline,
    pub(crate) timeline_sig: u64,
    pub(crate) frames: FrameCache,
    pub(crate) palette: Palette,
    pub(crate) rig: CameraRig,
    pub(crate) tool: Tool,
    pub(crate) selected_stage: u32,
    pub(crate) slice_offset_nm: f64,
    pub(crate) mouse: MouseState,
    pub(crate) hover: Option<crate::viewport_math::PlacePick>,
    pub(crate) inspect: Option<InspectReadout>,
    pub(crate) selection: Option<(Bounds, Material)>,
    pub(crate) display: Display,
    pub(crate) revision_counter: u64,
    pub(crate) gpu_ok: bool,
    pub(crate) logo: Option<egui::TextureHandle>,
    pub(crate) dark: bool,
    pub(crate) status: String,
    pub(crate) new_stage_kind: Stage,
    pub(crate) autoplay_pending: bool,
    pub(crate) viewport_aspect: f64,
    pub(crate) need_focus: bool,
    pub(crate) show_log: bool,
    pub(crate) show_tile_box: bool,
    pub(crate) theme_applied: Option<bool>,
    pub(crate) viewport_rect: egui::Rect,
    pub(crate) vv_report: Option<sim_core::vv::VvReport>,
    pub(crate) vv_rx: Option<std::sync::mpsc::Receiver<sim_core::vv::VvReport>>,
    pub(crate) gds_path: String,
    pub(crate) gds: Option<GdsState>,
    pub(crate) gds_msg: Option<(bool, String)>,
    pub(crate) gds_run_rx: Option<std::sync::mpsc::Receiver<Result<sim_core::gds_sweep::WholeFileRun, String>>>,
    pub(crate) gds_run_progress: Option<std::sync::Arc<sim_core::gds_sweep::RunProgress>>,
    pub(crate) gds_sweep_rx: Option<std::sync::mpsc::Receiver<Result<sim_core::gds_sweep::SweepPlan, String>>>,
    pub(crate) gds_sel: usize,
    pub(crate) gds_crop: [f64; 2],
    pub(crate) gds_view: crate::gds_view::GdsView,
    pub(crate) gds_view_open: bool,
    pub(crate) gds_view_fit: bool,
    pub(crate) gds_show_all: bool,
    /// Every layer goes into one mask, as a combined etch sees them (W-14).
    pub(crate) gds_all_layers: bool,
    /// The planned sweep, once the user asks for one.
    pub(crate) gds_sweep: Option<sim_core::gds_sweep::SweepPlan>,
    pub(crate) gds_canvas_rect: egui::Rect,
    pub(crate) scale_edit: (f64, u32),
    pub(crate) scale_synced: (f64, u32),
}

/// A loaded, flattened GDSII layout ready to pick a layer from.
pub(crate) struct GdsState {
    pub file: String,
    pub layout: Arc<sim_core::io::gds::GdsLayout>,
    pub top: String,
    pub tops: Vec<String>,
    pub flat: Arc<Vec<sim_core::io::gds::FlatPoly>>,
    pub layers: Vec<sim_core::io::gds::LayerSummary>,
    /// Per-polygon bounding boxes and per-layer index lists, for fast drawing and local statistics.
    pub bboxes: Vec<[f64; 4]>,
    pub by_layer: std::collections::BTreeMap<(i32, i32), Vec<u32>>,
}

impl GdsState {
    pub(crate) fn build(file: String, layout: Arc<sim_core::io::gds::GdsLayout>, top: String, tops: Vec<String>, flat: Vec<sim_core::io::gds::FlatPoly>) -> Self {
        let layers = sim_core::io::gds::summarize(&flat);
        let bboxes: Vec<[f64; 4]> = flat.iter().map(|p| p.pts.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)])).collect();
        let mut by_layer: std::collections::BTreeMap<(i32, i32), Vec<u32>> = Default::default();
        for (i, p) in flat.iter().enumerate() { by_layer.entry((p.layer, p.datatype)).or_default().push(i as u32); }
        Self { file, layout, top, tops, flat: Arc::new(flat), layers, bboxes, by_layer }
    }
}

pub(crate) fn col(c: [u8; 3]) -> Color32 { Color32::from_rgb(c[0], c[1], c[2]) }

pub(crate) fn status_color(s: &StageStatus) -> Color32 {
    match s {
        StageStatus::Committed => Color32::from_rgb(90, 200, 120), StageStatus::Stale => Color32::from_rgb(240, 180, 60),
        StageStatus::NotRun => Color32::from_rgb(140, 150, 165), StageStatus::Failed(_) => Color32::from_rgb(255, 80, 80), StageStatus::Disabled => Color32::from_rgb(90, 95, 105),
    }
}

pub(crate) fn status_word(s: &StageStatus) -> &'static str {
    match s { StageStatus::Committed => "committed", StageStatus::Stale => "stale", StageStatus::NotRun => "not run", StageStatus::Failed(_) => "failed", StageStatus::Disabled => "disabled" }
}

fn load_logo(ctx: &egui::Context, dark: bool) -> Option<egui::TextureHandle> {
    let bytes = include_bytes!("../assets/attofab_icon_1024.png");
    let img = image::load_from_memory(bytes).ok()?.to_rgba8();
    let size = [img.width() as usize, img.height() as usize];
    let mut px = img.into_raw();
    for p in px.chunks_exact_mut(4) {
        // The shipped icon is black artwork; for the dark theme keep its alpha and paint it white.
        if dark { p[0] = 245; p[1] = 248; p[2] = 252; }
    }
    Some(ctx.load_texture(if dark { "attofab-icon-dark" } else { "attofab-icon-light" }, egui::ColorImage::from_rgba_unmultiplied(size, &px), egui::TextureOptions::LINEAR))
}

impl AttoFabApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        let gpu_ok = SvoViewportCallback::install(cc);
        Self::build(&cc.egui_ctx, gpu_ok)
    }

    /// Construct without a GPU (tests, tooling): every panel, tool and overlay still runs.
    pub fn new_headless(ctx: &egui::Context) -> Self { Self::build(ctx, false) }

    fn build(ctx: &egui::Context, gpu_ok: bool) -> Self {
        let session = Session::with_default_recipe();
        let selected = session.state.recipe.stages[0].id;
        let wafer = session.committed_wafer();
        let origin = wafer.tile_origin();
        let nodes = Arc::new(flatten_for_gpu_relative(&wafer, origin));
        let ptr = Arc::as_ptr(&wafer) as usize;
        let mut app = Self {
            session, playback: Playback::default(), timeline: Timeline::default(), timeline_sig: 0, frames: FrameCache::default(), palette: Palette::new(Theme::Dark),
            rig: CameraRig::new(ViewMode::Perspective), tool: Tool::Navigate, selected_stage: selected, slice_offset_nm: 0.0, mouse: MouseState::default(), hover: None,
            inspect: None, selection: None,
            display: Display { wafer, overlay: Overlay::None, nodes, revision: 1, ptr, kind: DisplayKind::Committed, phase: None, stage: None, progress: 1.0 },
            revision_counter: 1, gpu_ok, logo: load_logo(ctx, true), dark: true, status: "Ready. Press Run Process, or select a stage and Preview it.".into(),
            new_stage_kind: Stage::Deposition, autoplay_pending: false, viewport_aspect: 1.6, need_focus: true, show_log: true, show_tile_box: true, theme_applied: None, viewport_rect: egui::Rect::NOTHING, vv_report: None, vv_rx: None, gds_path: String::new(), gds: None, gds_msg: None, gds_run_rx: None, gds_run_progress: None, gds_sweep_rx: None, gds_sel: 0, gds_crop: [0.0, 0.0], gds_view: Default::default(), gds_view_open: false, gds_view_fit: false, gds_show_all: true, gds_all_layers: false, gds_sweep: None, gds_canvas_rect: egui::Rect::NOTHING, scale_edit: (16.0, 7), scale_synced: (16.0, 7),
        };
        app.rig.mode = ViewMode::Perspective;
        app
    }

    // ----- helpers ---------------------------------------------------------------------

    pub(crate) fn transform(&self) -> WaferTransform { self.session.state.transform }
    pub(crate) fn wafer_radius(&self) -> f64 { self.session.state.initial.domain.radius_nm }
    pub(crate) fn surface_z(&self) -> f64 { self.session.state.initial.surface_z_nm }
    pub(crate) fn site(&self) -> [f64; 2] { self.session.state.recipe.window.center }
    /// Test/tooling helpers.
    pub fn stop_and_show_final(&mut self) { self.playback.state = PlayState::Stopped; self.playback.seek(self.timeline.segments.len() as f64); }
    pub fn set_layer_visible(&mut self, l: sim_core::legend::LegendLayer, v: bool) { self.palette.visible[l.index()] = v; }
    /// Show playback position `segment + fraction` (paused) — for contact sheets and tests.
    pub fn show_position(&mut self, segment: usize, fraction: f64) { self.playback.state = PlayState::Stopped; self.playback.seek(segment as f64 + fraction); }
    pub fn site_window(&self) -> sim_core::geom::ProcessWindow { self.session.state.recipe.window }
    pub fn surface_z_nm(&self) -> f64 { self.surface_z() }
    pub fn set_aspect(&mut self, a: f64) { self.viewport_aspect = a; }
    /// Frame a view for offscreen rendering: `site` = the process tile, otherwise the whole wafer.
    pub fn frame_view(&mut self, mode: ViewMode, site: bool) {
        self.rig.mode = mode;
        if site { self.focus_site(); } else { self.fit_wafer(); }
        if mode == ViewMode::Xz || mode == ViewMode::Yz { self.zoom_for_section(); }
    }
    /// Zoom a section view onto the tile (the cut face is only a few µm wide).
    pub fn zoom_for_section(&mut self) { self.rig.distance = CameraRig::fit_distance(self.rig.mode, self.session.state.spec.tile_edge_nm() * 0.55, self.viewport_aspect); }
    pub(crate) fn locked(&self) -> bool { self.playback.is_playing() }

    pub(crate) fn set_tool(&mut self, t: Tool) {
        if self.mouse.wafer_drag.is_some() { self.session.end_drag(); self.mouse.wafer_drag = None; }
        self.tool = t;
        self.status = t.hint().to_string();
    }

    pub(crate) fn site_world(&self) -> [f64; 3] { self.transform().local_to_world([self.site()[0], self.site()[1], self.surface_z()]) }

    pub(crate) fn focus_site(&mut self) {
        let edge = self.session.state.spec.tile_edge_nm();
        self.rig.target = self.site_world();
        self.rig.distance = CameraRig::fit_distance(ViewMode::Top, edge * 0.55, self.viewport_aspect);
    }

    pub(crate) fn fit_wafer(&mut self) {
        let t = self.transform();
        self.rig.target = t.local_to_world([0.0, 0.0, self.surface_z()]);
        self.rig.distance = CameraRig::fit_distance(self.rig.mode, self.wafer_radius(), self.viewport_aspect);
    }

    pub(crate) fn reset_view(&mut self) {
        self.rig = CameraRig::new(self.rig.mode);
        self.focus_site();
    }

    // ----- timeline / display ---------------------------------------------------------

    fn sync_timeline(&mut self) {
        let tl = self.session.timeline();
        let sig = tl.signature() ^ (tl.segments.len() as u64).wrapping_mul(0x9E37_79B9);
        if sig != self.timeline_sig {
            self.timeline_sig = sig;
            self.playback.set_segments(tl.segments.len());
            self.playback.seek(tl.segments.len() as f64);
            self.playback.state = PlayState::Stopped;
            if self.autoplay_pending && !tl.segments.is_empty() { self.playback.reset(); self.playback.play(); }
            self.autoplay_pending = false;
        }
        self.timeline = tl;
    }

    fn refresh_display(&mut self) {
        let (wafer, overlay, kind, phase, stage, progress) = if self.timeline.segments.is_empty() {
            (self.session.committed_wafer(), Overlay::None, DisplayKind::Committed, None, None, 1.0)
        } else {
            let (si, k) = self.playback.at_keyframe();
            let seg = self.timeline.segments[si.min(self.timeline.segments.len() - 1)].clone();
            match self.frames.get(&seg, k) {
                Ok(f) => {
                    let at_end = si + 1 == self.timeline.segments.len() && k >= FRAMES_PER_STAGE;
                    let kind = if self.session.preview.is_some() { DisplayKind::Preview } else if self.playback.is_playing() || !at_end { DisplayKind::AnimationFrame } else { DisplayKind::Committed };
                    (f.wafer.clone(), f.overlay.clone(), kind, Some(f.phase), Some(f.stage), f.progress)
                }
                Err(e) => { self.status = format!("Frame replay failed: {e}"); (self.session.committed_wafer(), Overlay::None, DisplayKind::Committed, None, None, 1.0) }
            }
        };
        let ptr = Arc::as_ptr(&wafer) as usize;
        if ptr != self.display.ptr {
            self.revision_counter += 1;
            self.display.nodes = Arc::new(flatten_for_gpu_relative(&wafer, wafer.tile_origin()));
            self.display.revision = self.revision_counter;
            self.display.ptr = ptr;
            self.selection = None;
        }
        self.display.wafer = wafer;
        self.display.overlay = overlay;
        self.display.kind = kind;
        self.display.phase = phase;
        self.display.stage = stage;
        self.display.progress = progress;
    }

    // ----- actions ---------------------------------------------------------------------

    pub(crate) fn do_preview(&mut self) {
        match self.session.preview_stage(self.selected_stage) {
            Ok(()) => { self.autoplay_pending = true; self.status = "Preview computed. Nothing is committed until you press Commit.".into(); }
            Err(e) => self.status = format!("Cannot preview: {e}"),
        }
    }
    pub(crate) fn do_commit(&mut self) {
        match self.session.commit_preview() { Ok(_) => self.status = "Committed.".into(), Err(e) => self.status = format!("Cannot commit: {e}") }
    }
    pub(crate) fn do_cancel(&mut self) { self.session.cancel_preview(); self.status = "Cancelled. Committed state is unchanged.".into(); }
    pub fn do_run(&mut self) {
        let r = self.session.run_recipe();
        self.status = match (&r.failed, r.ran.len()) {
            (Some((_, m)), n) => format!("Run stopped after {n} stage(s): {m}"),
            (None, 0) => "Everything is already committed. Edit a parameter to make downstream stages stale.".into(),
            (None, n) => { self.autoplay_pending = true; format!("Ran {n} stage(s); {} already up to date.", r.skipped_committed.len()) }
        };
    }

    fn shortcuts(&mut self, ctx: &egui::Context) {
        if ctx.wants_keyboard_input() { return; }
        let (undo, redo, space, esc, keys) = ctx.input(|i| {
            let cmd = i.modifiers.command;
            (cmd && !i.modifiers.shift && i.key_pressed(egui::Key::Z), (cmd && i.modifiers.shift && i.key_pressed(egui::Key::Z)) || (cmd && i.key_pressed(egui::Key::Y)),
             i.key_pressed(egui::Key::Space), i.key_pressed(egui::Key::Escape),
             [egui::Key::N, egui::Key::S, egui::Key::M, egui::Key::P, egui::Key::I, egui::Key::F].map(|k| !i.modifiers.command && i.key_pressed(k)))
        });
        if undo { self.session.undo(); }
        if redo { self.session.redo(); }
        if space { self.playback.toggle(); }
        if esc { self.do_cancel(); }
        for (k, t) in keys.iter().zip(Tool::ALL) { if *k { self.set_tool(t); } }
        if keys[5] { self.fit_wafer(); }
    }

    fn apply_theme(&mut self, ctx: &egui::Context) {
        if self.theme_applied == Some(self.dark) { return; }
        ctx.set_visuals(if self.dark { egui::Visuals::dark() } else { egui::Visuals::light() });
        let vis = self.palette.visible;
        self.palette = Palette::new(if self.dark { Theme::Dark } else { Theme::Light });
        self.palette.visible = vis;
        self.logo = load_logo(ctx, self.dark);
        self.theme_applied = Some(self.dark);
    }

    // ----- panels -----------------------------------------------------------------------

    fn toolbar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal_wrapped(|ui| {
            if let Some(l) = &self.logo { ui.add(egui::Image::new(l).fit_to_exact_size(Vec2::splat(24.0))); }
            ui.label(RichText::new("AttoFab").strong().size(16.0));
            ui.label(RichText::new("the earlier build").weak());
            ui.separator();
            for t in Tool::ALL {
                if ui.selectable_label(self.tool == t, t.label()).on_hover_text(t.hint()).clicked() { self.set_tool(t); }
            }
            ui.separator();
            if ui.add_enabled(self.session.can_undo() && !self.locked(), egui::Button::new("↶ Undo")).on_hover_text("Ctrl+Z").clicked() { self.session.undo(); }
            if ui.add_enabled(self.session.can_redo() && !self.locked(), egui::Button::new("↷ Redo")).on_hover_text("Ctrl+Y").clicked() { self.session.redo(); }
            ui.separator();
            for m in ViewMode::ALL {
                if ui.selectable_label(self.rig.mode == m, m.label()).clicked() { self.rig.mode = m; }
            }
            if self.rig.mode == ViewMode::Slice {
                ui.add(egui::DragValue::new(&mut self.slice_offset_nm).speed(2.0).clamp_range(-1000.0..=1000.0).suffix(" nm").prefix("slice Δz "));
            }
            ui.separator();
            if ui.button("Reset View").clicked() { self.reset_view(); }
            if ui.button("Center Wafer").on_hover_text("Move the wafer back to the world origin (undoable)").clicked() { self.session.center_wafer(); }
            if ui.button("Fit Wafer").on_hover_text("Frame the whole 300 mm wafer   [F]").clicked() { self.fit_wafer(); }
            if ui.button("Focus Site").on_hover_text("Zoom to the high-resolution process tile").clicked() { self.focus_site(); }
            ui.separator();
            let mut d = self.dark;
            if ui.selectable_label(d, if d { "🌙 Dark" } else { "☀ Light" }).clicked() { d = !d; }
            self.dark = d;
        });
        ui.horizontal_wrapped(|ui| {
            ui.label(RichText::new("Process Tool").strong());
            let sel = self.session.state.recipe.index_of(self.selected_stage).map(|i| self.session.state.recipe.stages[i].params.stage().label()).unwrap_or("—");
            ui.label(format!("stage: {sel}"));
            let has_prev = self.session.preview.is_some();
            if ui.add_enabled(!self.locked(), egui::Button::new("👁 Preview")).on_hover_text("Compute the selected stage on the committed state without changing it").clicked() { self.do_preview(); }
            if ui.add_enabled(has_prev, egui::Button::new("✔ Commit")).clicked() { self.do_commit(); }
            if ui.add_enabled(has_prev || self.session.pending_site.is_some(), egui::Button::new("✖ Cancel").fill(Color32::TRANSPARENT)).on_hover_text("Esc").clicked() { self.do_cancel(); }
            ui.separator();
            if ui.add_enabled(!self.locked(), egui::Button::new("▶ Run Process")).on_hover_text("Run every stage that is not already committed, then replay it").clicked() { self.do_run(); }
            if ui.add_enabled(!self.locked(), egui::Button::new("⟲ Reset Wafer")).on_hover_text("Back to the blank substrate (undoable)").clicked() { self.session.reset(); self.status = "Wafer reset to blank substrate.".into(); }
            if let Some(site) = self.session.pending_site {
                ui.separator();
                ui.label(RichText::new(format!("pending site ({:.0}, {:.0}) nm", site[0], site[1])).color(col(self.palette.layer(sim_core::legend::LegendLayer::SelectionPreview))));
                if ui.button("Apply site").on_hover_text("Re-centres the tile; every committed stage becomes stale").clicked() {
                    match self.session.commit_site() { Ok(()) => { self.need_focus = true; self.status = "Site applied. Downstream results are stale.".into(); } Err(e) => self.status = e.to_string() }
                }
            }
            ui.separator();
            ui.label(RichText::new(&self.status).weak());
        });
    }

    fn recipe_panel(&mut self, ui: &mut egui::Ui) {
        egui::ScrollArea::vertical().id_source("recipe_scroll").auto_shrink([false, false]).show(ui, |ui| {
        ui.heading("Recipe");
        let statuses = self.session.statuses();
        let committed = statuses.iter().filter(|s| s.1 == StageStatus::Committed).count();
        ui.label(RichText::new(format!("{committed}/{} committed", self.session.enabled_count())).weak());
        {
            let stages = self.session.state.recipe.stages.clone();
            for (i, (rs, (_, status))) in stages.iter().zip(statuses.iter()).enumerate() {
                ui.horizontal(|ui| {
                    let mut en = rs.enabled;
                    if ui.checkbox(&mut en, "").on_hover_text("Include in the recipe").changed() { let _ = self.session.set_enabled(rs.id, en); }
                    ui.label(RichText::new("●").color(status_color(status)));
                    let text = format!("{}. {}", i + 1, rs.params.stage().label());
                    let r = ui.selectable_label(self.selected_stage == rs.id, text).on_hover_text(match status { StageStatus::Failed(m) => format!("failed: {m}"), s => status_word(s).to_string() });
                    if r.clicked() { self.selected_stage = rs.id; }
                });
            }
        }
        ui.separator();
        ui.horizontal(|ui| {
            egui::ComboBox::from_id_source("add_stage").selected_text(self.new_stage_kind.label()).show_ui(ui, |ui| {
                for s in Stage::ALL { ui.selectable_value(&mut self.new_stage_kind, s, s.label()); }
            });
            if ui.button("＋ Add after").on_hover_text("Insert after the selected stage; cycles such as deposition → litho → etch → clean are allowed").clicked() {
                if let Ok(id) = self.session.add_stage_after(self.selected_stage, self.new_stage_kind) { self.selected_stage = id; }
            }
        });
        ui.horizontal(|ui| {
            if ui.button("Duplicate").clicked() { if let Ok(id) = self.session.duplicate_stage(self.selected_stage) { self.selected_stage = id; } }
            if ui.button("Delete").clicked() && self.session.state.recipe.stages.len() > 1 {
                let idx = self.session.state.recipe.index_of(self.selected_stage).unwrap_or(0);
                let _ = self.session.remove_stage(self.selected_stage);
                self.selected_stage = self.session.state.recipe.stages[idx.min(self.session.state.recipe.stages.len() - 1)].id;
            }
            if ui.button("▲").clicked() { self.session.move_stage(self.selected_stage, -1); }
            if ui.button("▼").clicked() { self.session.move_stage(self.selected_stage, 1); }
        });
        ui.separator();
        ui.label(RichText::new("Process site").strong());
        let w = self.session.state.recipe.window;
        ui.label(format!("centre ({:.0}, {:.0}) nm   wafer-local", w.center[0], w.center[1]));
        let mut size = w.size;
        let edge = self.session.max_window_nm();
        ui.horizontal(|ui| {
            ui.label("window");
            ui.add(egui::DragValue::new(&mut size[0]).speed(10.0).clamp_range(64.0..=edge).suffix(" nm"));
            ui.label("×");
            ui.add(egui::DragValue::new(&mut size[1]).speed(10.0).clamp_range(64.0..=edge).suffix(" nm"));
        });
        if size != w.size { let _ = self.session.set_window_size(size); }
        ui.label(RichText::new(format!("tile {:.0} nm, cell {:.0} nm", self.session.state.spec.tile_edge_nm(), self.session.state.spec.leaf_nm)).weak());
        ui.separator();
        self.scale_panel(ui);
        ui.separator();
        self.gds_panel(ui);
        });
    }

    fn gds_panel(&mut self, ui: &mut egui::Ui) {
        ui.label(RichText::new("GDSII mask").strong());
        #[cfg(not(target_arch = "wasm32"))]
        ui.horizontal(|ui| {
            let r = ui.add(egui::TextEdit::singleline(&mut self.gds_path).hint_text("C:\\path\\layout.gds  (or drop a file)").desired_width(ui.available_width() - 56.0));
            if ui.button("Load").clicked() || (r.lost_focus() && ui.input(|i| i.key_pressed(egui::Key::Enter))) { self.load_gds_path(); }
        });
        #[cfg(target_arch = "wasm32")]
        ui.label(RichText::new("Drop a .gds file on the window.").weak().small());
        if let Some((ok, msg)) = &self.gds_msg {
            let c = if *ok { Color32::from_rgb(90, 200, 120) } else { col(self.palette.layer(sim_core::legend::LegendLayer::WarningError)) };
            ui.label(RichText::new(msg).color(c).small());
        }
        let mut pick: Option<usize> = None;
        let mut new_top: Option<String> = None;
        let mut apply = false;
        let mut open_viewer = false;
        let mut finest = false;
        let mut nearest = false;
        let mut auto = false;
        let mut sweep = false;
        let mut run_file = false;
        if let Some(g) = &self.gds {
            if g.layout.structures.len() > 1 {
                ui.horizontal(|ui| {
                    ui.label("top cell");
                    egui::ComboBox::from_id_source("gds_top").selected_text(&g.top).show_ui(ui, |ui| { for t in g.layout.structures.keys() { if ui.selectable_label(*t == g.top, t).clicked() { new_top = Some(t.clone()); } } });
                });
            }
            ui.label(RichText::new(format!("{} — layers (click to select)", std::path::Path::new(&g.file).file_name().map(|n| n.to_string_lossy().into_owned()).unwrap_or_else(|| g.file.clone()))).small().weak());
            egui::ScrollArea::vertical().max_height(150.0).id_source("gds_layers").show(ui, |ui| {
                for (i, l) in g.layers.iter().enumerate() {
                    let (w, h) = l.size_nm();
                    let size = if w.max(h) >= 1.0e3 { format!("{:.2} × {:.2} µm", w / 1e3, h / 1e3) } else { format!("{w:.0} × {h:.0} nm") };
                    if ui.selectable_label(self.gds_sel == i, format!("{}/{}   {} poly   {}", l.layer, l.datatype, l.polygons, size)).clicked() { pick = Some(i); }
                }
            });
            ui.horizontal(|ui| {
                ui.label("crop centre");
                ui.add(egui::DragValue::new(&mut self.gds_crop[0]).speed(10.0).suffix(" nm"));
                ui.add(egui::DragValue::new(&mut self.gds_crop[1]).speed(10.0).suffix(" nm"));
            });
            ui.label(RichText::new("Only the area around the crop centre (tile size + margin) is simulated.").weak().small());
            ui.horizontal_wrapped(|ui| {
                if ui.button("Open layout viewer").on_hover_text("See the layout and click where to simulate").clicked() { open_viewer = true; }
                if ui.button("Finest feature").on_hover_text("Move the crop onto the smallest shape of the selected layer").clicked() { finest = true; }
                if ui.button("Nearest feature").on_hover_text("Move the crop onto the closest polygon of this layer").clicked() { nearest = true; }
                if ui.button("Use selected layer as Lithography mask").clicked() { apply = true; }
            });
            if ui.button("Auto-scale & use layer").on_hover_text("Choose cell size, tile and window from the local feature size and your thinnest film, then load the mask").clicked() { auto = true; }
            ui.separator();
            ui.horizontal_wrapped(|ui| {
                ui.checkbox(&mut self.gds_all_layers, "All layers (union mask)").on_hover_text(
                    "A combined etch sees every layer as resist at once, so a feature on any layer protects the same silicon.");
                if ui.button("Plan whole-file sweep").on_hover_text(
                    "The window is capped at 192 cells per axis, so a whole layout is simulated as a grid of tiles, not one run.").clicked() { sweep = true; }
                if ui.button("Run whole file (deduped)").on_hover_text(
                    "Simulates each distinct cell context once and stamps it at every placement.").clicked() { run_file = true; }
            });
        }
        if let Some(t) = new_top { self.set_gds_top(&t); }
        if let Some(i) = pick { self.select_gds_layer(i); }
        if apply { self.use_gds_layer(); }
        if finest { self.jump_to_finest_feature(); let _ = self.rebuild_mask_for_crop(); }
        if nearest { self.jump_to_nearest_feature(); let _ = self.rebuild_mask_for_crop(); }
        if auto { self.auto_scale_from_gds(); }
        if sweep { self.plan_whole_file_sweep(); }
        if run_file { self.run_whole_file(); }
        if open_viewer { self.gds_view_open = true; self.gds_view_fit = true; }
    }

    /// Accept what users actually paste: surrounding quotes ("Copy as path"), `file://` URLs, stray spaces.
    pub fn normalize_path(raw: &str) -> String {
        let t = raw.trim().trim_matches(|c| c == '"' || c == '\'').trim();
        let t = t.strip_prefix("file:///").or_else(|| t.strip_prefix("file://")).unwrap_or(t);
        t.to_string()
    }

    #[cfg(not(target_arch = "wasm32"))]
    pub fn load_gds_path(&mut self) {
        let path = Self::normalize_path(&self.gds_path);
        self.gds_path = path.clone();
        if path.is_empty() { self.gds_msg = Some((false, "Enter the path of a .gds file, or drop one on the window.".into())); return; }
        match std::fs::read(&path) {
            Ok(bytes) => self.load_gds_bytes(&path, &bytes),
            Err(e) => self.gds_msg = Some((false, format!("Cannot read '{path}': {e}"))),
        }
    }

    /// Parse a GDSII stream, flatten its hierarchy and list its layers.
    pub fn load_gds_bytes(&mut self, name: &str, bytes: &[u8]) {
        use sim_core::io::gds;
        let result = gds::parse(bytes).and_then(|lay| {
            let tops = lay.top_cells();
            let top = lay.default_top().ok_or(gds::GdsError::Empty)?;
            let flat = lay.flatten(Some(&top))?;
            Ok((lay, tops, top, flat))
        });
        match result {
            Ok((lay, tops, top, flat)) => {
                let layers = gds::summarize(&flat);
                if layers.is_empty() { self.gds_msg = Some((false, format!("'{name}' parsed but contains no polygon, box or path geometry in cell '{top}'."))); self.gds = None; return; }
                let best = layers.iter().enumerate().max_by_key(|(_, l)| l.polygons).map(|(i, _)| i).unwrap_or(0);
                self.gds_msg = Some((true, format!("Loaded {} polygon(s) on {} layer(s); top cell '{}'; database unit {:.3} nm.", flat.len(), layers.len(), top, lay.db_unit_m * 1e9)));
                self.gds_sel = best;
                self.gds_crop = [layers[best].anchor.0, layers[best].anchor.1];
                let layers_clone = layers.clone();
                let _ = layers_clone;
                self.gds = Some(GdsState::build(name.to_string(), Arc::new(lay), top, tops, flat));
                self.gds_view_open = true;
                self.gds_view_fit = true;
                self.status = "GDSII loaded. Pick a layer and press 'Use selected layer as Lithography mask'.".into();
            }
            Err(e) => { self.gds = None; self.gds_msg = Some((false, format!("Cannot load '{name}': {e}"))); }
        }
    }

    pub fn select_gds_layer(&mut self, i: usize) {
        if let Some(g) = &self.gds { if let Some(l) = g.layers.get(i) { self.gds_sel = i; self.gds_crop = [l.anchor.0, l.anchor.1]; } }
    }

    /// Switch the flattened top cell (when the file has several roots).
    pub fn set_gds_top(&mut self, top: &str) {
        let Some(g) = &self.gds else { return };
        match g.layout.flatten(Some(top)) {
            Ok(flat) => {
                let layers = sim_core::io::gds::summarize(&flat);
                let (file, layout, tops) = (g.file.clone(), g.layout.clone(), g.tops.clone());
                self.gds_sel = 0;
                if let Some(l) = layers.first() { self.gds_crop = [l.anchor.0, l.anchor.1]; }
                self.gds_msg = Some((true, format!("Top cell '{top}': {} polygon(s) on {} layer(s).", flat.len(), layers.len())));
                self.gds = Some(GdsState::build(file, layout, top.to_string(), tops, flat));
                self.gds_view_fit = true;
            }
            Err(e) => self.gds_msg = Some((false, format!("Cannot flatten '{top}': {e}"))),
        }
    }

    /// Put the selected layer (cropped around the crop centre) into the selected Lithography stage.
    /// Run the whole layout with hierarchical dedup, on a worker thread.
    ///
    /// Everything that can be slow — the hierarchy walk, the cap check, the etch — happens
    /// inside the thread. Doing the walk on the UI thread first made the button look dead:
    /// a 13k-polygon layout takes minutes to enumerate and the window simply froze.
    pub fn run_whole_file(&mut self) {
        if self.gds_run_rx.is_some() {
            self.gds_msg = Some((false, "A whole-file run is already going.".into()));
            return;
        }
        let Some(g) = &self.gds else { self.gds_msg = Some((false, "Load a GDSII file first.".into())); return };
        let leaf = self.session.state.spec.leaf_nm;
        let spec = sim_core::gds_sweep::CellSimSpec {
            cell_nm: leaf,
            core_nm: (leaf * 16.0, leaf * 16.0),
            margin_nm: leaf * 2.0,
            substrate_nm: 500.0,
        };
        let top = g.top.clone();
        let layout = g.layout.clone();
        let radius = spec.margin_nm;
        let (tx, rx) = std::sync::mpsc::channel();
        let progress = std::sync::Arc::new(sim_core::gds_sweep::RunProgress::default());
        self.gds_run_rx = Some(rx);
        self.gds_run_progress = Some(progress.clone());
        self.gds_msg = Some((true, "Whole-file run started — walking the hierarchy…".into()));
        self.status = "Whole-file run started.".into();
        std::thread::spawn(move || {
            const CAP: usize = 20_000;
            let n = layout.placement_count(&top);
            progress.placements.store(n, std::sync::atomic::Ordering::Relaxed);
            if n > CAP {
                progress.phase.store(sim_core::gds_sweep::RunProgress::PHASE_STOPPED, std::sync::atomic::Ordering::Relaxed);
                let _ = tx.send(Err(format!(
                    "{n} placement(s) — more than the {CAP} this run will simulate. Narrow the crop first.")));
                return;
            }
            let etch = sim_core::gds_sweep::default_etch_spec();
            let out = sim_core::gds_sweep::run_whole_file_with_progress(&layout, &top, &spec, &etch, radius, 0.0, Some(&progress));
            progress.phase.store(
                if out.is_some() { sim_core::gds_sweep::RunProgress::PHASE_DONE } else { sim_core::gds_sweep::RunProgress::PHASE_STOPPED },
                std::sync::atomic::Ordering::Relaxed);
            let _ = tx.send(match out {
                Some(r) => Ok(r),
                None => Err("The layout has no area to run, or the canvas was too large.".to_string()),
            });
        });
    }

    /// Collect a finished whole-file run, if one has arrived. Called once per frame.
    fn poll_whole_file(&mut self) {
        let Some(rx) = &self.gds_run_rx else { return };
        let Ok(res) = rx.try_recv() else { return };
        self.gds_run_rx = None;
        match res {
            Ok(r) => {
                let etched = r.canvas.depth_nm.iter().filter(|v| **v > 0.0).count();
                self.gds_msg = Some((true, format!(
                    "Whole file: {} simulation(s) for {} placement(s), {} failed. Canvas {}x{} at {:.0} nm, {} column(s) etched.",
                    r.simulations, r.placements, r.failed, r.canvas.w, r.canvas.h, r.canvas.cell_nm, etched)));
                self.status = "Whole-file run complete.".into();
            }
            Err(why) => {
                self.gds_msg = Some((false, why.clone()));
                self.status = format!("Whole-file run stopped: {why}");
            }
        }
    }

    /// The loading window, drawn while a whole-file run is going.
    ///
    /// It reads the shared counters once per frame. The window never blocks: closing it
    /// leaves the run going, and the run never waits for the window.
    fn gds_progress_window(&mut self, ctx: &egui::Context) {
        let Some(p) = self.gds_run_progress.clone() else { return };
        let mut close = false;
        egui::Window::new("Whole-file run")
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.set_min_width(360.0);
                ui.label(p.label());
                ui.add(egui::ProgressBar::new(p.fraction()).show_percentage());
                let (d, t) = (
                    p.done.load(std::sync::atomic::Ordering::Relaxed),
                    p.total.load(std::sync::atomic::Ordering::Relaxed),
                );
                if t > 0 {
                    ui.label(RichText::new(format!("{d} of {t} distinct context(s) simulated")).small().weak());
                }
                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Hide").clicked() { close = true; }
                    ui.label(RichText::new("The run continues while hidden.").small().weak());
                });
            });
        if close || (self.gds_run_rx.is_none() && self.gds_sweep_rx.is_none()) {
            self.gds_run_progress = None;
        }
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
    }

    /// Plan a sweep of the whole layout: every layer, tiled to fit the grid cap.
    ///
    /// The flatten and the plan both run on a worker: on a 13k-polygon layout the flatten
    /// alone takes long enough that the button looked dead.
    pub fn plan_whole_file_sweep(&mut self) {
        if self.gds_sweep_rx.is_some() || self.gds_run_rx.is_some() {
            self.gds_msg = Some((false, "A job is already going.".into()));
            return;
        }
        let Some(g) = &self.gds else { self.gds_msg = Some((false, "Load a GDSII file first.".into())); return };
        let layout = g.layout.clone();
        let top = g.top.clone();
        let cell = self.session.state.spec.leaf_nm;
        let cells = 192.min(sim_core::geom::MAX_GRID_AXIS);
        let (tx, rx) = std::sync::mpsc::channel();
        let progress = std::sync::Arc::new(sim_core::gds_sweep::RunProgress::default());
        self.gds_sweep_rx = Some(rx);
        self.gds_run_progress = Some(progress.clone());
        self.gds_msg = Some((true, "Planning the sweep — flattening the layout…".into()));
        self.status = "Sweep planning started.".into();
        std::thread::spawn(move || {
            let polys = match layout.flatten(Some(&top)) {
                Ok(v) => v,
                Err(e) => { let _ = tx.send(Err(format!("Cannot flatten the layout: {e}"))); return; }
            };
            progress.phase.store(sim_core::gds_sweep::RunProgress::PHASE_SIM, std::sync::atomic::Ordering::Relaxed);
            match sim_core::gds_sweep::plan(&polys, cell, cells) {
                Some(p) => {
                    progress.phase.store(sim_core::gds_sweep::RunProgress::PHASE_DONE, std::sync::atomic::Ordering::Relaxed);
                    let _ = tx.send(Ok(p));
                }
                None => {
                    progress.phase.store(sim_core::gds_sweep::RunProgress::PHASE_STOPPED, std::sync::atomic::Ordering::Relaxed);
                    let _ = tx.send(Err("The layout has no area to sweep.".into()));
                }
            }
        });
    }

    /// Collect a finished sweep plan, if one has arrived. Called once per frame.
    fn poll_sweep(&mut self) {
        let Some(rx) = &self.gds_sweep_rx else { return };
        let Ok(res) = rx.try_recv() else { return };
        self.gds_sweep_rx = None;
        match res {
            Ok(p) => {
                let (w, h) = (p.extent.width(), p.extent.height());
                let summary = format!(
                    "Whole-file sweep: {:.1} × {:.1} µm, {} tiles of {:.1} µm — {} hold geometry, {} are empty and skipped. {} layer(s), {} cells per tile.",
                    w / 1e3, h / 1e3, p.tiles(), p.grid.tile_nm / 1e3, p.active, p.empty, p.layer_count, p.cells);
                self.gds_sweep = Some(p);
                if let Some(g) = &self.gds {
                    let d = g.layout.dedup();
                    self.gds_msg = Some((true, format!(
                        "{summary} Hierarchy: {} cell(s), {} distinct content, {} placement(s) — {} saved by dedup ({:.0}%).",
                        d.cells, d.unique, d.placements, d.saved, d.saved_fraction() * 100.0)));
                } else {
                    self.gds_msg = Some((true, summary));
                }
                self.status = "Sweep planned. Each tile is one simulation; the window is capped at 192 cells per axis.".into();
            }
            Err(why) => {
                self.gds_msg = Some((false, why.clone()));
                self.status = format!("Sweep planning stopped: {why}");
            }
        }
    }

    pub fn use_gds_layer(&mut self) {
        let Some(g) = &self.gds else { self.gds_msg = Some((false, "Load a GDSII file first.".into())); return };
        let Some(l) = g.layers.get(self.gds_sel).cloned() else { return };
        let Some(i) = self.session.state.recipe.index_of(self.selected_stage) else { return };
        let mut p = self.session.state.recipe.stages[i].params.clone();
        let sim_core::params::StageParams::Lithography(lp) = &mut p else { self.gds_msg = Some((false, "Select a Lithography stage in the recipe first.".into())); return };
        let clear_field = match &lp.mask { sim_core::geom::MaskPattern::LineSpace { clear_field, .. } | sim_core::geom::MaskPattern::Rects { clear_field, .. } | sim_core::geom::MaskPattern::Polygons { clear_field, .. } => *clear_field };
        let half = self.session.state.recipe.window.size[0].max(self.session.state.recipe.window.size[1]) * 0.5 + 500.0;
        // Expand only the neighbourhood of the crop centre (instances far away are skipped without being expanded).
        let region = [self.gds_crop[0] - half, self.gds_crop[1] - half, self.gds_crop[0] + half, self.gds_crop[1] + half];
        let near = match g.layout.flatten_region(Some(&g.top), Some(region)) { Ok(v) => v, Err(e) => { self.gds_msg = Some((false, format!("Cannot read the crop region: {e}"))); return } };
        let all = self.gds_all_layers;
        let built = if all {
            sim_core::geom::MaskPattern::from_layout_all(&near, (self.gds_crop[0], self.gds_crop[1]), half, clear_field)
        } else {
            sim_core::geom::MaskPattern::from_layout(&near, l.layer, Some(l.datatype), (self.gds_crop[0], self.gds_crop[1]), half, clear_field)
        };
        match built {
            Ok((m, kept)) => {
                lp.mask = m;
                let _ = self.session.edit_params(self.selected_stage, p);
                let (w, h) = l.size_nm();
                let note = if w > 2.0 * half || h > 2.0 * half { format!(" The layer spans {:.1} × {:.1} µm but only ±{:.0} nm around the crop centre is used.", w / 1e3, h / 1e3, half) } else { String::new() };
                self.gds_msg = Some((true, format!("Mask = layer {}/{}: {kept} polygon(s) near ({:.0}, {:.0}) nm.{note}", l.layer, l.datatype, self.gds_crop[0], self.gds_crop[1])));
                self.status = "GDS mask applied to the Lithography stage; results are stale until you Run.".into();
            }
            Err(e) => self.gds_msg = Some((false, e.to_string())),
        }
    }

    /// Start the fast verification suite (32/16 nm cells) without blocking the UI.
    pub(crate) fn start_verification(&mut self) {
        if self.vv_rx.is_some() { return; }
        let cfg = sim_core::vv::VvConfig { leaf_sizes_nm: vec![32.0, 16.0], ..Default::default() };
        self.status = "Running verification (conservation, exact solutions, FEM)…".into();
        #[cfg(not(target_arch = "wasm32"))]
        {
            let (tx, rx) = std::sync::mpsc::channel();
            std::thread::spawn(move || { let _ = tx.send(sim_core::vv::run_suite(&cfg)); });
            self.vv_rx = Some(rx);
        }
        #[cfg(target_arch = "wasm32")]
        { self.vv_report = Some(sim_core::vv::run_suite(&cfg)); }
    }

    fn timeline_panel(&mut self, ui: &mut egui::Ui) {
        let n = self.timeline.segments.len();
        // Stage chips: simulation state (derived) for every recipe stage.
        ui.horizontal_wrapped(|ui| {
            let stages = self.session.state.recipe.stages.clone();
            let statuses = self.session.statuses();
            for (i, (rs, (_, st))) in stages.iter().zip(statuses.iter()).enumerate() {
                let sel = self.selected_stage == rs.id;
                let text = RichText::new(format!("{} {}", i + 1, rs.params.stage().label())).color(status_color(st));
                if ui.selectable_label(sel, text).on_hover_text(status_word(st)).clicked() { self.selected_stage = rs.id; }
                if i + 1 < stages.len() { ui.label("›"); }
            }
        });
        ui.horizontal(|ui| {
            let has = n > 0;
            if ui.add_enabled(has, egui::Button::new("⏮")).on_hover_text("Reset to the start").clicked() { self.playback.reset(); }
            if ui.add_enabled(has, egui::Button::new("⏪")).on_hover_text("Step back one keyframe").clicked() { self.playback.step_back(); }
            let label = match self.playback.state { PlayState::Playing => "⏸ Pause", PlayState::Paused => "▶ Resume", PlayState::Stopped => "▶ Play" };
            if ui.add_enabled(has, egui::Button::new(label)).on_hover_text("Space").clicked() { self.playback.toggle(); }
            if ui.add_enabled(has, egui::Button::new("⏩")).on_hover_text("Step forward one keyframe").clicked() { self.playback.step_forward(); }
            if ui.add_enabled(has, egui::Button::new("⏭")).clicked() { self.playback.seek(n as f64); self.playback.state = PlayState::Stopped; }
            ui.add(egui::Slider::new(&mut self.playback.speed, 0.25..=8.0).logarithmic(true).text("speed ×"));
            ui.separator();
            let total = self.timeline.total_duration_s();
            let cur = self.timeline.sim_time_at(self.playback.position);
            ui.label(format!("process time {} / {}", fmt_time(cur), fmt_time(total)));
        });
        ui.horizontal(|ui| {
            let mut pos = self.playback.position;
            ui.spacing_mut().slider_width = (ui.available_width() - 20.0).max(100.0);
            if ui.add_enabled(n > 0, egui::Slider::new(&mut pos, 0.0..=(n.max(1)) as f64).show_value(false)).changed() { self.playback.seek(pos); if self.playback.is_playing() { self.playback.pause(); } }
        });
        // Simulation progress vs animation progress are different things; show both.
        let enabled = self.session.enabled_count().max(1);
        let committed = self.session.committed_count();
        ui.columns(2, |c| {
            c[0].label(RichText::new("Simulation (computed & committed)").small());
            c[0].add(egui::ProgressBar::new(committed as f32 / enabled as f32).text(format!("{committed}/{enabled} stages · rev {}", self.session.state.revision)));
            c[1].label(RichText::new("Animation (replay of recorded stages)").small());
            let frac = if n == 0 { 0.0 } else { (self.playback.position / n as f64) as f32 };
            let phase = self.display.phase.map(|p| format!("{} [{}]", p.name, p.kind.label())).unwrap_or_else(|| "—".into());
            c[1].add(egui::ProgressBar::new(frac).text(phase));
        });
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.show_log, "Process log");
            ui.label(RichText::new("Frames are replayed deterministically from the recorded inputs; the last frame of a stage is its committed geometry.").weak().small());
        });
        if self.show_log {
            egui::ScrollArea::vertical().stick_to_bottom(true).max_height(80.0).show(ui, |ui| {
                for l in self.session.log.iter() { ui.label(RichText::new(l).monospace().small()); }
            });
        }
    }
}

fn fmt_time(s: f64) -> String { let t = s.max(0.0) as u64; format!("{:02}:{:02}:{:02}", t / 3600, (t / 60) % 60, t % 60) }

impl eframe::App for AttoFabApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) { self.poll_whole_file(); self.poll_sweep(); self.gds_progress_window(ctx); self.ui(ctx); }
}

impl AttoFabApp {
    /// One full UI pass (what `eframe::App::update` runs).
    pub fn ui(&mut self, ctx: &egui::Context) {
        if let Some(rx) = &self.vv_rx {
            match rx.try_recv() {
                Ok(r) => { self.status = format!("Verification finished: {} case(s), {} failed.", r.cases.len(), r.failed().len()); self.vv_report = Some(r); self.vv_rx = None; }
                Err(std::sync::mpsc::TryRecvError::Disconnected) => { self.status = "Verification worker stopped unexpectedly.".into(); self.vv_rx = None; }
                Err(std::sync::mpsc::TryRecvError::Empty) => ctx.request_repaint_after(std::time::Duration::from_millis(200)),
            }
        }
        let dropped = ctx.input(|i| i.raw.dropped_files.clone());
        for f in dropped {
            if let Some(bytes) = &f.bytes { self.load_gds_bytes(&f.name, bytes); }
            #[cfg(not(target_arch = "wasm32"))]
            if f.bytes.is_none() { if let Some(p) = &f.path { self.gds_path = p.display().to_string(); self.load_gds_path(); } }
        }
        self.apply_theme(ctx);
        self.shortcuts(ctx);
        let dt = ctx.input(|i| i.stable_dt) as f64;
        if self.playback.advance(dt) || self.playback.is_playing() { ctx.request_repaint(); }
        self.sync_timeline();
        self.refresh_display();
        if self.session.state.spec.leaf_nm != self.scale_synced.0 || self.session.state.spec.tile_levels != self.scale_synced.1 { self.scale_synced = (self.session.state.spec.leaf_nm, self.session.state.spec.tile_levels); self.scale_edit = self.scale_synced; }

        egui::TopBottomPanel::top("toolbar").show(ctx, |ui| self.toolbar(ui));
        egui::TopBottomPanel::bottom("timeline").resizable(true).default_height(190.0).show(ctx, |ui| self.timeline_panel(ui));
        egui::SidePanel::left("recipe").resizable(true).default_width(230.0).show(ctx, |ui| self.recipe_panel(ui));
        egui::SidePanel::right("inspector").resizable(true).default_width(340.0).show(ctx, |ui| self.inspector_panel(ui));
        egui::CentralPanel::default().frame(egui::Frame::none()).show(ctx, |ui| self.viewport(ui));
        self.gds_viewer(ctx);
    }
}

#[cfg(test)]
mod ui_tests {
    //! End-to-end input tests: synthetic egui pointer events drive the real app
    //! (tool modes, drag capture, locks, undo) with no GPU.
    use super::*;
    use eframe::egui::{Event, Modifiers, PointerButton, Pos2, RawInput, Rect};
    use crate::viewport_math::pointer_to_ndc;

    struct Harness { ctx: egui::Context, app: AttoFabApp, t: f64, pos: Pos2, mods: Modifiers }

    impl Harness {
        fn new() -> Self {
            let ctx = egui::Context::default();
            let app = AttoFabApp::new_headless(&ctx);
            let mut h = Self { ctx, app, t: 0.0, pos: Pos2::new(10.0, 10.0), mods: Modifiers::NONE };
            for _ in 0..3 { h.frame(vec![]); }
            h
        }
        fn frame(&mut self, events: Vec<Event>) {
            self.t += 1.0 / 60.0;
            let input = RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1500.0, 950.0))), time: Some(self.t), predicted_dt: 1.0 / 60.0, events, modifiers: self.mods, ..Default::default() };
            let app = &mut self.app;
            let _ = self.ctx.run(input, |ctx| app.ui(ctx));
        }
        fn mv(&mut self, p: Pos2) { self.pos = p; self.frame(vec![Event::PointerMoved(p)]); }
        fn press(&mut self, b: PointerButton) { self.frame(vec![Event::PointerButton { pos: self.pos, button: b, pressed: true, modifiers: Modifiers::NONE }]); }
        fn release(&mut self, b: PointerButton) { self.frame(vec![Event::PointerButton { pos: self.pos, button: b, pressed: false, modifiers: Modifiers::NONE }]); }
        fn drag(&mut self, b: PointerButton, from: Pos2, to: Pos2, steps: usize) {
            self.mv(from); self.press(b);
            for i in 1..=steps { let f = i as f32 / steps as f32; self.mv(from + (to - from) * f); }
            self.release(b); self.frame(vec![]);
        }
        fn centre(&self) -> Pos2 { self.app.viewport_rect.center() }
    }

    #[test]
    fn every_tool_view_and_panel_renders_without_panicking() {
        let mut h = Harness::new();
        h.app.do_run();
        for _ in 0..5 { h.frame(vec![]); }
        for m in ViewMode::ALL { for t in Tool::ALL {
            h.app.rig.mode = m; h.app.set_tool(t);
            for _ in 0..2 { h.frame(vec![]); }
        } }
        for k in 0..=FRAMES_PER_STAGE { h.app.playback.seek(2.0 + k as f64 / FRAMES_PER_STAGE as f64); h.frame(vec![]); }
        h.app.dark = false; h.frame(vec![]); h.frame(vec![]);
        assert!(h.app.viewport_rect.width() > 300.0);
    }

    #[test]
    fn int01_move_wafer_mode_drag_follows_the_pointer_and_never_rewrites_the_svo() {
        let mut h = Harness::new();
        h.app.rig.mode = ViewMode::Top;
        h.app.fit_wafer();
        h.frame(vec![]);
        h.app.set_tool(Tool::MoveWafer);
        let wafer_ptr = Arc::as_ptr(&h.app.session.committed_wafer());
        let fp = h.app.session.committed_wafer().fingerprint();
        let depth = h.app.session.undo_depth();
        let c = h.centre();
        // Grab a point left of the wafer centre so it is on the disk.
        let from = c + Vec2::new(-60.0, 20.0);
        let to = from + Vec2::new(140.0, -70.0);
        h.mv(from); h.press(PointerButton::Primary);
        let aspect = h.app.viewport_aspect;
        let (mn, sz) = ([h.app.viewport_rect.min.x, h.app.viewport_rect.min.y], [h.app.viewport_rect.width(), h.app.viewport_rect.height()]);
        let grab_local = {
            let n = pointer_to_ndc([from.x, from.y], mn, sz).unwrap();
            let plane = h.app.work_plane_for_test();
            let hit = crate::viewport_math::pick_on_wafer(&h.app.rig, n, aspect, &plane, &h.app.transform(), h.app.wafer_radius()).unwrap();
            hit.local
        };
        for i in 1..=12 { let f = i as f32 / 12.0; h.mv(from + (to - from) * f); }
        // While still dragging, the grabbed wafer point is under the pointer.
        let under = h.app.rig.project(h.app.transform().local_to_world(grab_local), aspect).unwrap();
        let want = pointer_to_ndc([to.x, to.y], mn, sz).unwrap();
        assert!((under[0] - want[0]).abs() < 1e-6 && (under[1] - want[1]).abs() < 1e-6, "grabbed point drifted: {under:?} vs {want:?}");
        h.release(PointerButton::Primary); h.frame(vec![]);
        assert!(h.app.transform().position[0].abs() + h.app.transform().position[1].abs() > 1.0e5, "the wafer moved");
        assert_eq!(h.app.transform().position[2], 0.0, "stays horizontal / at height");
        assert_eq!(Arc::as_ptr(&h.app.session.committed_wafer()), wafer_ptr, "SVO not rewritten");
        assert_eq!(h.app.session.committed_wafer().fingerprint(), fp);
        assert_eq!(h.app.session.undo_depth(), depth + 1);
        assert!(h.app.session.undo());
        assert_eq!(h.app.transform(), WaferTransform::IDENTITY);
    }

    #[test]
    fn int02_navigate_mode_drag_moves_the_camera_not_the_wafer() {
        let mut h = Harness::new();
        h.app.session.set_transform(WaferTransform { position: [3.0e6, 1.0e6, 0.0], yaw_rad: 0.2 });
        h.frame(vec![]);
        let (t0, yaw0) = (h.app.transform(), h.app.rig.yaw);
        h.app.set_tool(Tool::Navigate);
        let c = h.centre();
        h.drag(PointerButton::Primary, c, c + Vec2::new(90.0, 30.0), 8);
        assert_eq!(h.app.transform(), t0);
        assert!(h.app.rig.yaw != yaw0, "3D left-drag orbits");
        h.app.rig.mode = ViewMode::Top;
        let target0 = h.app.rig.target;
        h.drag(PointerButton::Primary, c, c + Vec2::new(60.0, 0.0), 6);
        assert!(h.app.rig.target != target0, "top-view left-drag pans");
        assert_eq!(h.app.transform(), t0);
        h.drag(PointerButton::Middle, c, c + Vec2::new(-40.0, 20.0), 5);
        assert_eq!(h.app.transform(), t0, "middle-drag pans in any tool");
    }

    #[test]
    fn place_mode_drag_never_moves_the_wafer_and_only_proposes_a_site() {
        let mut h = Harness::new();
        h.app.rig.mode = ViewMode::Top;
        h.app.fit_wafer();
        h.frame(vec![]);
        h.app.set_tool(Tool::Place);
        let c = h.centre();
        h.drag(PointerButton::Primary, c + Vec2::new(-30.0, 10.0), c + Vec2::new(50.0, 40.0), 6);
        assert_eq!(h.app.transform(), WaferTransform::IDENTITY);
        assert!(h.app.session.pending_site.is_some(), "the click/drag proposed a site");
        assert_eq!(h.app.session.committed_count(), 0);
        h.frame(vec![Event::Key { key: egui::Key::Escape, physical_key: None, pressed: true, repeat: false, modifiers: Modifiers::NONE }]);
        assert!(h.app.session.pending_site.is_none(), "Esc cancels the preview");
    }

    #[test]
    fn dragging_outside_the_viewport_or_while_a_process_plays_never_moves_the_wafer() {
        let mut h = Harness::new();
        h.app.rig.mode = ViewMode::Top; h.app.fit_wafer(); h.frame(vec![]);
        h.app.set_tool(Tool::MoveWafer);
        // Press on the right-hand inspector panel (a UI panel), drag across the viewport.
        let start = Pos2::new(1400.0, 300.0);
        assert!(!h.app.viewport_rect.contains(start));
        h.drag(PointerButton::Primary, start, h.centre(), 8);
        assert_eq!(h.app.transform(), WaferTransform::IDENTITY, "a drag that starts on a panel must not move the wafer");
        // While a process animation plays the wafer is locked.
        h.app.do_run();
        for _ in 0..3 { h.frame(vec![]); }
        assert!(h.app.playback.is_playing());
        let c = h.centre();
        h.drag(PointerButton::Primary, c, c + Vec2::new(80.0, 40.0), 6);
        assert_eq!(h.app.transform(), WaferTransform::IDENTITY, "wafer movement is locked while running");
        h.app.playback.pause();
        h.drag(PointerButton::Primary, c, c + Vec2::new(80.0, 40.0), 6);
        assert!(h.app.transform() != WaferTransform::IDENTITY, "unlocked once paused");
    }

    fn real_gds_path() -> String { format!("{}/../sim_core/tests/data/attofab_master_wafer.gds", env!("CARGO_MANIFEST_DIR")) }

    #[test]
    fn real_layout_loads_shows_in_the_viewer_and_a_click_chooses_the_site() {
        let mut h = Harness::new();
        h.app.gds_path = real_gds_path();
        h.app.load_gds_path();
        let (n_layers, polys) = { let g = h.app.gds.as_ref().unwrap_or_else(|| panic!("load failed: {:?}", h.app.gds_msg)); (g.layers.len(), g.flat.len()) };
        assert_eq!((n_layers, polys > 1000), (18, true));
        assert!(h.app.gds_view_open, "loading opens the layout viewer");
        for _ in 0..3 { h.frame(vec![]); }
        let canvas = h.app.gds_canvas_rect;
        assert!(canvas.width() > 150.0 && canvas.height() > 100.0, "viewer canvas was laid out: {canvas:?}");
        // Click the middle of the canvas: the crop centre becomes the world point there.
        let before = h.app.gds_crop;
        let want = h.app.gds_view.to_world([canvas.center().x, canvas.center().y], [canvas.center().x, canvas.center().y]);
        h.mv(canvas.center()); h.press(PointerButton::Primary); h.release(PointerButton::Primary); h.frame(vec![]);
        let tol = h.app.gds_view.nm_per_px * 3.0;
        assert!((h.app.gds_crop[0] - want[0]).abs() <= tol && (h.app.gds_crop[1] - want[1]).abs() <= tol, "crop {:?} vs clicked {:?} (was {:?})", h.app.gds_crop, want, before);
        // Zooming with the wheel keeps working and never panics; fit/zoom-to-crop buttons exist via state flags.
        h.frame(vec![Event::MouseWheel { unit: egui::MouseWheelUnit::Point, delta: Vec2::new(0.0, 120.0), modifiers: Modifiers::NONE }]);
        h.frame(vec![]);
    }

    #[test]
    fn auto_scale_picks_a_fine_cell_for_nm_features_and_a_coarse_one_for_um_features_and_runs() {
        let mut h = Harness::new();
        h.app.gds_path = real_gds_path();
        h.app.load_gds_path();
        let lith = h.app.session.state.recipe.stages.iter().find(|s| s.params.stage() == Stage::Lithography).unwrap().id;
        h.app.selected_stage = lith;
        let pick = |h: &mut Harness, layer: i32, dt: i32| { let i = h.app.gds.as_ref().unwrap().layers.iter().position(|l| l.layer == layer && l.datatype == dt).unwrap(); h.app.select_gds_layer(i); };
        pick(&mut h, 2, 0);
        assert!(h.app.jump_to_finest_feature());
        h.app.auto_scale_from_gds();
        let fine = h.app.session.state.spec.leaf_nm;
        assert!(matches!(&h.app.session.state.recipe.stages.iter().find(|s| s.id == lith).unwrap().params, sim_core::params::StageParams::Lithography(l) if matches!(l.mask, sim_core::geom::MaskPattern::Polygons { .. })), "{:?}", h.app.gds_msg);
        pick(&mut h, 1, 4);
        h.app.auto_scale_from_gds();
        let coarse = h.app.session.state.spec.leaf_nm;
        assert!(fine < coarse && fine <= 16.0 && coarse >= 25.0 && coarse <= 64.0, "cell for 50 nm features {fine} nm, for µm features {coarse} nm");
        assert!(h.app.session.state.recipe.window.size[0] <= h.app.session.max_window_nm() + 1e-9);
        // Run only up to lithography (keeps the test fast) and require a real pattern.
        let ids: Vec<u32> = h.app.session.state.recipe.stages.iter().map(|s| s.id).collect();
        for id in &ids[3..] { h.app.session.set_enabled(*id, false).unwrap(); }
        let r = h.app.session.run_recipe();
        assert!(r.failed.is_none(), "{:?}", r.failed);
        let rec = &h.app.session.state.records[&lith];
        let (opened, retained) = (rec.outcome.metric("Opened columns").unwrap(), rec.outcome.metric("Retained columns").unwrap());
        assert!(opened > 0.0 && retained > 0.0, "opened {opened}, retained {retained}");
        for _ in 0..3 { h.frame(vec![]); }
    }

    /// Load the master wafer, pick layer 2/0, jump to its finest feature and auto-scale, so the
    /// lithography stage holds a polygon mask built from the layout. Returns the lithography stage id.
    fn harness_with_gds_mask() -> (Harness, u32) {
        let mut h = Harness::new();
        h.app.gds_path = real_gds_path();
        h.app.load_gds_path();
        let lith = h.app.session.state.recipe.stages.iter().find(|s| s.params.stage() == Stage::Lithography).unwrap().id;
        h.app.selected_stage = lith;
        let i = h.app.gds.as_ref().unwrap().layers.iter().position(|l| l.layer == 2 && l.datatype == 0).unwrap();
        h.app.select_gds_layer(i);
        assert!(h.app.jump_to_finest_feature());
        h.app.auto_scale_from_gds();
        (h, lith)
    }

    fn litho_mask(h: &Harness, lith: u32) -> sim_core::geom::MaskPattern {
        match &h.app.session.state.recipe.stages.iter().find(|s| s.id == lith).unwrap().params {
            sim_core::params::StageParams::Lithography(l) => l.mask.clone(),
            other => panic!("stage {lith} is not lithography: {other:?}"),
        }
    }

    #[test]
    fn moving_the_crop_off_the_pattern_says_the_mask_was_not_rebuilt_instead_of_leaving_it_silently() {
        let (mut h, lith) = harness_with_gds_mask();
        let before = litho_mask(&h, lith);
        assert!(matches!(before, sim_core::geom::MaskPattern::Polygons { .. }), "{:?}", h.app.gds_msg);
        h.app.gds_crop = [1.0e9, 1.0e9];
        let note = h.app.rebuild_mask_for_crop().expect("a rebuild was attempted, so the outcome is reported");
        assert!(note.contains("NOT rebuilt"), "{note}");
        assert!(matches!(&h.app.gds_msg, Some((false, m)) if m.contains("NOT rebuilt")), "{:?}", h.app.gds_msg);
        assert_eq!(litho_mask(&h, lith), before, "a failed rebuild leaves the previous mask untouched");
    }

    #[test]
    fn a_rebuilt_mask_that_changes_the_recipe_says_so_and_lithography_is_never_reported_up_to_date() {
        let (mut h, lith) = harness_with_gds_mask();
        let ids: Vec<u32> = h.app.session.state.recipe.stages.iter().map(|s| s.id).collect();
        for id in &ids[3..] { h.app.session.set_enabled(*id, false).unwrap(); }
        let r = h.app.session.run_recipe();
        assert!(r.failed.is_none(), "{:?}", r.failed);
        assert_eq!(h.app.session.status_of(lith), Some(sim_core::session::StageStatus::Committed));
        h.app.gds_all_layers = true;
        let note = h.app.rebuild_mask_for_crop().expect("a rebuild was attempted, so the outcome is reported");
        assert!(matches!(&h.app.gds_msg, Some((true, _))), "{:?}", h.app.gds_msg);
        let stale = h.app.session.status_of(lith) == Some(sim_core::session::StageStatus::Stale);
        assert_eq!(note.contains("stale"), stale, "the note and the stage status must agree: {note}");
        let r = h.app.session.run_recipe();
        if stale { assert!(!r.skipped_committed.contains(&lith), "a stale lithography stage must be recomputed, never 'already up to date'"); }
    }

    #[test]
    fn left_panel_with_the_full_gds_section_renders_inside_a_scroll_area() {
        let mut h = Harness::new();
        h.app.gds_path = real_gds_path();
        h.app.load_gds_path();
        for _ in 0..4 { h.frame(vec![]); }
        // The recipe + site + scale + GDS sections are taller than the window area; they live in a scroll area (no panic, no overlap).
        assert!(h.app.gds_msg.as_ref().map(|m| m.0).unwrap_or(false));
    }

    #[test]
    fn select_tool_picks_a_leaf_and_the_marker_is_large_enough_to_see() {
        let mut h = Harness::new();
        h.app.do_run();
        h.app.playback.state = PlayState::Stopped; h.app.playback.seek(7.0);
        h.app.rig.mode = ViewMode::Top; h.app.focus_site(); h.frame(vec![]);
        h.app.set_tool(Tool::Select);
        let c = h.centre();
        h.mv(c); h.press(PointerButton::Primary); h.release(PointerButton::Primary); h.frame(vec![]);
        let (b, m) = h.app.selection.expect("a leaf was selected");
        assert!(m != Material::Air && b.max.0 > b.min.0);
        assert!(h.app.status.starts_with("Selected"), "{}", h.app.status);
        h.frame(vec![]); h.frame(vec![]); // overlay + inspector section render without panicking
        h.app.select_at_ndc([5.0, 5.0]);
        assert!(h.app.selection.is_none(), "clicking empty space clears the selection");
    }

    #[test]
    fn inspect_reports_the_committed_stack_under_the_cursor() {
        let mut h = Harness::new();
        h.app.do_run();
        h.app.playback.state = PlayState::Stopped; h.app.playback.seek(7.0);
        h.app.rig.mode = ViewMode::Top; h.app.focus_site(); h.frame(vec![]);
        h.app.set_tool(Tool::Inspect);
        // Click the window centre: between lines the stack is oxide-free silicon with dopants, or oxide under resist-free line.
        let c = h.centre();
        h.mv(c); h.press(PointerButton::Primary); h.release(PointerButton::Primary); h.frame(vec![]);
        let i = h.app.inspect.as_ref().expect("inspect readout");
        assert!(i.in_tile && !i.runs.is_empty());
        assert!(i.runs.iter().any(|r| r.material == Material::Silicon || matches!(r.material, Material::Dopant(_))));
    }

    #[test]
    fn undo_redo_shortcuts_and_pause_resume_keep_animation_and_state_in_sync() {
        let mut h = Harness::new();
        h.app.do_run();
        for _ in 0..4 { h.frame(vec![]); }
        assert!(h.app.playback.is_playing());
        h.app.playback.pause();
        let pos = h.app.playback.position;
        for _ in 0..10 { h.frame(vec![]); }
        assert_eq!(h.app.playback.position, pos, "paused animation does not advance");
        assert_eq!(h.app.display.kind, DisplayKind::AnimationFrame);
        assert_eq!(h.app.session.committed_count(), 7, "the simulation state is unaffected by the replay position");
        h.app.playback.resume();
        for _ in 0..10 { h.frame(vec![]); }
        assert!(h.app.playback.position > pos);
        h.app.playback.seek(7.0); h.app.playback.state = PlayState::Stopped; h.frame(vec![]);
        assert_eq!(h.app.display.kind, DisplayKind::Committed);
        h.mods = Modifiers::COMMAND;
        h.frame(vec![Event::Key { key: egui::Key::Z, physical_key: None, pressed: true, repeat: false, modifiers: Modifiers::COMMAND }]);
        assert_eq!(h.app.session.committed_count(), 0, "Ctrl+Z undoes the whole run");
        h.frame(vec![Event::Key { key: egui::Key::Y, physical_key: None, pressed: true, repeat: false, modifiers: Modifiers::COMMAND }]);
        assert_eq!(h.app.session.committed_count(), 7);
    }

    #[test]
    fn verification_runs_off_the_ui_thread_and_reports_back() {
        let mut h = Harness::new();
        h.app.start_verification();
        assert!(h.app.vv_rx.is_some(), "worker started");
        for _ in 0..600 {
            h.frame(vec![]);
            if h.app.vv_report.is_some() { break; }
            std::thread::sleep(std::time::Duration::from_millis(50));
        }
        let r = h.app.vv_report.as_ref().expect("verification finished");
        assert!(r.failed().is_empty(), "{:?}", r.failed().iter().map(|c| &c.id).collect::<Vec<_>>());
        assert!(h.app.vv_rx.is_none());
        h.frame(vec![]);
    }

    #[test]
    fn gds_hierarchical_file_lists_layers_and_loads_as_the_lithography_mask() {
        let mut h = Harness::new();
        let bytes = sim_core::io::gds::tests_support::sample_stream();
        let path = std::env::temp_dir().join("attofab test mask.gds");
        std::fs::write(&path, &bytes).unwrap();
        // A pasted Windows-style quoted path with stray spaces must work.
        h.app.gds_path = format!("  \"{}\"  ", path.display());
        h.app.load_gds_path();
        let g = h.app.gds.as_ref().unwrap_or_else(|| panic!("load failed: {:?}", h.app.gds_msg));
        assert_eq!(g.top, "TOP");
        assert!(g.layers.iter().any(|l| l.layer == 5 && l.polygons == 7), "hierarchy (SREF + AREF) flattened: {:?}", g.layers.iter().map(|l| (l.layer, l.polygons)).collect::<Vec<_>>());
        assert!(g.layers.iter().any(|l| l.layer == 7) && g.layers.iter().any(|l| l.layer == 9), "paths and boxes are listed");
        let idx = g.layers.iter().position(|l| l.layer == 5).unwrap();
        h.app.select_gds_layer(idx);
        let lith = h.app.session.state.recipe.stages.iter().find(|s| s.params.stage() == Stage::Lithography).unwrap().id;
        h.app.selected_stage = lith;
        h.app.use_gds_layer();
        assert!(matches!(&h.app.session.state.recipe.stages.iter().find(|s| s.id == lith).unwrap().params, sim_core::params::StageParams::Lithography(l) if matches!(l.mask, sim_core::geom::MaskPattern::Polygons { .. })), "{:?}", h.app.gds_msg);
        h.frame(vec![]); h.frame(vec![]);
        // Errors are explicit, never silent.
        h.app.gds_path = "C:\\definitely\\missing.gds".into(); h.app.load_gds_path();
        assert!(matches!(&h.app.gds_msg, Some((false, m)) if m.contains("Cannot read")));
        h.app.load_gds_bytes("notes.txt", b"this is not a layout at all, sorry");
        assert!(matches!(&h.app.gds_msg, Some((false, m)) if m.contains("not a GDSII")), "{:?}", h.app.gds_msg);
        assert_eq!(AttoFabApp::normalize_path("\"C:\\a b\\c.gds\""), "C:\\a b\\c.gds");
    }

    #[test]
    fn preview_shows_an_uncommitted_badge_and_cancel_restores_the_committed_view() {
        let mut h = Harness::new();
        h.app.do_preview();
        for _ in 0..2 { h.frame(vec![]); }
        h.app.playback.state = PlayState::Stopped; h.app.playback.seek(1.0); h.frame(vec![]);
        assert_eq!(h.app.display.kind, DisplayKind::Preview);
        assert_eq!(h.app.session.committed_count(), 0);
        h.app.do_cancel(); h.frame(vec![]); h.frame(vec![]);
        assert_eq!(h.app.display.kind, DisplayKind::Committed);
        assert_eq!(h.app.display.wafer.fingerprint(), h.app.session.committed_wafer().fingerprint());
    }
}
