// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The the crates panels, drawn with egui.
//!
//! **This file has never been compiled.** `eframe`/`egui` need X11 and Wayland
//! headers that cannot be installed without root in the environment where this
//! was written. The rendering is a proposal; the *logic* behind it is not —
//! every value shown comes from `atto_ui`, which has 40 tests and no toolkit
//! dependency.
//!
//! The shape is deliberate. Each panel is a function of `&Session` and the `ui`,
//! so nothing here holds state. A panel that cannot derive what it needs from the
//! session has to say so rather than caching a copy that can go stale.
//!
//! Written against `eframe 0.27` (the pin in the workspace): `App::update(&mut
//! self, ctx, frame)`, `egui::Grid`, `ProgressBar`, and `allocate_painter` for
//! the plots, so no extra plotting crate is needed.

use atto_ui::doe::DoePanel;
use atto_ui::inspector::{Inspector, InspectorTab};
use atto_ui::plots::{cv_series, iv_series, Series};
use atto_ui::stack::{solid_thickness_nm, stack_rows};
use atto_ui::{measure_rows, palette, Colour, Progress, Session};
use eframe::egui;

pub struct Mk2Shell {
    session: Session,
    journal: atto_ui::journal::Journal,
    inspector: Inspector,
    progress: Option<Progress>,
    doe: Option<DoePanel>,
}

impl Default for Mk2Shell {
    fn default() -> Self {
        Mk2Shell {
            session: Session::new(),
            journal: atto_ui::journal::Journal::new(),
            inspector: Inspector::new(),
            progress: None,
            doe: None,
        }
    }
}

impl eframe::App for Mk2Shell {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::SidePanel::left("recipe").resizable(true).show(ctx, |ui| {
            recipe_panel(ui, &mut self.session);
        });
        egui::TopBottomPanel::bottom("status").show(ctx, |ui| {
            status_bar(ui, &self.session, self.progress.as_ref());
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            inspector_tabs(ui, &mut self.inspector, &self.session);
        });
    }
}

fn colour_of(c: Colour) -> egui::Color32 {
    match c {
        Colour::Grey => egui::Color32::GRAY,
        Colour::Green => egui::Color32::from_rgb(80, 160, 80),
        Colour::Amber => egui::Color32::from_rgb(200, 160, 40),
        Colour::Red => egui::Color32::from_rgb(200, 80, 80),
    }
}

/// The recipe list. Status is derived from the session, never stored (G-02), and
/// always shows the word beside the colour (G-04).
pub fn recipe_panel(ui: &mut egui::Ui, session: &mut Session) {
    ui.heading("Recipe");
    egui::ScrollArea::vertical().show(ui, |ui| {
        for row in session.recipe_rows() {
            ui.horizontal(|ui| {
                let text = format!("{}. {}", row.index + 1, row.label);
                if ui.selectable_label(row.selected, text).clicked() {
                    let _ = session.select(row.index);
                }
                ui.colored_label(colour_of(row.status.colour()), row.status.word());
                // G-07: a failure shows its own message, not just a red word
                if let Some(why) = session.failure_reason(row.index) {
                    ui.label("?").on_hover_text(why);
                }
            });
        }
    });
}

/// The five inspector tabs, each keeping its own selection and scroll (G-11).
pub fn inspector_tabs(ui: &mut egui::Ui, inspector: &mut Inspector, session: &Session) {
    ui.horizontal(|ui| {
        for (label, current) in inspector.tabs() {
            let tab = InspectorTab::ALL.iter().find(|t| t.label() == label).copied();
            if let Some(tab) = tab {
                if ui.selectable_label(current, label).clicked() {
                    inspector.select_tab(tab);
                }
            }
        }
    });
    ui.separator();
    egui::ScrollArea::vertical().vertical_scroll_offset(inspector.state().scroll as f32).show(ui, |ui| {
        match inspector.current() {
            InspectorTab::Stage => stage_tab(ui, session),
            InspectorTab::Stack => stack_tab(ui, session),
            InspectorTab::Measure => measure_tab(ui, session),
            InspectorTab::Device => device_tab(ui),
            InspectorTab::Evidence => evidence_tab(ui, session),
        }
    });
}

fn stage_tab(ui: &mut egui::Ui, session: &Session) {
    let i = session.selected();
    match session.stages().get(i) {
        Some(spec) => {
            ui.heading(spec.id.label());
            egui::Grid::new("params").striped(true).show(ui, |ui| {
                for (k, v) in &spec.params {
                    ui.label(k);
                    ui.label(format!("{v}"));
                    ui.end_row();
                }
            });
        }
        None => {
            ui.label("No stage selected.");
        }
    }
}

fn stack_tab(ui: &mut egui::Ui, session: &Session) {
    ui.heading("Stack at the cursor");
    // The cursor position comes from the viewport; with no viewport yet this
    // panel says so rather than showing a stack from nowhere.
    ui.label("The viewport is not open, so there is no cursor to report.");
    let _ = (stack_rows, solid_thickness_nm, session);
}

fn measure_tab(ui: &mut egui::Ui, session: &Session) {
    match session.result(session.selected()) {
        Some(res) => {
            egui::Grid::new("measure").striped(true).show(ui, |ui| {
                for row in measure_rows(res) {
                    ui.label(&row.label);
                    match row.value {
                        Some(v) => ui.label(format!("{v:.2} {}", row.unit)),
                        None => ui.colored_label(
                            egui::Color32::from_rgb(180, 120, 40),
                            format!("not measured — {}", row.reason.clone().unwrap_or_else(|| "no reason".into())),
                        ),
                    };
                    // G-03: every number is one hover from its definition
                    ui.label("?").on_hover_text(row.provenance.summary());
                    ui.end_row();
                }
            });
        }
        None => {
            ui.label("Nothing to show: this stage has not been committed.");
        }
    }
}

fn device_tab(ui: &mut egui::Ui) {
    ui.heading("Device");
    ui.label("No cross-section is attached to this session yet.");
    let _ = (cv_series, iv_series);
}

fn evidence_tab(ui: &mut egui::Ui, session: &Session) {
    ui.heading("Evidence");
    ui.label(format!("state {:016x}", session.fingerprint()));
    ui.label("A number without a model, a parameter list and an uncertainty is not evidence.");
}

/// Draw a series with the core painter, so no plotting crate is needed (G-14).
pub fn draw_series(ui: &mut egui::Ui, series: &Series, height: f32) {
    if series.is_empty() {
        ui.label(format!("{}: no data", series.title));
        return;
    }
    let (response, painter) = ui.allocate_painter(egui::vec2(ui.available_width(), height), egui::Sense::hover());
    let rect = response.rect;
    let (x0, x1) = series.x_range().unwrap_or((0.0, 1.0));
    let (y0, y1) = series.y_range().unwrap_or((0.0, 1.0));
    let sx = |x: f64| rect.left() + ((x - x0) / (x1 - x0).max(f64::MIN_POSITIVE)) as f32 * rect.width();
    let sy = |y: f64| rect.bottom() - ((y - y0) / (y1 - y0).max(f64::MIN_POSITIVE)) as f32 * rect.height();
    for w in series.points.windows(2) {
        painter.line_segment(
            [egui::pos2(sx(w[0].0), sy(w[0].1)), egui::pos2(sx(w[1].0), sy(w[1].1))],
            egui::Stroke::new(1.5_f32, egui::Color32::LIGHT_BLUE),
        );
    }
    painter.text(rect.left_bottom(), egui::Align2::LEFT_TOP, &series.x_label, egui::FontId::default(), egui::Color32::GRAY);
    painter.text(rect.left_top(), egui::Align2::LEFT_BOTTOM, &series.y_label, egui::FontId::default(), egui::Color32::GRAY);
}

/// Always-visible state: fingerprint, history depth and any named progress (G-23).
pub fn status_bar(ui: &mut egui::Ui, session: &Session, progress: Option<&Progress>) {
    ui.horizontal(|ui| {
        ui.label(format!("state {:016x}", session.fingerprint()));
        ui.separator();
        ui.label(format!("history {}/{}", session.history_depth(), session.history_limit()));
        if let Some(p) = progress {
            ui.separator();
            ui.add(egui::ProgressBar::new(p.fraction() as f32).text(p.text()));
        }
        ui.separator();
        // G-19: the palette is the same table the shortcuts come from
        ui.menu_button("Commands", |ui| {
            for entry in palette() {
                let label = if entry.binding.is_empty() {
                    entry.label
                } else {
                    format!("{}  ({})", entry.label, entry.binding)
                };
                ui.label(label);
            }
        });
    });
}
