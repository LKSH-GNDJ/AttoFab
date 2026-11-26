// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Right-hand process inspector: active stage, parameters, validation, results,
//! metrology report, material legend (with visibility + colour configuration),
//! wafer transform and the inspect readout.

use eframe::egui::{self, Color32, RichText, Sense, Vec2};
use sim_core::geom::{EtchChemistry, MaskPattern};
use sim_core::legend::{layers_for_stage, LegendLayer};
use sim_core::measure::{MetrologyReport, Status};
use sim_core::params::*;
use sim_core::physics::Stage;
use sim_core::pipeline::{Fidelity, StageOutcome};
use sim_core::session::StageStatus;
use sim_core::transform::WaferTransform;
use sim_core::wafer::Material;

use crate::app::{col, status_color, status_word, AttoFabApp};

fn dv(ui: &mut egui::Ui, label: &str, v: &mut f64, range: std::ops::RangeInclusive<f64>, speed: f64, suffix: &str) -> bool {
    ui.horizontal(|ui| { ui.label(label); ui.add(egui::DragValue::new(v).speed(speed).clamp_range(range).suffix(suffix)).changed() }).inner
}

fn slider(ui: &mut egui::Ui, label: &str, v: &mut f64, range: std::ops::RangeInclusive<f64>, log: bool) -> bool {
    ui.add(egui::Slider::new(v, range).logarithmic(log).text(label)).changed()
}

const FILMS: [Material; 5] = [Material::SiliconDioxide, Material::SiliconCarbide, Material::GalliumNitride, Material::Metal, Material::Silicon];

fn mat_combo(ui: &mut egui::Ui, id: &str, m: &mut Material) -> bool {
    let before = *m;
    egui::ComboBox::from_id_source(id).selected_text(mat_name(*m)).show_ui(ui, |ui| { for f in FILMS { ui.selectable_value(m, f, mat_name(f)); } });
    before != *m
}

fn edit_params(ui: &mut egui::Ui, p: &mut StageParams) -> bool {
    let mut c = false;
    match p {
        StageParams::SubstratePrep(x) => {
            c |= dv(ui, "RCA dilution", &mut x.rca_dilution, 0.0..=1000.0, 0.1, "");
            c |= dv(ui, "Roughness Ra", &mut x.ra_nm, 0.0..=1000.0, 0.05, " nm");
            c |= dv(ui, "γ solid-gas", &mut x.gamma_sg, 0.0..=1000.0, 0.5, " mN/m");
            c |= dv(ui, "γ solid-liquid", &mut x.gamma_sl, 0.0..=1000.0, 0.5, " mN/m");
            c |= dv(ui, "γ liquid-gas", &mut x.gamma_lg, 0.001..=1000.0, 0.5, " mN/m");
            c |= ui.checkbox(&mut x.hf_dip, "HF dip (strip native oxide)").changed();
            c |= slider(ui, "particle removal", &mut x.particle_removal_eff, 0.0..=1.0, false);
        }
        StageParams::Deposition(x) => {
            ui.horizontal(|ui| { ui.label("Material"); c |= mat_combo(ui, "dep_mat", &mut x.material); });
            c |= slider(ui, "thickness (nm)", &mut x.thickness_nm, 1.0..=600.0, true);
            c |= slider(ui, "step coverage", &mut x.step_coverage, 0.0..=1.0, false);
            ui.label(RichText::new("0 = directional (PVD), 1 = fully conformal (CVD/ALD). Film grows above the surface.").weak().small());
            c |= dv(ui, "Wafer bow radius", &mut x.bow_radius_m, -1.0e9..=1.0e9, 1.0, " m");
            c |= ui.checkbox(&mut x.flux3d, "3-D transport solve")
                .on_hover_text("Solve the arriving flux per column from the geometry, instead of one step-coverage number for the whole grid.")
                .changed();
            if x.flux3d {
                c |= slider(ui, "sticking", &mut x.sticking, 0.0..=1.0, false);
                ui.label(RichText::new("What the surface does with what arrives. Separate from step coverage, which is what arrives.").weak().small());
                let source_before = x.source;
                egui::ComboBox::from_id_source("dep_source")
                    .selected_text(match x.source {
                        DepositionSource::Cosine => "Cosine (CVD)",
                        DepositionSource::Beam { .. } => "Beam (PVD)",
                    })
                    .show_ui(ui, |ui| {
                        if ui.selectable_label(matches!(x.source, DepositionSource::Cosine), "Cosine (CVD)").clicked() {
                            x.source = DepositionSource::Cosine;
                        }
                        if ui.selectable_label(matches!(x.source, DepositionSource::Beam { .. }), "Beam (PVD)").clicked() {
                            x.source = DepositionSource::Beam { tilt_deg: 0.0, spread_deg: 10.0 };
                        }
                    });
                if source_before != x.source { c = true; }
                if let DepositionSource::Beam { tilt_deg, spread_deg } = &mut x.source {
                    c |= dv(ui, "Beam tilt", tilt_deg, -89.0..=89.0, 1.0, " °");
                    c |= dv(ui, "Beam spread", spread_deg, 0.0..=90.0, 1.0, " °");
                }
                c |= ui.add(egui::Slider::new(&mut x.flux_iterations, 1..=64).text("flux iterations")).changed();
                c |= slider(ui, "flux tolerance", &mut x.flux_tolerance, 0.0..=0.5, false);
                ui.label(RichText::new("Stop early once the mean flux stops changing by this fraction. 0 runs every iteration.").weak().small());
                ui.label(RichText::new("More iterations re-solve as the film grows, so a narrowing mouth is accounted for.").weak().small());
            }
        }
        StageParams::Lithography(x) => {
            c |= dv(ui, "Resist thickness", &mut x.resist_thickness_nm, 1.0..=1.0e5, 1.0, " nm");
            c |= edit_mask(ui, &mut x.mask);
            c |= dv(ui, "Wavelength", &mut x.wavelength_nm, 1.0..=1.0e4, 1.0, " nm");
            c |= dv(ui, "NA", &mut x.na, 0.01..=2.0, 0.01, "");
            c |= dv(ui, "k1", &mut x.k1, 0.0..=2.0, 0.01, "");
            c |= slider(ui, "dose / dose-to-clear", &mut x.dose_ratio, 0.2..=10.0, true);
            c |= ui.checkbox(&mut x.positive_tone, "Positive-tone resist").changed();
            c |= dv(ui, "Alignment offset X", &mut x.align_offset_nm[0], -1.0e5..=1.0e5, 1.0, " nm");
            c |= dv(ui, "Alignment offset Y", &mut x.align_offset_nm[1], -1.0e5..=1.0e5, 1.0, " nm");
            c |= dv(ui, "Soft bake", &mut x.soft_bake_c, 20.0..=300.0, 1.0, " °C");
            c |= dv(ui, "Post-exposure bake", &mut x.peb_c, 20.0..=300.0, 1.0, " °C");
            c |= dv(ui, "Develop time", &mut x.develop_time_min, 0.01..=60.0, 0.05, " min");
            c |= dv(ui, "Planarity error", &mut x.planarity_error_nm, 0.0..=1000.0, 1.0, " nm");
            ui.label(RichText::new("Field non-planarity from an earlier polish; prints as defocus.").weak().small());
        }
        StageParams::Etching(x) => {
            let before = x.chemistry;
            egui::ComboBox::from_id_source("etch_chem").selected_text(x.chemistry.label()).show_ui(ui, |ui| { for k in EtchChemistry::ALL { ui.selectable_value(&mut x.chemistry, k, k.label()); } });
            if before != x.chemistry { if x.chemistry != EtchChemistry::Custom { x.rates = x.chemistry.default_rates(); } c = true; }
            c |= dv(ui, "Etch time", &mut x.time_min, 0.0..=1.0e4, 0.05, " min");
            c |= slider(ui, "anisotropy", &mut x.anisotropy, 0.0..=1.0, false);
            egui::CollapsingHeader::new("Etch rates (nm/min)").show(ui, |ui| {
                let mut r = false;
                r |= dv(ui, "Si / dopant", &mut x.rates.si, 0.0..=1.0e5, 0.5, ""); r |= dv(ui, "SiO₂", &mut x.rates.sio2, 0.0..=1.0e5, 0.5, "");
                r |= dv(ui, "SiC", &mut x.rates.sic, 0.0..=1.0e5, 0.5, ""); r |= dv(ui, "GaN", &mut x.rates.gan, 0.0..=1.0e5, 0.5, "");
                r |= dv(ui, "Photoresist", &mut x.rates.resist, 0.0..=1.0e5, 0.5, ""); r |= dv(ui, "Metal", &mut x.rates.metal, 0.0..=1.0e5, 0.5, "");
                if r { x.chemistry = EtchChemistry::Custom; c = true; }
            });
            c |= ui.checkbox(&mut x.strip_after, "Strip resist at the end of this stage").changed();
        }
        StageParams::Doping(x) => {
            egui::ComboBox::from_id_source("species").selected_text(x.species.label()).show_ui(ui, |ui| { for s in Species::ALL { c |= ui.selectable_value(&mut x.species, s, s.label()).changed(); } });
            c |= slider(ui, "energy (keV)", &mut x.energy_kev, 1.0..=400.0, true);
            c |= slider(ui, "dose (cm⁻²)", &mut x.dose_cm2, 1.0e11..=1.0e16, true);
            c |= dv(ui, "Anneal temperature", &mut x.anneal_temp_c, 25.0..=1400.0, 5.0, " °C");
            c |= dv(ui, "Anneal time", &mut x.anneal_min, 0.0..=1.0e4, 0.5, " min");
            c |= slider(ui, "background (cm⁻³)", &mut x.background_cm3, 1.0e12..=1.0e19, true);
            let mut ov = x.rp_override_nm.is_some();
            if ui.checkbox(&mut ov, "Override projected range").changed() { x.rp_override_nm = if ov { Some(100.0) } else { None }; c = true; }
            if let Some(r) = &mut x.rp_override_nm { c |= dv(ui, "Rp", r, 0.1..=1.0e5, 1.0, " nm"); }
        }
        StageParams::FinalClean(x) => {
            c |= slider(ui, "strip efficiency", &mut x.strip_efficiency, 0.0..=1.0, false);
            c |= ui.checkbox(&mut x.chemical_clean, "Chemical clean (RCA) + rinse + dry").changed();
            c |= slider(ui, "particle removal", &mut x.particle_removal_eff, 0.0..=1.0, false);
        }
        StageParams::Cmp(x) => {
            ui.horizontal(|ui| { ui.label("Polishing"); c |= mat_combo(ui, "cmp_mat", &mut x.material); });
            c |= dv(ui, "Pressure", &mut x.pressure_kpa, 0.1..=200.0, 0.5, " kPa");
            c |= dv(ui, "Pad velocity", &mut x.velocity_mps, 0.01..=5.0, 0.01, " m/s");
            c |= dv(ui, "Polish time", &mut x.time_min, 0.01..=60.0, 0.05, " min");
            c |= slider(ui, "over-polish", &mut x.over_polish, 0.0..=1.0, false);
            ui.label(RichText::new("Removal is Preston's k p v t, applied uniformly: dishing and erosion are not modelled here.").weak().small());
        }
        StageParams::Metrology(x) => {
            let mut sites = x.sites as f64;
            if slider(ui, "wafer-map sites per side", &mut sites, 1.0..=15.0, false) { x.sites = sites.round() as usize; c = true; }
            c |= dv(ui, "Repeatability noise (seeded)", &mut x.noise_nm, 0.0..=1000.0, 0.05, " nm");
            c |= dv(ui, "Ellipsometer λ", &mut x.wavelength_nm, 100.0..=2000.0, 1.0, " nm");
            c |= dv(ui, "Incidence angle", &mut x.incidence_deg, 1.0..=89.0, 0.5, "°");
            ui.label(RichText::new("Targets").strong());
            let mut remove = None;
            for (i, t) in x.targets.iter_mut().enumerate() {
                ui.horizontal(|ui| {
                    ui.label(t.quantity.label());
                    c |= ui.add(egui::DragValue::new(&mut t.nominal).speed(0.5)).changed();
                    ui.label("±");
                    c |= ui.add(egui::DragValue::new(&mut t.tolerance).speed(0.5).clamp_range(0.0..=1.0e9)).changed();
                    ui.label(t.quantity.unit());
                    if ui.small_button("✖").clicked() { remove = Some(i); }
                });
            }
            if let Some(i) = remove { x.targets.remove(i); c = true; }
            ui.menu_button("＋ Add target", |ui| {
                let presets: [(Quantity, f64, f64); 10] = [
                    (Quantity::FilmThickness(Material::SiliconDioxide), 100.0, 10.0), (Quantity::FilmThickness(Material::Metal), 100.0, 10.0), (Quantity::FilmThickness(Material::SiliconCarbide), 100.0, 10.0),
                    (Quantity::LineWidth(Material::SiliconDioxide), 150.0, 40.0), (Quantity::LineWidth(Material::Photoresist), 150.0, 40.0), (Quantity::StepHeight, 100.0, 20.0),
                    (Quantity::JunctionDepth, 250.0, 100.0), (Quantity::TopographyRms, 0.0, 50.0), (Quantity::ResistResidue, 0.0, 1.0), (Quantity::PeakFilmStress, 200.0, 150.0),
                ];
                for (q, n, t) in presets { if ui.button(q.label()).clicked() { x.targets.push(Target { quantity: q, nominal: n, tolerance: t }); c = true; ui.close_menu(); } }
            });
        }
    }
    c
}

fn edit_mask(ui: &mut egui::Ui, m: &mut MaskPattern) -> bool {
    let mut c = false;
    let kind = |m: &MaskPattern| match m { MaskPattern::LineSpace { vertical: true, .. } => 0, MaskPattern::LineSpace { .. } => 1, MaskPattern::Rects { .. } => 2, MaskPattern::Polygons { .. } => 3 };
    let mut k = kind(m);
    let names = ["Vertical lines", "Horizontal lines", "Single opening", "GDSII polygons"];
    ui.horizontal(|ui| {
        ui.label("Mask");
        egui::ComboBox::from_id_source("mask_kind").selected_text(names[k]).show_ui(ui, |ui| { for (i, n) in names.iter().enumerate() { ui.selectable_value(&mut k, i, *n); } });
    });
    if k != kind(m) {
        *m = match k {
            0 | 1 => { let (l, p, n, len, cf) = match m { MaskPattern::LineSpace { line_nm, pitch_nm, count, length_nm, clear_field, .. } => (*line_nm, *pitch_nm, *count, *length_nm, *clear_field), _ => (150.0, 320.0, 3, 1000.0, true) }; MaskPattern::LineSpace { line_nm: l, pitch_nm: p, count: n, length_nm: len, vertical: k == 0, clear_field: cf } }
            2 => MaskPattern::Rects { rects: vec![[-200.0, -200.0, 200.0, 200.0]], clear_field: false },
            _ => m.clone(),
        };
        c = true;
    }
    match m {
        MaskPattern::LineSpace { line_nm, pitch_nm, count, length_nm, clear_field, .. } => {
            c |= dv(ui, "Line width", line_nm, 1.0..=1.0e5, 1.0, " nm"); c |= dv(ui, "Pitch", pitch_nm, 1.0..=1.0e5, 1.0, " nm");
            let mut n = *count as f64; if dv(ui, "Lines", &mut n, 1.0..=50.0, 0.2, "") { *count = n.round().max(1.0) as u32; c = true; }
            c |= dv(ui, "Line length", length_nm, 1.0..=1.0e5, 5.0, " nm");
            c |= ui.checkbox(clear_field, "Clear-field mask (lines are opaque)").changed();
        }
        MaskPattern::Polygons { polys, clear_field } => {
            ui.label(format!("{} polygon(s) from a GDSII layer", polys.len()));
            c |= ui.checkbox(clear_field, "Clear-field mask").changed();
            ui.label(RichText::new("Load a layer with the GDS box below the stage list (native builds).").weak().small());
        }
        MaskPattern::Rects { rects, clear_field } => {
            if let Some(r) = rects.first_mut() {
                let (mut w, mut h) = (r[2] - r[0], r[3] - r[1]);
                let a = dv(ui, "Opening width", &mut w, 1.0..=1.0e5, 1.0, " nm"); let b = dv(ui, "Opening height", &mut h, 1.0..=1.0e5, 1.0, " nm");
                if a || b { *r = [-w * 0.5, -h * 0.5, w * 0.5, h * 0.5]; c = true; }
            }
            c |= ui.checkbox(clear_field, "Clear-field mask").changed();
        }
    }
    c
}

fn fmt_val(v: f64) -> String {
    if v != 0.0 && (v.abs() >= 1.0e5 || v.abs() < 1.0e-2) { format!("{v:.3e}") } else { format!("{v:.2}") }
}

impl AttoFabApp {
    pub(crate) fn inspector_panel(&mut self, ui: &mut egui::Ui) {
        if self.session.state.recipe.index_of(self.selected_stage).is_none() { self.selected_stage = self.session.state.recipe.stages.first().map(|s| s.id).unwrap_or(0); }
        let id = self.selected_stage;
        egui::ScrollArea::vertical().auto_shrink([false, false]).show(ui, |ui| {
            let Some(idx) = self.session.state.recipe.index_of(id) else { ui.label("No stages in the recipe."); return };
            let rs = self.session.state.recipe.stages[idx].clone();
            let status = self.session.status_of(id).unwrap_or(StageStatus::NotRun);
            let stage = rs.params.stage();

            ui.horizontal(|ui| { ui.heading(format!("{}. {}", idx + 1, stage.label())); ui.label(RichText::new(format!("● {}", status_word(&status))).color(status_color(&status))); });
            if let StageStatus::Failed(m) = &status { ui.colored_label(col(self.palette.layer(LegendLayer::WarningError)), format!("✖ {m}")); }
            if status == StageStatus::Stale { ui.colored_label(Color32::from_rgb(240, 180, 60), "An upstream change invalidated this result. Run Process recomputes it."); }

            egui::CollapsingHeader::new("Parameters").default_open(true).show(ui, |ui| {
                let mut p = rs.params.clone();
                if edit_params(ui, &mut p) { let _ = self.session.edit_params(id, p); }
            });

            let issues = rs.params.validate();
            if !issues.is_empty() {
                for i in &issues { ui.colored_label(if i.severity == Severity::Error { col(self.palette.layer(LegendLayer::WarningError)) } else { Color32::from_rgb(240, 180, 60) }, format!("{} {}", if i.severity == Severity::Error { "✖" } else { "⚠" }, i.message)); }
            } else { ui.label(RichText::new("✔ parameters valid").color(Color32::from_rgb(90, 200, 120)).small()); }

            // Outcome to show: the preview if it is this stage, otherwise the committed (or stale) record.
            let preview = self.session.preview.as_ref().filter(|p| p.stage_id == id).map(|p| (p.record.outcome.clone(), true));
            let rec = self.session.state.records.get(&id).map(|r| (r.outcome.clone(), false));
            if let Some((out, is_preview)) = preview.or(rec) { self.outcome_section(ui, &out, is_preview, &status); }
            else { egui::CollapsingHeader::new("Sub-steps").show(ui, |ui| { for s in sim_core::pipeline::substeps(stage) { ui.label(RichText::new(format!("• {}  [{}]", s.name, s.kind.label())).small()); } }); }

            self.verification_section(ui);
            self.legend_section(ui, stage);
            self.transform_section(ui);
            self.selection_section(ui);
            self.inspect_section(ui);

            ui.separator();
            if ui.button("📋 Copy process report (Markdown)").on_hover_text("Recipe, parameters, input/output fingerprints, fidelity, metrics and metrology").clicked() {
                let md = sim_core::report::markdown_report(&self.session);
                ui.output_mut(|o| o.copied_text = md);
                self.status = "Process report copied to the clipboard.".into();
            }
        });
    }

    fn outcome_section(&mut self, ui: &mut egui::Ui, out: &StageOutcome, is_preview: bool, status: &StageStatus) {
        ui.separator();
        ui.horizontal(|ui| {
            ui.label(RichText::new(if is_preview { "Preview result (uncommitted)" } else if *status == StageStatus::Stale { "Last result (stale)" } else { "Committed result" }).strong());
            let (txt, c) = match out.fidelity { Fidelity::Calibrated => ("calibrated", Color32::from_rgb(90, 200, 120)), Fidelity::Simulated => ("simulated", Color32::from_rgb(110, 170, 255)), Fidelity::Approximated => ("approximated", Color32::from_rgb(240, 180, 60)) };
            ui.label(RichText::new(format!("[{txt}]")).color(c)).on_hover_text("Simulated: measured from simulated geometry. Approximated: physically motivated reference model, not calibrated to a tool. Calibrated: fitted to tool data.");
        });
        ui.label(RichText::new(format!("{} · nominal process time {:.0} s", if out.geometry_changed { "geometry changed" } else { "no net geometry change" }, out.duration_s)).weak().small());
        egui::CollapsingHeader::new("Sub-steps").show(ui, |ui| { for s in &out.substeps { ui.label(RichText::new(format!("• {}  [{}]", s.name, s.kind.label())).small()); } });
        if !out.metrics.is_empty() {
            egui::Grid::new("metrics").num_columns(2).striped(true).show(ui, |ui| {
                for m in &out.metrics { ui.label(&m.name); ui.label(RichText::new(format!("{} {}", fmt_val(m.value), m.unit)).monospace()); ui.end_row(); }
            });
        }
        for i in &out.issues { ui.colored_label(if i.severity == Severity::Error { col(self.palette.layer(LegendLayer::WarningError)) } else { Color32::from_rgb(240, 180, 60) }, format!("{} {}", if i.severity == Severity::Error { "✖" } else { "⚠" }, i.message)); }
        if let Some(r) = &out.report { self.report_section(ui, r); }
    }

    fn report_section(&self, ui: &mut egui::Ui, r: &MetrologyReport) {
        ui.separator();
        let (txt, c) = match r.pass { Some(true) => ("RECIPE PASS", Color32::from_rgb(90, 200, 120)), Some(false) => ("RECIPE FAIL", col(self.palette.layer(LegendLayer::WarningError))), None => ("no measurable targets", Color32::GRAY) };
        ui.label(RichText::new(txt).strong().color(c));
        egui::Grid::new("meas").num_columns(4).striped(true).show(ui, |ui| {
            ui.label(RichText::new("measurement").small()); ui.label(RichText::new("value").small()); ui.label(RichText::new("target ± tol").small()); ui.label(RichText::new("").small()); ui.end_row();
            for m in &r.measurements {
                ui.label(&m.name);
                ui.label(RichText::new(m.value.map(fmt_val).unwrap_or("—".into())).monospace());
                ui.label(RichText::new(format!("{} ± {} {}", fmt_val(m.nominal), fmt_val(m.tolerance), m.unit)).monospace().small());
                let (s, c) = match m.status { Status::Pass => ("pass", Color32::from_rgb(90, 200, 120)), Status::Fail => ("FAIL", col(self.palette.layer(LegendLayer::WarningError))), Status::NotMeasured => ("n/m", Color32::GRAY) };
                ui.label(RichText::new(s).color(c)); ui.end_row();
            }
        });
        if let Some(o) = &r.optical { ui.label(RichText::new(format!("Ellipsometry on {} ({:.0} nm): Ψ = {:.2}°, Δ = {:.2}°", mat_name(o.film), o.film_nm, o.psi_deg, o.delta_deg)).small()); }
        ui.label(RichText::new(format!("measured wafer state: {:016x}", r.wafer_fingerprint)).monospace().small());
        // Wafer map.
        let wm = &r.wafer_map;
        ui.label(RichText::new(&wm.label).small());
        let cell = (ui.available_width().min(220.0) / wm.n as f32).clamp(10.0, 34.0);
        let (rect, resp) = ui.allocate_exact_size(Vec2::splat(cell * wm.n as f32), Sense::hover());
        let (lo, hi) = wm.values.iter().flatten().fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), v| (a.min(*v), b.max(*v)));
        for j in 0..wm.n { for i in 0..wm.n {
            let k = j * wm.n + i;
            let cr = egui::Rect::from_min_size(rect.min + Vec2::new(i as f32 * cell, (wm.n - 1 - j) as f32 * cell), Vec2::splat(cell - 1.0));
            let c = match wm.values[k] {
                None => Color32::from_gray(60),
                Some(_) if wm.out_of_spec(k) => col(self.palette.layer(LegendLayer::WarningError)),
                Some(v) => { let t = if hi > lo { ((v - lo) / (hi - lo)) as f32 } else { 0.5 }; let m = col(self.palette.layer(LegendLayer::Metrology)); m.gamma_multiply(0.35 + 0.65 * t) }
            };
            ui.painter().rect_filled(cr, 2.0, c);
            if let (Some(p), Some(v)) = (resp.hover_pos(), wm.values[k]) { if cr.contains(p) { egui::show_tooltip_at_pointer(ui.ctx(), egui::Id::new("wm_tip"), |ui| { ui.label(format!("site ({i}, {j}): {v:.1}{}", if wm.out_of_spec(k) { "  OUT OF SPEC" } else { "" })); }); } }
        } }
    }

    fn legend_section(&mut self, ui: &mut egui::Ui, stage: Stage) {
        egui::CollapsingHeader::new("Material legend & layers").default_open(true).show(ui, |ui| {
            let active = layers_for_stage(stage);
            for l in LegendLayer::ALL {
                ui.horizontal(|ui| {
                    let mut v = self.palette.visible[l.index()];
                    if ui.checkbox(&mut v, "").on_hover_text("Show / hide this layer").changed() { self.palette.visible[l.index()] = v; }
                    let mut c = self.palette.layer(l);
                    if ui.color_edit_button_srgb(&mut c).changed() { self.palette.set_layer(l, c); }
                    let t = RichText::new(l.label());
                    ui.label(if active.contains(&l) { t.strong() } else { t.weak() }).on_hover_text(l.meaning());
                    if !l.is_material() { ui.label(RichText::new("overlay").weak().small()); }
                });
            }
            ui.horizontal_wrapped(|ui| {
                ui.label(RichText::new("Films:").small());
                for (name, which) in [("SiO₂", 0), ("SiC", 1), ("GaN", 2), ("Metal", 3)] {
                    let slot = match which { 0 => &mut self.palette.sio2, 1 => &mut self.palette.sic, 2 => &mut self.palette.gan, _ => &mut self.palette.metal };
                    ui.color_edit_button_srgb(slot); ui.label(RichText::new(name).small());
                }
            });
            ui.label(RichText::new("Dopant concentration (one step per decade):").small());
            ui.horizontal(|ui| {
                for b in 1..=6u8 {
                    let (r, _) = ui.allocate_exact_size(Vec2::new(26.0, 14.0), Sense::hover());
                    ui.painter().rect_filled(r, 2.0, col(self.palette.dopant_rgb(b)));
                    ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, format!("1e{}", 16 + b as i32), egui::FontId::monospace(8.0), Color32::BLACK);
                }
            });
            if ui.small_button("Reset colours").clicked() { let vis = self.palette.visible; self.palette = sim_core::legend::Palette::new(self.palette.theme); self.palette.visible = vis; }
        });
    }

    fn verification_section(&mut self, ui: &mut egui::Ui) {
        egui::CollapsingHeader::new("Verification (V&V)").show(ui, |ui| {
            ui.label(RichText::new("Checks the process operators and FEM against exact solutions, conservation and mesh refinement. This is verification, not validation against a real fab.").weak().small());
            let busy = self.vv_rx.is_some();
            if ui.add_enabled(!busy, egui::Button::new(if busy { "Running…" } else { "Run fast verification" })).on_hover_text("Cells 32 and 16 nm; the full 3-level study with sign-off is `cargo run --release -p sim_core --bin signoff_pipeline`").clicked() { self.start_verification(); }
            if let Some(r) = &self.vv_report {
                let failed = r.failed();
                let (txt, c) = if failed.is_empty() { (format!("✔ {} cases verified", r.cases.len()), Color32::from_rgb(90, 200, 120)) } else { (format!("✖ {} of {} cases failed", failed.len(), r.cases.len()), col(self.palette.layer(LegendLayer::WarningError))) };
                ui.label(RichText::new(txt).strong().color(c));
                for f in failed { ui.colored_label(col(self.palette.layer(LegendLayer::WarningError)), format!("{}: {}", f.id, f.message)); }
                egui::CollapsingHeader::new("All cases").show(ui, |ui| {
                    for c in &r.cases { ui.label(RichText::new(format!("{} {}  {:.3} ≤ {:.3} {}", if c.status == sim_core::validation::CheckStatus::Pass { "✔" } else { "✖" }, c.id, c.value, c.limit, c.unit)).monospace().small()).on_hover_text(&c.message); }
                });
            }
        });
    }

    fn transform_section(&mut self, ui: &mut egui::Ui) {
        egui::CollapsingHeader::new("Wafer transform (environment)").show(ui, |ui| {
            let t = self.transform();
            let (mut x, mut y, mut z, mut yaw) = (t.position[0] / 1.0e6, t.position[1] / 1.0e6, t.position[2] / 1.0e6, t.yaw_rad.to_degrees());
            let mut ch = false;
            ch |= dv(ui, "X", &mut x, -1000.0..=1000.0, 0.1, " mm");
            ch |= dv(ui, "Y", &mut y, -1000.0..=1000.0, 0.1, " mm");
            ch |= dv(ui, "Z", &mut z, -1000.0..=1000.0, 0.1, " mm");
            ch |= dv(ui, "Yaw", &mut yaw, -360.0..=360.0, 0.5, "°");
            if ch && !self.locked() { self.session.set_transform(WaferTransform { position: [x * 1.0e6, y * 1.0e6, z * 1.0e6], yaw_rad: yaw.to_radians() }); }
            ui.label(RichText::new("The wafer stays circular and horizontal; only this transform changes. SVO nodes are never rewritten.").weak().small());
            if ui.button("Center Wafer").clicked() { self.session.center_wafer(); }
        });
    }

    fn selection_section(&mut self, ui: &mut egui::Ui) {
        let Some((b, m)) = self.selection else { return };
        let mut clear = false;
        egui::CollapsingHeader::new("Selection").default_open(true).show(ui, |ui| {
            ui.label(RichText::new(format!("{} leaf", match m { Material::Dopant(k) => format!("Dopant ≥1e{}", 16 + k as i32), other => mat_name(other).to_string() })).strong());
            ui.label(format!("size {:.1} × {:.1} × {:.1} nm", b.max.0 - b.min.0, b.max.1 - b.min.1, b.max.2 - b.min.2));
            let s = self.surface_z();
            ui.label(format!("centre ({:.1}, {:.1}) nm, z {:+.1} nm vs surface", (b.min.0 + b.max.0) * 0.5, (b.min.1 + b.max.1) * 0.5, (b.min.2 + b.max.2) * 0.5 - s));
            if ui.small_button("Clear selection").clicked() { clear = true; }
        });
        if clear { self.selection = None; }
    }

    fn inspect_section(&mut self, ui: &mut egui::Ui) {
        let mut clear = false;
        let Some(i) = &self.inspect else { return };
        egui::CollapsingHeader::new("Inspect").default_open(true).show(ui, |ui| {
            ui.label(format!("wafer-local ({:.1}, {:.1}) nm{}", i.x, i.y, if i.in_tile { "" } else { "  (bulk, outside the high-resolution tile)" }));
            let top = i.runs.iter().rev().find(|r| r.material != Material::Air);
            if let Some(t) = top { ui.label(format!("top surface: {} at {:+.1} nm vs original surface", mat_name(t.material), t.z1 - i.surface_z)); }
            egui::Grid::new("stack").striped(true).show(ui, |ui| {
                ui.label(RichText::new("material").small()); ui.label(RichText::new("z range (vs surface)").small()); ui.label(RichText::new("thickness").small()); ui.end_row();
                for r in i.runs.iter().rev().filter(|r| r.material != Material::Air).take(12) {
                    let label = match r.material { Material::Dopant(b) => format!("Dopant ≥1e{}", 16 + b as i32), m => mat_name(m).to_string() };
                    ui.label(label); ui.label(RichText::new(format!("{:+.0} … {:+.0} nm", r.z0 - i.surface_z, r.z1 - i.surface_z)).monospace().small()); ui.label(format!("{:.1} nm", r.z1 - r.z0)); ui.end_row();
                }
            });
            if let Some(b) = i.bucket { ui.label(RichText::new(format!("dopant concentration bucket {b} (≥ 1e{} cm⁻³)", 16 + b as i32)).small()); }
            if ui.small_button("Clear").clicked() { clear = true; }
        });
        if clear { self.inspect = None; }
    }
}
