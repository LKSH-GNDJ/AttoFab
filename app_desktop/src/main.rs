// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
fn main() -> eframe::Result<()> {
    let (rgba, width, height) = sim_ui::app_icon_rgba();
    let options = eframe::NativeOptions {
        renderer: eframe::Renderer::Wgpu,
        viewport: eframe::egui::ViewportBuilder::default()
            .with_inner_size([1440.0, 900.0])
            .with_min_inner_size([980.0, 640.0])
            .with_title("AttoFab — Semiconductor Fabrication Simulator")
            .with_icon(std::sync::Arc::new(eframe::egui::IconData { rgba, width, height })),
        ..Default::default()
    };

    eframe::run_native(
        "AttoFab",
        options,
        Box::new(|cc| Box::new(sim_ui::app::AttoFabApp::new(cc))),
    )
}
