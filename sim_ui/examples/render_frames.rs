// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Contact sheet of animation keyframes for every stage: `cargo run --release -p sim_ui --example render_frames -- out.png`
use eframe::egui::{self, Pos2, RawInput, Rect, Vec2};
use eframe::egui_wgpu::wgpu;
use sim_ui::app::AttoFabApp;
use sim_ui::svo_view::render_offscreen;
use sim_ui::viewport_math::ViewMode;

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "frames.png".into());
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor { backends: wgpu::Backends::VULKAN, ..Default::default() });
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions::default())).expect("adapter");
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default(), None)).expect("device");
    let ctx = egui::Context::default();
    let mut app = AttoFabApp::new_headless(&ctx);
    let mut t = 0.0;
    let mut frame = |app: &mut AttoFabApp| { t += 1.0 / 60.0; let _ = ctx.run(RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1500.0, 950.0))), time: Some(t), predicted_dt: 1.0 / 60.0, ..Default::default() }, |c| app.ui(c)); };
    for _ in 0..3 { frame(&mut app); }
    app.do_run();
    for _ in 0..4 { frame(&mut app); }
    let (w, h) = (320u32, 200u32);
    // Cross-section view shows growth/removal most clearly.
    app.set_aspect(w as f64 / h as f64);
    app.frame_view(ViewMode::Xz, true);
    let fractions = [0.0, 0.2, 0.45, 0.7, 0.85, 1.0];
    let segs = [1usize, 2, 3, 4, 5];   // deposition, litho, etch, doping, clean
    let mut sheet = image::RgbaImage::new(w * fractions.len() as u32, h * segs.len() as u32);
    for (r, seg) in segs.iter().enumerate() {
        for (c, f) in fractions.iter().enumerate() {
            app.show_position(*seg, *f);
            frame(&mut app);
            let (nodes, u) = app.offscreen_inputs(w as f32, h as f32);
            let px = render_offscreen(&device, &queue, &nodes, u, w, h, [0.043, 0.063, 0.094, 1.0]);
            let img = image::RgbaImage::from_raw(w, h, px).unwrap();
            image::imageops::replace(&mut sheet, &img, (c as u32 * w) as i64, (r as u32 * h) as i64);
        }
    }
    sheet.save(&out).unwrap();
    println!("wrote {out}");
}
