//! Headless render of the real shader (software Vulkan works): `cargo run --release -p sim_ui --example render_views -- out_dir`
use eframe::egui::{self, Pos2, RawInput, Rect, Vec2};
use eframe::egui_wgpu::wgpu;
use sim_ui::app::AttoFabApp;
use sim_ui::svo_view::render_offscreen;
use sim_ui::viewport_math::ViewMode;

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "render_out".into());
    std::fs::create_dir_all(&out).unwrap();
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor { backends: wgpu::Backends::VULKAN, ..Default::default() });
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions { force_fallback_adapter: false, ..Default::default() })).expect("no Vulkan adapter");
    println!("adapter: {}", adapter.get_info().name);
    let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default(), None)).expect("device");
    let ctx = egui::Context::default();
    let mut app = AttoFabApp::new_headless(&ctx);
    let mut t = 0.0;
    let mut frame = |app: &mut AttoFabApp| { t += 1.0 / 60.0; let _ = ctx.run(RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1500.0, 950.0))), time: Some(t), predicted_dt: 1.0 / 60.0, ..Default::default() }, |c| app.ui(c)); };
    for _ in 0..3 { frame(&mut app); }
    app.do_run();
    for _ in 0..4 { frame(&mut app); }
    app.stop_and_show_final();
    for _ in 0..3 { frame(&mut app); }
    let (w, h) = (960u32, 600u32);
    let views: [(&str, ViewMode, &str); 7] = [
        ("3d_site", ViewMode::Perspective, "site"), ("top_site", ViewMode::Top, "site"), ("xz_cut", ViewMode::Xz, "site"),
        ("yz_cut", ViewMode::Yz, "site"), ("slice", ViewMode::Slice, "site"), ("3d_wafer", ViewMode::Perspective, "wafer"), ("top_wafer", ViewMode::Top, "wafer"),
    ];
    for (name, mode, focus) in views {
        app.set_aspect(w as f64 / h as f64);
        app.frame_view(mode, focus == "site");
        frame(&mut app);
        let (nodes, u) = app.offscreen_inputs(w as f32, h as f32);
        let px = render_offscreen(&device, &queue, &nodes, u, w, h, [0.043, 0.063, 0.094, 1.0]);
        let img = image::RgbaImage::from_raw(w, h, px).unwrap();
        let nonbg = img.pixels().filter(|p| (p[0] as i32 - 11).abs() > 3 || (p[1] as i32 - 16).abs() > 3).count();
        println!("{name}: {} non-background pixels", nonbg);
        img.save(format!("{out}/{name}.png")).unwrap();
    }
}
