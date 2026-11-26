//! Inset-viewport render with the CPU overlay (process window + origin) drawn on top: `-- out.png`
use eframe::egui::{self, Pos2, RawInput, Rect, Vec2};
use eframe::egui_wgpu::wgpu;
use sim_ui::app::AttoFabApp;
use sim_ui::svo_view::render_offscreen_in;
use sim_ui::viewport_math::ViewMode;

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "alignment.png".into());
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::default());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions::default())).expect("adapter");
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
    let (w, h, ox, oy, tw, th) = (900u32, 560u32, 288u32, 59u32, 1500u32, 700u32);
    let mut canvas = image::RgbaImage::from_pixel(tw, th, image::Rgba([30, 33, 40, 255]));
    for (k, (mode, name)) in [(ViewMode::Top, "Top"), (ViewMode::Perspective, "3D")].into_iter().enumerate() {
        app.set_aspect(w as f64 / h as f64);
        app.frame_view(mode, true);
        frame(&mut app);
        let (nodes, u) = app.offscreen_inputs(w as f32, h as f32);
        let px = render_offscreen_in(&device, &queue, &nodes, u, (tw, th), (ox, oy, w, h), [0.043, 0.063, 0.094, 1.0]);
        let img = image::RgbaImage::from_raw(tw, th, px).unwrap();
        let mut img = img;
        let z = app.surface_z_nm();
        let win = app.site_window();
        let (hx, hy) = (win.size[0] / 2.0, win.size[1] / 2.0);
        let pts: Vec<[f32; 2]> = [[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy], [-hx, -hy]].iter().filter_map(|d| app.project_local_px([win.center[0] + d[0], win.center[1] + d[1], z], w as f32, h as f32)).map(|p| [p[0] + ox as f32, p[1] + oy as f32]).collect();
        let mut line = |a: [f32; 2], b: [f32; 2], c: [u8; 4]| { for i in 0..=200 { let f = i as f32 / 200.0; let (x, y) = ((a[0] + (b[0] - a[0]) * f) as i64, (a[1] + (b[1] - a[1]) * f) as i64); for dx in 0..2 { for dy in 0..2 { let (px, py) = (x + dx, y + dy); if px >= 0 && py >= 0 && (px as u32) < tw && (py as u32) < th { img.put_pixel(px as u32, py as u32, image::Rgba(c)); } } } } };
        for w2 in pts.windows(2) { line(w2[0], w2[1], [255, 60, 60, 255]); }
        if let Some(o) = app.project_local_px([0.0, 0.0, z], w as f32, h as f32) { let o = [o[0] + ox as f32, o[1] + oy as f32]; line([o[0] - 10.0, o[1]], [o[0] + 10.0, o[1]], [255, 255, 255, 255]); line([o[0], o[1] - 10.0], [o[0], o[1] + 10.0], [255, 255, 255, 255]); }
        let _ = (k, name);
        if k == 0 { canvas = img; } else { image::imageops::replace(&mut canvas, &image::imageops::crop_imm(&img, ox, oy, w, h).to_image(), 0, 0); }
    }
    canvas.save(&out).unwrap();
    println!("wrote {out}");
}
