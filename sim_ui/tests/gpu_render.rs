// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Renders the real WGSL shader offscreen and checks what it draws.
//! Skips (passes with a note) when no Vulkan/DX12/Metal adapter exists. CI installs lavapipe.

use eframe::egui::{self, Pos2, RawInput, Rect, Vec2};
use eframe::egui_wgpu::wgpu;
use sim_core::legend::{Palette, Theme};
use sim_ui::app::AttoFabApp;
use sim_ui::svo_view::{render_offscreen, render_offscreen_in};
use sim_ui::viewport_math::ViewMode;

const W: u32 = 640;
const H: u32 = 400;

fn gpu() -> Option<(wgpu::Device, wgpu::Queue)> {
    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::default());
    let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions::default()))?;
    pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default(), None)).ok()
}

struct Scene { ctx: egui::Context, app: AttoFabApp, t: f64 }
impl Scene {
    fn new() -> Self {
        let ctx = egui::Context::default();
        let app = AttoFabApp::new_headless(&ctx);
        let mut s = Self { ctx, app, t: 0.0 };
        for _ in 0..3 { s.frame(); }
        s.app.do_run();
        for _ in 0..4 { s.frame(); }
        s.app.stop_and_show_final();
        for _ in 0..3 { s.frame(); }
        s
    }
    fn frame(&mut self) {
        self.t += 1.0 / 60.0;
        let app = &mut self.app;
        let _ = self.ctx.run(RawInput { screen_rect: Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1500.0, 950.0))), time: Some(self.t), predicted_dt: 1.0 / 60.0, ..Default::default() }, |c| app.ui(c));
    }
    fn render(&mut self, gpu: &(wgpu::Device, wgpu::Queue), mode: ViewMode, site: bool) -> Vec<[u8; 3]> {
        self.app.set_aspect(W as f64 / H as f64);
        self.app.frame_view(mode, site);
        self.frame();
        let (nodes, u) = self.app.offscreen_inputs(W as f32, H as f32);
        let px = render_offscreen(&gpu.0, &gpu.1, &nodes, u, W, H, [11.0 / 255.0, 16.0 / 255.0, 24.0 / 255.0, 1.0]);
        px.chunks_exact(4).map(|p| [p[0], p[1], p[2]]).collect()
    }
}

fn near(p: [u8; 3], c: [u8; 3], tol: f64) -> bool { (0..3).map(|i| (p[i] as f64 - c[i] as f64).powi(2)).sum::<f64>().sqrt() < tol }
fn count(img: &[[u8; 3]], c: [u8; 3]) -> usize { img.iter().filter(|p| near(**p, c, 30.0)).count() }

#[test]
fn shader_draws_tile_film_dopants_and_bulk_with_correct_geometry() {
    let Some(gpu) = gpu() else { eprintln!("no GPU adapter: skipping"); return };
    let pal = Palette::new(Theme::Dark);
    let (film, sub, bg) = (pal.sio2, pal.layer(sim_core::legend::LegendLayer::Substrate), pal.background);
    let dopant: Vec<[u8; 3]> = (1..=6).map(|b| pal.dopant_rgb(b)).collect();
    let mut s = Scene::new();

    // 3D at the site: oxide lines, implanted trenches and silicon bulk; nothing left as background.
    let img = s.render(&gpu, ViewMode::Perspective, true);
    assert!(count(&img, film) > 2000, "SiO2 film pixels: {}", count(&img, film));
    assert!(dopant.iter().map(|d| count(&img, *d)).sum::<usize>() > 2000, "dopant pixels");
    assert!(count(&img, sub) > 10_000, "silicon bulk fills the surroundings");
    assert!(count(&img, bg) < img.len() / 50, "the wafer fills the frame when zoomed onto the tile");

    // Whole wafer from above: a filled disk with background in the corners.
    let img = s.render(&gpu, ViewMode::Top, false);
    let at = |x: u32, y: u32| img[(y * W + x) as usize];
    assert!(near(at(W / 2, H / 2), sub, 30.0) && near(at(4, 4), bg, 8.0) && near(at(W - 5, H - 5), bg, 8.0), "disk with empty corners");
    let disk = count(&img, sub) as f64 / (W * H) as f64;
    assert!(disk > 0.3 && disk < 0.7, "disk fraction {disk}");

    // XZ cross-section: air above, film strictly above the original surface, dopant below it.
    let img = s.render(&gpu, ViewMode::Xz, true);
    let surface_row = (0..H).find(|y| !near(at_img(&img, 2, *y), bg, 8.0)).expect("bulk reaches the left edge");
    let film_rows: Vec<u32> = (0..H).filter(|y| (0..W).any(|x| near(at_img(&img, x, *y), film, 25.0))).collect();
    assert!(!film_rows.is_empty() && film_rows.iter().all(|r| *r < surface_row + 2), "film must sit above the substrate surface (row {surface_row}): {:?}", (film_rows.first(), film_rows.last()));
    let dop_below = (surface_row + 2..H).any(|y| (0..W).any(|x| dopant.iter().any(|d| near(at_img(&img, x, y), *d, 25.0))));
    assert!(dop_below, "implanted silicon appears below the surface in the cut face");
    assert!(count(&img, bg) > img.len() / 5, "air above the surface is visible in the section");
}

fn at_img(img: &[[u8; 3]], x: u32, y: u32) -> [u8; 3] { img[(y * W + x) as usize] }

#[test]
fn hiding_a_layer_removes_it_from_the_render() {
    let Some(gpu) = gpu() else { eprintln!("no GPU adapter: skipping"); return };
    let pal = Palette::new(Theme::Dark);
    let mut s = Scene::new();
    let before = count(&s.render(&gpu, ViewMode::Perspective, true), pal.sio2);
    s.app.set_layer_visible(sim_core::legend::LegendLayer::DepositedFilm, false);
    let after = count(&s.render(&gpu, ViewMode::Perspective, true), pal.sio2);
    assert!(before > 2000 && after < before / 20, "film pixels {before} -> {after}");
}

/// Regression for "rendered geometry is displaced from the overlays": the live viewport is inset in the
/// window, so the shader must subtract the viewport origin. Rendering through an inset viewport must give
/// exactly the same picture as rendering at the origin.
#[test]
fn an_inset_viewport_renders_the_same_picture_as_one_at_the_origin() {
    let Some(gpu) = gpu() else { eprintln!("no GPU adapter: skipping"); return };
    let mut s = Scene::new();
    s.app.set_aspect(W as f64 / H as f64);
    s.app.frame_view(ViewMode::Perspective, true);
    s.frame();
    let (nodes, u) = s.app.offscreen_inputs(W as f32, H as f32);
    let bg = [11.0 / 255.0, 16.0 / 255.0, 24.0 / 255.0, 1.0];
    let base = render_offscreen(&gpu.0, &gpu.1, &nodes, u, W, H, bg);
    let (ox, oy, tw, th) = (288u32, 59u32, W + 600, H + 250);
    let big = render_offscreen_in(&gpu.0, &gpu.1, &nodes, u, (tw, th), (ox, oy, W, H), bg);
    let mut diff = 0usize;
    for y in 0..H { for x in 0..W {
        let (a, b) = (((y * W + x) * 4) as usize, (((y + oy) * tw + x + ox) * 4) as usize);
        if (0..3).any(|k| (base[a + k] as i32 - big[b + k] as i32).abs() > 2) { diff += 1; }
    } }
    assert!(diff < (W * H / 500) as usize, "{diff} of {} pixels differ: the shader is not viewport-relative", W * H);
}

/// Select tool: clicking geometry must produce a visible highlight in the real render.
#[test]
fn selecting_a_leaf_highlights_it_in_the_render() {
    let Some(gpu) = gpu() else { eprintln!("no GPU adapter: skipping"); return };
    let mut s = Scene::new();
    s.render(&gpu, ViewMode::Perspective, true);
    // Probe a few points near the middle of the frame until one hits geometry.
    let mut hit = false;
    for (x, y) in [(0.0, 0.0), (0.1, -0.1), (-0.1, 0.05), (0.2, 0.0), (0.0, -0.2)] { if s.app.select_at_ndc([x, y]) { hit = true; break; } }
    assert!(hit, "the Select tool found no geometry near the centre of the tile view");
    let img = s.render(&gpu, ViewMode::Perspective, true);
    let gold = img.iter().filter(|p| near(**p, [255, 209, 46], 70.0)).count();
    // One 16 nm leaf is only a few pixels at tile zoom; the viewport adds a ring marker (see the UI test) on top.
    assert!(gold >= 3, "the selected leaf must be drawn in the highlight colour (gold pixels: {gold})");
}

/// The symptom reported from the live app: rendered geometry sat far from the CPU-drawn overlays
/// (process window, origin, metrology sites). Render through an inset viewport exactly like the window does and
/// require the film pixels to fall inside the projected process window.
#[test]
fn rendered_geometry_lands_inside_the_projected_process_window() {
    let Some(gpu) = gpu() else { eprintln!("no GPU adapter: skipping"); return };
    let pal = Palette::new(Theme::Dark);
    let mut s = Scene::new();
    s.app.set_aspect(W as f64 / H as f64);
    s.app.frame_view(ViewMode::Top, true);
    s.frame();
    let (nodes, u) = s.app.offscreen_inputs(W as f32, H as f32);
    let (ox, oy, tw, th) = (288u32, 59u32, W + 700, H + 300);
    let big = render_offscreen_in(&gpu.0, &gpu.1, &nodes, u, (tw, th), (ox, oy, W, H), [11.0 / 255.0, 16.0 / 255.0, 24.0 / 255.0, 1.0]);
    let (mut lo, mut hi, mut n) = ([f32::INFINITY; 2], [f32::NEG_INFINITY; 2], 0usize);
    for y in 0..th { for x in 0..tw {
        let i = ((y * tw + x) * 4) as usize;
        if near([big[i], big[i + 1], big[i + 2]], pal.sio2, 25.0) {
            n += 1;
            lo = [lo[0].min(x as f32), lo[1].min(y as f32)]; hi = [hi[0].max(x as f32), hi[1].max(y as f32)];
        }
    } }
    assert!(n > 500, "film pixels found: {n}");
    // Overlay corners of the process window at the surface (viewport-local px -> window px by adding the inset).
    let w = s.app.site_window();
    let z = s.app.surface_z_nm();
    let (hx, hy) = (w.size[0] * 0.5, w.size[1] * 0.5);
    let corners: Vec<[f32; 2]> = [[-hx, -hy], [hx, -hy], [hx, hy], [-hx, hy]].iter()
        .map(|d| s.app.project_local_px([w.center[0] + d[0], w.center[1] + d[1], z], W as f32, H as f32).expect("projects")).map(|p| [p[0] + ox as f32, p[1] + oy as f32]).collect();
    let (wlo, whi) = (corners.iter().fold([f32::INFINITY; 2], |a, c| [a[0].min(c[0]), a[1].min(c[1])]), corners.iter().fold([f32::NEG_INFINITY; 2], |a, c| [a[0].max(c[0]), a[1].max(c[1])]));
    let tol = 6.0;
    assert!(lo[0] >= wlo[0] - tol && lo[1] >= wlo[1] - tol && hi[0] <= whi[0] + tol && hi[1] <= whi[1] + tol,
        "film pixels {lo:?}..{hi:?} must lie inside the overlay window {wlo:?}..{whi:?} (the old bug put them ~{} px away)", ox);
    // And the film is centred on the window centre, not off to one side.
    let (cx, cy) = ((lo[0] + hi[0]) * 0.5, (lo[1] + hi[1]) * 0.5);
    let (wcx, wcy) = ((wlo[0] + whi[0]) * 0.5, (wlo[1] + whi[1]) * 0.5);
    assert!((cx - wcx).abs() < 25.0 && (cy - wcy).abs() < 25.0, "film centre ({cx:.0},{cy:.0}) vs window centre ({wcx:.0},{wcy:.0})");
}
