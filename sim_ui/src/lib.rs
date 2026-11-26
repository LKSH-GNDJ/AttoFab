// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! `sim_ui`: AttoFab's presentation layer.
//!
//! This crate is UI-toolkit-facing but toolkit-agnostic at the type level:
//! each panel module exposes a plain-data "view model" derived from
//! `sim_core` state, plus the interaction/query logic a real immediate-mode
//! UI (egui) would drive every frame. The actual `egui`/`egui_plot` widget
//! calls live in the platform shells' render loop and are intentionally
//! thin — they just walk these view models and paint them.
//!
//! (This scaffold omits the `egui`/`egui_plot`/`wgpu` crate dependencies
//! themselves to keep the reference build fast and network-light; see
//! `Cargo.toml` for how to re-enable them for a full graphical build.)

pub mod app;
pub mod shell;
mod app_gds;
mod app_inspector;
mod app_viewport;
pub mod gds_view;
pub mod svo_view;
pub mod gds_canvas;
pub mod plot_view;
pub mod process_builder;
pub mod table_view;
pub mod viewport_math;

/// The application icon (256×256 RGBA) for window / taskbar use.
pub fn app_icon_rgba() -> (Vec<u8>, u32, u32) {
    let img = image::load_from_memory(include_bytes!("../assets/attofab_app_icon.png")).expect("bundled icon decodes").to_rgba8();
    let (w, h) = (img.width(), img.height());
    (img.into_raw(), w, h)
}

#[cfg(test)]
mod icon_tests {
    use super::app_icon_rgba;

    #[test]
    fn bundled_icon_is_a_256_square_with_rounded_transparent_corners_and_a_bright_mark() {
        let (px, w, h) = app_icon_rgba();
        assert_eq!((w, h), (256, 256));
        let at = |x: u32, y: u32| { let i = ((y * w + x) * 4) as usize; [px[i], px[i + 1], px[i + 2], px[i + 3]] };
        assert_eq!(at(0, 0)[3], 0, "rounded corner is transparent");
        assert_eq!(at(255, 255)[3], 0);
        let bg = at(128, 8);
        assert!(bg[3] == 255 && bg[0] < 60 && bg[1] < 70 && bg[2] < 90, "dark tile: {bg:?}");
        let bright = (0..h).flat_map(|y| (0..w).map(move |x| (x, y))).filter(|&(x, y)| { let p = at(x, y); p[3] == 255 && p[0] > 220 && p[1] > 220 }).count();
        assert!(bright > 6000, "the white logo mark is present ({bright} px)");
    }
}
