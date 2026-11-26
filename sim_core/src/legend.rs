//! One material palette for the viewport shader, the egui panels and reports
//! (PRD §4). Colours are configurable, distinguishable and contrast-checked
//! against the viewport background in both themes.

use crate::physics::Stage;
use crate::render_wgpu::{DOPANT_ID_BASE, PALETTE_SLOTS};
use crate::wafer::Material;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Theme { Dark, Light }

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum LegendLayer { Substrate, DepositedFilm, Photoresist, MaskOpenings, EtchedRegion, DopantOverlay, Metrology, SelectionPreview, WarningError }

impl LegendLayer {
    pub const ALL: [LegendLayer; 9] = [LegendLayer::Substrate, LegendLayer::DepositedFilm, LegendLayer::Photoresist, LegendLayer::MaskOpenings,
        LegendLayer::EtchedRegion, LegendLayer::DopantOverlay, LegendLayer::Metrology, LegendLayer::SelectionPreview, LegendLayer::WarningError];
    pub fn index(self) -> usize { Self::ALL.iter().position(|l| *l == self).unwrap() }
    pub fn label(self) -> &'static str {
        match self {
            LegendLayer::Substrate => "Substrate", LegendLayer::DepositedFilm => "Deposited film", LegendLayer::Photoresist => "Photoresist",
            LegendLayer::MaskOpenings => "Mask / openings", LegendLayer::EtchedRegion => "Etched region", LegendLayer::DopantOverlay => "Dopant overlay",
            LegendLayer::Metrology => "Metrology", LegendLayer::SelectionPreview => "Selection / preview", LegendLayer::WarningError => "Warning / error",
        }
    }
    pub fn meaning(self) -> &'static str {
        match self {
            LegendLayer::Substrate => "Silicon or selected wafer base material",
            LegendLayer::DepositedFilm => "Material identity and film thickness",
            LegendLayer::Photoresist => "Resist-covered regions",
            LegendLayer::MaskOpenings => "Protected and exposed regions",
            LegendLayer::EtchedRegion => "Removed volume, trench or cavity",
            LegendLayer::DopantOverlay => "Concentration (one step per decade ≥ 1e17 cm⁻³)",
            LegendLayer::Metrology => "Measurement points, values, out-of-spec areas",
            LegendLayer::SelectionPreview => "Temporary selection and uncommitted geometry",
            LegendLayer::WarningError => "Invalid parameters, failed stages, unmet tolerances",
        }
    }
    /// Layers that exist as 3D material in the raymarcher (the rest are overlays).
    pub fn is_material(self) -> bool { matches!(self, LegendLayer::Substrate | LegendLayer::DepositedFilm | LegendLayer::Photoresist | LegendLayer::DopantOverlay) }
}

/// Layers relevant to a recipe stage (shown highlighted in the legend).
pub fn layers_for_stage(stage: Stage) -> &'static [LegendLayer] {
    use LegendLayer::*;
    match stage {
        Stage::SubstratePrep => &[Substrate, Metrology, WarningError],
        Stage::Deposition => &[Substrate, DepositedFilm, Metrology, SelectionPreview],
        Stage::Lithography => &[Substrate, DepositedFilm, Photoresist, MaskOpenings, SelectionPreview],
        Stage::Etching => &[Substrate, DepositedFilm, Photoresist, EtchedRegion, MaskOpenings],
        Stage::Doping => &[Substrate, DepositedFilm, Photoresist, DopantOverlay],
        Stage::FinalClean => &[Substrate, DepositedFilm, Photoresist, WarningError],
        Stage::Metrology => &[Substrate, DepositedFilm, EtchedRegion, DopantOverlay, Metrology, WarningError],
        Stage::Cmp => &[Substrate, DepositedFilm, Metrology, WarningError],
    }
}

pub type Rgb = [u8; 3];

#[derive(Debug, Clone)]
pub struct Palette {
    pub theme: Theme,
    pub background: Rgb,
    pub layers: [Rgb; 9],
    pub sio2: Rgb,
    pub sic: Rgb,
    pub gan: Rgb,
    pub metal: Rgb,
    pub visible: [bool; 9],
}

impl Palette {
    pub fn new(theme: Theme) -> Self {
        match theme {
            Theme::Dark => Self {
                theme, background: [11, 16, 24],
                layers: [[140, 160, 190], [80, 200, 240], [255, 165, 50], [245, 245, 245], [255, 110, 70], [230, 60, 200], [0, 230, 190], [60, 130, 255], [255, 40, 40]],
                sio2: [80, 200, 240], sic: [170, 130, 255], gan: [110, 220, 110], metal: [240, 230, 100], visible: [true; 9],
            },
            Theme::Light => Self {
                theme, background: [238, 242, 247],
                layers: [[80, 95, 120], [0, 130, 180], [200, 100, 0], [20, 20, 20], [220, 70, 60], [180, 0, 150], [0, 140, 120], [30, 80, 220], [200, 0, 0]],
                sio2: [0, 130, 180], sic: [110, 60, 200], gan: [30, 140, 50], metal: [140, 120, 0], visible: [true; 9],
            },
        }
    }

    pub fn layer(&self, l: LegendLayer) -> Rgb { self.layers[l.index()] }
    pub fn set_layer(&mut self, l: LegendLayer, c: Rgb) {
        self.layers[l.index()] = c;
        if l == LegendLayer::DepositedFilm { self.sio2 = c; }
    }
    pub fn is_visible(&self, l: LegendLayer) -> bool { self.visible[l.index()] }

    pub fn material_rgb(&self, m: Material) -> Rgb {
        match m {
            Material::Air => self.background, Material::Silicon => self.layer(LegendLayer::Substrate), Material::SiliconDioxide => self.sio2,
            Material::SiliconCarbide => self.sic, Material::GalliumNitride => self.gan, Material::Photoresist => self.layer(LegendLayer::Photoresist),
            Material::Metal => self.metal, Material::Dopant(b) => self.dopant_rgb(b),
        }
    }

    /// Dopant concentration ramp: one step per decade, fading toward the background at low concentration.
    pub fn dopant_rgb(&self, bucket: u8) -> Rgb {
        let base = self.layer(LegendLayer::DopantOverlay);
        let s = 0.45 + 0.55 * (bucket.min(6) as f64 / 6.0);
        let mix = |b: u8, bg: u8| (bg as f64 + (b as f64 - bg as f64) * s).round().clamp(0.0, 255.0) as u8;
        [mix(base[0], self.background[0]), mix(base[1], self.background[1]), mix(base[2], self.background[2])]
    }

    pub fn material_layer(m: Material) -> Option<LegendLayer> {
        match m {
            Material::Air => None, Material::Silicon => Some(LegendLayer::Substrate), Material::Photoresist => Some(LegendLayer::Photoresist),
            Material::Dopant(_) => Some(LegendLayer::DopantOverlay), _ => Some(LegendLayer::DepositedFilm),
        }
    }

    /// Shader palette: one RGBA per slot, plus a per-slot visibility flag (1.0 / 0.0).
    pub fn shader_tables(&self) -> ([[f32; 4]; PALETTE_SLOTS], [f32; PALETTE_SLOTS]) {
        let mut pal = [[0.0f32; 4]; PALETTE_SLOTS];
        let mut vis = [1.0f32; PALETTE_SLOTS];
        let mats = [Material::Air, Material::Silicon, Material::SiliconDioxide, Material::SiliconCarbide, Material::GalliumNitride, Material::Photoresist, Material::Metal];
        let alpha = |m: Material| match m { Material::Air => 0.0, Material::SiliconDioxide => 0.82, Material::Photoresist => 0.85, _ => 1.0 };
        for (slot, m) in mats.iter().enumerate() {
            let c = self.material_rgb(*m);
            pal[slot] = [c[0] as f32 / 255.0, c[1] as f32 / 255.0, c[2] as f32 / 255.0, alpha(*m)];
            if let Some(l) = Self::material_layer(*m) { vis[slot] = if self.is_visible(l) { 1.0 } else { 0.0 }; }
        }
        for b in 0..=8u8 {
            let slot = (DOPANT_ID_BASE + b as i32) as usize;
            let c = self.dopant_rgb(b);
            pal[slot] = [c[0] as f32 / 255.0, c[1] as f32 / 255.0, c[2] as f32 / 255.0, 1.0];
            vis[slot] = if self.is_visible(LegendLayer::DopantOverlay) { 1.0 } else { 0.0 };
        }
        (pal, vis)
    }
}

fn lin(c: u8) -> f64 { let v = c as f64 / 255.0; if v <= 0.03928 { v / 12.92 } else { ((v + 0.055) / 1.055).powf(2.4) } }
pub fn luminance(c: Rgb) -> f64 { 0.2126 * lin(c[0]) + 0.7152 * lin(c[1]) + 0.0722 * lin(c[2]) }
/// WCAG contrast ratio (1..21).
pub fn contrast_ratio(a: Rgb, b: Rgb) -> f64 { let (la, lb) = (luminance(a), luminance(b)); (la.max(lb) + 0.05) / (la.min(lb) + 0.05) }
pub fn rgb_distance(a: Rgb, b: Rgb) -> f64 { (0..3).map(|i| (a[i] as f64 - b[i] as f64).powi(2)).sum::<f64>().sqrt() }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_layer_is_contrast_safe_against_its_theme_background() {
        for th in [Theme::Dark, Theme::Light] {
            let p = Palette::new(th);
            for l in LegendLayer::ALL { assert!(contrast_ratio(p.layer(l), p.background) >= 3.0, "{th:?} {l:?}: {:.2}", contrast_ratio(p.layer(l), p.background)); }
            for c in [p.sio2, p.sic, p.gan, p.metal] { assert!(contrast_ratio(c, p.background) >= 3.0, "{th:?} film colour {c:?}"); }
            for b in 1..=6 { assert!(contrast_ratio(p.dopant_rgb(b), p.background) >= 1.5, "{th:?} dopant bucket {b}"); }
        }
    }

    #[test]
    fn layers_and_film_materials_are_mutually_distinguishable() {
        for th in [Theme::Dark, Theme::Light] {
            let p = Palette::new(th);
            let mut all: Vec<(String, Rgb)> = LegendLayer::ALL.iter().map(|l| (format!("{l:?}"), p.layer(*l))).collect();
            all.push(("SiC".into(), p.sic)); all.push(("GaN".into(), p.gan)); all.push(("Metal".into(), p.metal));
            all.retain(|(n, _)| n != "DepositedFilm");
            all.push(("SiO2".into(), p.sio2));
            for i in 0..all.len() { for j in i + 1..all.len() {
                assert!(rgb_distance(all[i].1, all[j].1) >= 45.0, "{th:?}: {} vs {} too similar ({:.0})", all[i].0, all[j].0, rgb_distance(all[i].1, all[j].1));
            } }
        }
    }

    #[test]
    fn dopant_ramp_is_monotonic_in_distance_from_background() {
        let p = Palette::new(Theme::Dark);
        let d = |b| rgb_distance(p.dopant_rgb(b), p.background);
        for b in 1..6 { assert!(d(b + 1) > d(b)); }
    }

    #[test]
    fn colours_are_configurable_and_reach_the_shader_tables() {
        let mut p = Palette::new(Theme::Dark);
        p.set_layer(LegendLayer::Photoresist, [1, 2, 3]);
        let (pal, vis) = p.shader_tables();
        assert_eq!(pal[5][..3], [1.0 / 255.0, 2.0 / 255.0, 3.0 / 255.0]);
        assert_eq!(vis[5], 1.0);
        p.visible[LegendLayer::Photoresist.index()] = false;
        p.visible[LegendLayer::DopantOverlay.index()] = false;
        let (_, vis) = p.shader_tables();
        assert_eq!((vis[5], vis[7], vis[10], vis[1]), (0.0, 0.0, 0.0, 1.0));
        assert_eq!(pal.len(), PALETTE_SLOTS);
    }

    #[test]
    fn every_stage_names_its_active_layers() {
        for s in Stage::ALL { assert!(!layers_for_stage(s).is_empty()); }
        assert!(layers_for_stage(Stage::Doping).contains(&LegendLayer::DopantOverlay));
        assert!(layers_for_stage(Stage::Etching).contains(&LegendLayer::EtchedRegion));
    }
}
