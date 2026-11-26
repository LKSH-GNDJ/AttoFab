// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Mask geometry and the aerial image.
//!
//! Two audit defects are fixed here. The mask is sampled with 4x4 area
//! integration per pixel rather than at the pixel centre (audit F-22), and the
//! develop threshold is applied to a normalised image with the mask polarity
//! carried explicitly, so a dark line on a clear field and a clear line on a
//! dark field are not interchangeable (audit F-01).

use atto_units::{Nm};
use crate::column::ColumnGrid;

#[derive(Debug, Clone, PartialEq)]
/// A mask, either a drawn pattern or a uniform field.
pub enum MaskPattern {
    /// Equal lines centred on the mask origin, repeating along X when `vertical`.
    LineSpace { line_nm: f64, pitch_nm: f64, count: u32, vertical: bool, clear_field: bool },
}

impl MaskPattern {
    /// Whether this is a clear field rather than a dark one.
    pub fn clear_field(&self) -> bool {
        match self { MaskPattern::LineSpace { clear_field, .. } => *clear_field }
    }

    /// Check the pattern can be imaged, returning why it cannot.
    pub fn validate(&self) -> Result<(), String> {
        match self {
            MaskPattern::LineSpace { line_nm, pitch_nm, count, .. } => {
                if !(line_nm.is_finite() && *line_nm > 0.0) { return Err("line width must be > 0".into()); }
                if !(pitch_nm.is_finite() && *pitch_nm > 0.0) { return Err("pitch must be > 0".into()); }
                if *line_nm > *pitch_nm { return Err("line width exceeds pitch".into()); }
                if *count == 0 { return Err("mask needs at least one line".into()); }
            }
        }
        Ok(())
    }

    /// Is this point covered by an opaque feature?
    pub fn is_feature(&self, u: f64, v: f64) -> bool {
        match self {
            MaskPattern::LineSpace { line_nm, pitch_nm, count, vertical, .. } => {
                let across = if *vertical { u } else { v };
                let k = (across / pitch_nm + (*count as f64 - 1.0) * 0.5).round();
                if k < 0.0 || k > *count as f64 - 1.0 { return false; }
                let centre = (k - (*count as f64 - 1.0) * 0.5) * pitch_nm;
                (across - centre).abs() <= line_nm * 0.5
            }
        }
    }

    /// Transmission at a point: 1 where light passes.
    pub fn transmission(&self, u: f64, v: f64) -> f64 {
        let f = self.is_feature(u, v);
        if self.clear_field() { if f { 0.0 } else { 1.0 } } else if f { 1.0 } else { 0.0 }
    }

    /// Area-integrated transmission over a pixel of edge `px`, 4x4 samples.
    pub fn transmission_area(&self, u: f64, v: f64, px: f64) -> f64 {
        const N: usize = 4;
        let mut acc = 0.0;
        for i in 0..N {
            for j in 0..N {
                let du = (i as f64 + 0.5) / N as f64 * px - px * 0.5;
                let dv = (j as f64 + 0.5) / N as f64 * px - px * 0.5;
                acc += self.transmission(u + du, v + dv);
            }
        }
        acc / (N * N) as f64
    }
}

/// Normalised aerial image on the grid: 0 is fully dark, 1 is the open field.
/// A separable Gaussian of `sigma_nm` stands in for the optical system; the
/// production model replaces this with Fourier imaging (build PRD L-01).
pub fn aerial_image(g: &ColumnGrid, mask: &MaskPattern, origin: [f64; 2], sigma_nm: f64) -> Result<Vec<f32>, String> {
    mask.validate()?;
    if !(sigma_nm.is_finite() && sigma_nm >= 0.0) { return Err("blur must be finite and >= 0".into()); }
    let pad = if sigma_nm > 0.0 { (3.0 * sigma_nm / g.pitch).ceil() as usize } else { 0 };
    let (wx, wy) = (g.nx + 2 * pad, g.ny + 2 * pad);
    let mut t = vec![0.0f64; wx * wy];
    for j in 0..wy {
        for i in 0..wx {
            let x = g.x0 + (i as f64 - pad as f64 + 0.5) * g.pitch - origin[0];
            let y = g.y0 + (j as f64 - pad as f64 + 0.5) * g.pitch - origin[1];
            t[j * wx + i] = mask.transmission_area(x, y, g.pitch);
        }
    }
    if pad > 0 {
        let sig = sigma_nm / g.pitch;
        let kern: Vec<f64> = (-(pad as isize)..=pad as isize)
            .map(|k| (-(k * k) as f64 / (2.0 * sig * sig)).exp())
            .collect();
        let ksum: f64 = kern.iter().sum();
        let mut tmp = t.clone();
        for j in 0..wy {
            for i in 0..wx {
                let mut acc = 0.0;
                for (m, kv) in kern.iter().enumerate() {
                    let ii = (i as isize + m as isize - pad as isize).clamp(0, wx as isize - 1) as usize;
                    acc += kv * t[j * wx + ii];
                }
                tmp[j * wx + i] = acc / ksum;
            }
        }
        let mut out = tmp.clone();
        for j in 0..wy {
            for i in 0..wx {
                let mut acc = 0.0;
                for (m, kv) in kern.iter().enumerate() {
                    let jj = (j as isize + m as isize - pad as isize).clamp(0, wy as isize - 1) as usize;
                    acc += kv * tmp[jj * wx + i];
                }
                out[j * wx + i] = acc / ksum;
            }
        }
        t = out;
    }
    let mut img = vec![0.0f32; g.nx * g.ny];
    for j in 0..g.ny {
        for i in 0..g.nx {
            img[j * g.nx + i] = t[(j + pad) * wx + (i + pad)] as f32;
        }
    }
    Ok(img)
}

/// Columns whose resist is removed. `dose_ratio` is the exposure dose divided by
/// the dose to clear, so the threshold on the normalised image is 1/dose_ratio.
pub fn develop_map(intensity: &[f32], dose_ratio: f64, positive: bool) -> Vec<bool> {
    let thr = (1.0 / dose_ratio.max(1e-6)) as f32;
    intensity.iter().map(|&i| if positive { i >= thr } else { i < thr }).collect()
}

/// Rayleigh half-pitch in nm.
pub fn rayleigh_nm(k1: f64, lambda: Nm, na: f64) -> Result<f64, String> {
    let lambda_nm = lambda.value();
    if !(k1.is_finite() && k1 > 0.0) { return Err("k1 must be > 0".into()); }
    if !(lambda_nm.is_finite() && lambda_nm > 0.0) { return Err("wavelength must be > 0".into()); }
    if !(na.is_finite() && na > 0.0) { return Err("NA must be > 0".into()); }
    Ok(k1 * lambda_nm / na)
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;
    use crate::column::ProcessWindow;
    use crate::svo::Wafer;

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        ColumnGrid::sample(&w, &ProcessWindow { center: [0.0, 0.0], size: [1200.0, 1200.0] }).unwrap()
    }

    #[test]
    fn area_sampling_antialiases_an_edge() {
        let m = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, vertical: true, clear_field: true };
        // a point half a pixel inside the edge must read between 0 and 1
        let edge = 75.0;
        let v = m.transmission_area(edge, 0.0, 16.0);
        assert!(v > 0.0 && v < 1.0, "an edge inside a pixel must integrate to a fraction, got {v}");
        let centre = m.transmission_area(0.0, 0.0, 16.0);
        assert!(centre < 0.01, "the middle of a 150 nm line is fully covered, got {centre}");
    }

    #[test]
    fn dark_line_on_clear_field_retains_the_line() {
        let g = grid();
        let m = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, vertical: true, clear_field: true };
        let r = rayleigh_nm(0.4, Nm::new(193.0), 0.85).unwrap();
        let sigma = r / 2.355;
        let img = aerial_image(&g, &m, [0.0, 0.0], sigma).unwrap();
        let removed = develop_map(&img, 2.0, true);
        let row = g.ny / 2;
        let at = |x_nm: f64| g.idx(((x_nm - g.x0) / g.pitch) as usize, row);
        assert!(!removed[at(0.0)], "resist retained under the opaque line");
        assert!(removed[at(160.0)], "the space between lines opens");
    }

    #[test]
    fn printed_cd_tracks_the_drawn_width() {
        let g = grid();
        let r = rayleigh_nm(0.4, Nm::new(193.0), 0.85).unwrap();
        let sigma = r / 2.355;
        let mut cds = Vec::new();
        for drawn in [120.0, 150.0, 200.0] {
            let m = MaskPattern::LineSpace { line_nm: drawn, pitch_nm: 400.0, count: 1, vertical: true, clear_field: true };
            let img = aerial_image(&g, &m, [0.0, 0.0], sigma).unwrap();
            let removed = develop_map(&img, 2.0, true);
            let row = g.ny / 2;
            let n = (0..g.nx).filter(|i| !removed[g.idx(*i, row)]).count();
            cds.push((drawn, n as f64 * g.pitch));
        }
        for (drawn, got) in cds {
            assert!((got - drawn).abs() <= g.pitch, "drawn {drawn} printed {got}");
        }
    }

    #[test]
    fn polarity_changes_the_result() {
        let g = grid();
        let r = rayleigh_nm(0.4, Nm::new(193.0), 0.85).unwrap();
        let sigma = r / 2.355;
        let clear = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, vertical: true, clear_field: true };
        let dark = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, vertical: true, clear_field: false };
        let a = aerial_image(&g, &clear, [0.0, 0.0], sigma).unwrap();
        let b = aerial_image(&g, &dark, [0.0, 0.0], sigma).unwrap();
        let row = g.ny / 2;
        let k = g.idx(g.nx / 2, row);
        assert!((a[k] - (1.0 - b[k])).abs() < 1e-5, "the two masks are complementary");
    }

    #[test]
    fn higher_dose_widens_the_opening() {
        let g = grid();
        let m = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, vertical: true, clear_field: true };
        let img = aerial_image(&g, &m, [0.0, 0.0], 20.0).unwrap();
        let n = |d: f64| develop_map(&img, d, true).iter().filter(|r| **r).count();
        assert!(n(3.0) >= n(1.5), "more dose clears more resist");
    }
}
