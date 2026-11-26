// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Fourier-optics imaging and the Bossung curve (build PRD L-01, L-04, L-07).
//!
//! A periodic mask is imaged by summing its diffraction orders inside the pupil,
//! each carrying a defocus phase, and averaging the resulting intensity over the
//! source points that make up partial coherence. That averaging is what produces
//! the focus sensitivity: with a single coherent source point, defocus merely
//! shifts the fringes and the contrast never falls, which is the mistake a naive
//! model makes.

use atto_units::{Nm, Um};
use std::f64::consts::PI;

#[derive(Debug, Clone, PartialEq)]
/// Why an image or an exposure was refused.
pub enum LithoError {
    BadInput(String),
}

impl std::fmt::Display for LithoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            LithoError::BadInput(s) => write!(f, "{s}"),
        }
    }
}
impl std::error::Error for LithoError {}

/// The illumination and projection system.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Optics {
    pub lambda_nm: f64,
    pub na: f64,
    /// Partial coherence: source radius as a fraction of the pupil.
    pub sigma: f64,
    /// Distance from best focus along the axis, nm.
    pub defocus_nm: f64,
}

impl Default for Optics {
    fn default() -> Self {
        Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 }
    }
}

/// A one-dimensional line/space grating.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Grating {
    pub pitch_nm: f64,
    /// Width of the opaque (or clear) line within the pitch.
    pub line_nm: f64,
}

impl Grating {
    /// Fourier coefficient of order `m` for a 1:1-style grating with duty `d`.
    fn coefficient(&self, m: i32) -> f64 {
        let duty = (self.line_nm / self.pitch_nm).clamp(0.0, 1.0);
        if m == 0 {
            duty
        } else {
            (PI * m as f64 * duty).sin() / (PI * m as f64)
        }
    }
}

impl Optics {
    fn validate(&self, g: &Grating) -> Result<(), LithoError> {
        if !(self.lambda_nm > 0.0) || !(self.na > 0.0 && self.na <= 1.0) {
            return Err(LithoError::BadInput("wavelength and NA must be positive".into()));
        }
        if !(g.pitch_nm > 0.0) || !(g.line_nm > 0.0) {
            return Err(LithoError::BadInput("pitch and line width must be positive".into()));
        }
        if !(0.0..=1.0).contains(&self.sigma) {
            return Err(LithoError::BadInput("sigma must be between 0 and 1".into()));
        }
        Ok(())
    }

    /// Highest diffraction order that fits in the pupil for a source point `s`.
    fn max_order(&self, g: &Grating, s: f64) -> i32 {
        let cutoff = self.na / self.lambda_nm;
        let mut m = 0;
        loop {
            let next = m + 1;
            if (next as f64 / g.pitch_nm + s * cutoff).abs() > cutoff {
                return m;
            }
            m = next;
            if m > 60 {
                return m;
            }
        }
    }

    /// Intensity across one period, normalised to the clear-field intensity.
    pub fn image(&self, g: &Grating, samples: usize) -> Result<Vec<f64>, LithoError> {
        self.validate(g)?;
        if samples < 8 {
            return Err(LithoError::BadInput("need at least 8 samples".into()));
        }
        let cutoff = self.na / self.lambda_nm;
        // source points across the partial-coherence disc (1-D, on axis)
        let ns = if self.sigma > 0.0 { 17 } else { 1 };
        let mut out = vec![0.0; samples];
        let mut norm = 0.0;
        for k in 0..ns {
            let frac = if ns == 1 { 0.0 } else { (k as f64 / (ns - 1) as f64) * 2.0 - 1.0 };
            let s = frac * self.sigma;
            let m_max = self.max_order(g, s);
            // pupil transmission of each order, with its defocus phase
            let mut orders: Vec<(i32, f64, f64)> = Vec::new();
            for m in -m_max..=m_max {
                let f = m as f64 / g.pitch_nm + s * cutoff;
                if f.abs() > cutoff {
                    continue;
                }
                // Fresnel defocus phase: pi * lambda * z * f^2
                let phase = PI * self.lambda_nm * self.defocus_nm * f * f;
                orders.push((m, g.coefficient(m), phase));
            }
            for (i, o) in out.iter_mut().enumerate() {
                let x = i as f64 / samples as f64;
                let (mut re, mut im) = (0.0, 0.0);
                for (m, a, ph) in &orders {
                    let th = 2.0 * PI * (*m as f64) * x + ph;
                    re += a * th.cos();
                    im += a * th.sin();
                }
                *o += re * re + im * im;
            }
            norm += 1.0;
        }
        for v in out.iter_mut() {
            *v /= norm;
        }
        Ok(out)
    }

    /// Image contrast, `(Imax - Imin) / (Imax + Imin)`.
    pub fn contrast(&self, g: &Grating) -> Result<f64, LithoError> {
        let img = self.image(g, 1024)?;
        let mx = img.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let mn = img.iter().cloned().fold(f64::INFINITY, f64::min);
        if mx + mn <= 0.0 {
            return Ok(0.0);
        }
        Ok((mx - mn) / (mx + mn))
    }

    /// Normalised image log slope at the point where the image crosses
    /// `threshold`, per nm.
    pub fn nils(&self, g: &Grating, threshold: f64) -> Result<f64, LithoError> {
        let n = 512;
        let img = self.image(g, n)?;
        let dx = g.pitch_nm / n as f64;
        for i in 1..n {
            if (img[i - 1] - threshold) * (img[i] - threshold) < 0.0 {
                let slope = (img[i] - img[i - 1]) / dx;
                return Ok((slope / threshold).abs() * g.pitch_nm);
            }
        }
        Ok(0.0)
    }
}

/// An imaged layout: the mask, the aerial image on the same grid, and what
/// prints.
#[derive(Debug, Clone, PartialEq)]
pub struct LayoutImage {
    /// Mask transmission, 1 where the film is clear.
    pub mask: Vec<f64>,
    /// Aerial image intensity, normalised so a clear field is 1.0.
    pub intensity: Vec<f64>,
    pub contrast: f64,
    /// Total printed width across the field (the sum of every feature), nm.
    /// Divide by the feature count for a per-feature CD.
    pub printed_nm: f64,
    /// Field width, nm.
    pub field_nm: f64,
}

/// Image a layout layer through the projection system (build PRD L-02, L-03).
///
/// The mask row is Fourier transformed, each spatial frequency is multiplied by
/// the pupil (cut off at NA/lambda, with a Fresnel defocus phase), and the
/// result is transformed back. Partial coherence is handled by averaging over
/// source points, each shifting the pupil, which is what makes focus matter.
#[allow(clippy::too_many_arguments)]
pub fn image_layout(
    layout: &atto_geom::gds::Layout,
    layer: i16,
    optics: &Optics,
    field_nm: f64,
    y_nm: f64,
    nx: usize,
    threshold: f64,
) -> Result<LayoutImage, LithoError> {
    if !(field_nm > 0.0) || nx < 16 || !(optics.lambda_nm > 0.0) || !(optics.na > 0.0) {
        return Err(LithoError::BadInput("invalid field or optics".into()));
    }
    let pitch_um = field_nm / nx as f64 / 1000.0;
    let x0_um = 0.0;
    let y_um = (y_nm - pitch_um * 500.0) / 1000.0;
    let mask = layout.rasterise(layer, x0_um, y_um, atto_units::Um::new(pitch_um), nx, 1);

    // forward transform of the mask row
    let mut re = vec![0.0; nx];
    let mut im = vec![0.0; nx];
    for k in 0..nx {
        let fk = bin_frequency(k, nx, field_nm);
        let phase = PI * optics.lambda_nm * optics.defocus_nm * fk * fk;
        let (mut sr, mut si) = (0.0, 0.0);
        for (i, m) in mask.iter().enumerate() {
            let th = 2.0 * PI * (k as f64) * (i as f64) / (nx as f64);
            sr += m * th.cos();
            si -= m * th.sin();
        }
        re[k] = sr * phase.cos() - si * phase.sin();
        im[k] = sr * phase.sin() + si * phase.cos();
    }

    let cutoff = optics.na / optics.lambda_nm;
    let ns = if optics.sigma > 0.0 { 17 } else { 1 };
    let mut intensity = vec![0.0; nx];
    let mut norm = 0.0;
    for sidx in 0..ns {
        let frac = if ns == 1 { 0.0 } else { (sidx as f64 / (ns - 1) as f64) * 2.0 - 1.0 };
        let shift = frac * optics.sigma * cutoff;
        for i in 0..nx {
            let (mut ar, mut ai) = (0.0, 0.0);
            for k in 0..nx {
                let fk = bin_frequency(k, nx, field_nm);
                if (fk - shift).abs() > cutoff {
                    continue;
                }
                let th = 2.0 * PI * (k as f64) * (i as f64) / (nx as f64);
                ar += re[k] * th.cos() - im[k] * th.sin();
                ai += re[k] * th.sin() + im[k] * th.cos();
            }
            intensity[i] += ar * ar + ai * ai;
        }
        norm += 1.0;
    }
    // normalise so a fully clear mask reads 1.0
    let scale = (nx as f64) * (nx as f64) * norm;
    for v in intensity.iter_mut() {
        *v /= scale;
    }
    let mx = intensity.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mn = intensity.iter().cloned().fold(f64::INFINITY, f64::min);
    let contrast = if mx + mn > 0.0 { (mx - mn) / (mx + mn) } else { 0.0 };
    let over = intensity.iter().filter(|v| **v > threshold).count();
    let printed_nm = over as f64 * (field_nm / nx as f64);
    Ok(LayoutImage { mask, intensity, contrast, printed_nm, field_nm })
}

/// Expose a layout through the optics and return the resist pattern.
///
/// `true` means resist remains (a positive resist, so where the image is below
/// the threshold). The image is taken along one row and applied to every row,
/// which is exact for a line/space layout and an approximation for anything
/// that varies in y — stated rather than hidden.
#[allow(clippy::too_many_arguments)]
pub fn expose_layout(
    layout: &atto_geom::gds::Layout,
    layer: i16,
    optics: &Optics,
    field_nm: f64,
    y_nm: f64,
    nx: usize,
    threshold: f64,
    positive: bool,
) -> Result<Vec<bool>, LithoError> {
    let img = image_layout(layout, layer, optics, field_nm, y_nm, nx, threshold)?;
    Ok(img
        .intensity
        .iter()
        .map(|v| if positive { *v < threshold } else { *v > threshold })
        .collect())
}

/// Write a resist pattern onto a grid as a film of `thickness_nm`.
///
/// The pattern is sampled onto the grid's own columns, so a grid coarser than
/// the image is resampled rather than mismatched.
pub fn apply_resist(
    grid: &mut atto_geom::column::ColumnGrid,
    pattern: &[bool],
    field_nm: f64,
    thickness_nm: f64,
) {
    use atto_geom::column::set_interval;
    use atto_geom::svo::Material;
    if pattern.is_empty() || !(field_nm > 0.0) || !(thickness_nm > 0.0) {
        return;
    }
    for k in 0..grid.cols.len() {
        let ix = k % grid.nx;
        let x_nm = grid.x0 + (ix as f64 + 0.5) * grid.pitch;
        let frac = (x_nm / field_nm).clamp(0.0, 0.999999);
        let idx = (frac * pattern.len() as f64) as usize;
        if !pattern.get(idx).copied().unwrap_or(false) {
            continue;
        }
        let top = match atto_geom::column::solid_top(&grid.cols[k]) {
            Some((t, _)) => t,
            None => continue,
        };
        set_interval(&mut grid.cols[k], top, top + thickness_nm, Material::Photoresist);
    }
}

/// Proximity-effect correction parameters (build PRD L-05, audit F-06).
///
/// Electrons scattered in the resist deposit energy over a range far wider than
/// the beam, so a dense array receives extra background dose that an isolated
/// feature does not. the earlier build corrected this with one global constant; here the
/// background is the local pattern density convolved with a Gaussian, and the
/// dose is scaled down where that density is high.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pec {
    /// Range of the backscattered energy, nm.
    pub sigma_nm: f64,
    /// Ratio of backscattered to forward energy.
    pub eta: f64,
}

impl Default for Pec {
    fn default() -> Self {
        Pec { sigma_nm: 1000.0, eta: 0.5 }
    }
}

/// Local pattern density: the mask convolved with a Gaussian of `sigma_nm`.
///
/// The kernel is normalised, so a fully clear field gives 1.0 and an empty one
/// gives 0.0 regardless of the range.
pub fn proximity_density(mask: &[f64], field: Nm, sigma_nm: f64) -> Vec<f64> {
    let field_nm = field.value();
    let n = mask.len();
    if n == 0 || !(field_nm > 0.0) || !(sigma_nm > 0.0) {
        return vec![0.0; n];
    }
    let dx = field_nm / n as f64;
    let r = ((3.0 * sigma_nm) / dx).ceil() as isize;
    let mut kern: Vec<f64> = (-r..=r)
        .map(|i| {
            let d = i as f64 * dx;
            (-(d * d) / (2.0 * sigma_nm * sigma_nm)).exp()
        })
        .collect();
    let ksum: f64 = kern.iter().sum();
    for k in kern.iter_mut() {
        *k /= ksum;
    }
    let mut out = vec![0.0; n];
    for i in 0..n {
        // Renormalise by the weight that actually lands inside the field: the
        // kernel is truncated at the edges, and dividing by its full sum makes a
        // clear field read 0.7 at the boundary instead of 1.0.
        let mut acc = 0.0;
        let mut used = 0.0;
        for (j, w) in kern.iter().enumerate() {
            let x = i as isize + j as isize - r;
            if x >= 0 && (x as usize) < n {
                acc += w * mask[x as usize];
                used += w;
            }
        }
        out[i] = if used > 0.0 { (acc / used).clamp(0.0, 1.0) } else { 0.0 };
    }
    out
}

/// Dose multiplier that compensates the background, `1 / (1 + eta*rho)`.
pub fn pec_dose(density: &[f64], eta: f64) -> Vec<f64> {
    let e = eta.max(0.0);
    density
        .iter()
        .map(|rho| {
            let d = 1.0 + e * rho.clamp(0.0, 1.0);
            if d > 0.0 { 1.0 / d } else { 1.0 }
        })
        .collect()
}

/// Spatial frequency of transform bin `k`, cycles per nm.
fn bin_frequency(k: usize, nx: usize, field_nm: f64) -> f64 {
    let signed = if k <= nx / 2 { k as f64 } else { k as f64 - nx as f64 };
    signed / field_nm
}

/// Critical dimension from a threshold develop model, nm.
///
/// The printed width is the span of the image that lies above (positive resist)
/// or below (negative resist) the threshold.
pub fn cd_from_threshold(img: &[f64], pitch: Nm, threshold: f64, positive: bool) -> f64 {
    let pitch_nm = pitch.value();
    let n = img.len();
    if n < 8 {
        return 0.0;
    }
    let dx = pitch_nm / n as f64;
    let over: Vec<bool> = img.iter().map(|v| if positive { *v > threshold } else { *v < threshold }).collect();
    let count = over.iter().filter(|b| **b).count();
    // a single contiguous run is the usual case; the count is the width either way
    let _ = over;
    count as f64 * dx
}

#[cfg(test)]
mod pec_tests {
    use atto_units::{MTorr, Mm, Watts};
    use super::*;

    #[test]
    fn the_kernel_is_normalised() {
        // a fully clear mask has density 1 everywhere, whatever the range
        let mask = vec![1.0; 64];
        for sigma in [100.0, 500.0, 2000.0] {
            let d = proximity_density(&mask, Nm::new(6400.0), sigma);
            assert!(d.iter().all(|v| (v - 1.0).abs() < 1e-9), "sigma {sigma}: {:?}", &d[..3]);
        }
        // and an empty one has none
        let empty = vec![0.0; 64];
        assert!(proximity_density(&empty, Nm::new(6400.0), 500.0).iter().all(|v| *v < 1e-12));
    }

    #[test]
    fn a_dense_array_has_more_background_than_an_isolated_line() {
        let mut sparse = vec![0.0; 128];
        sparse[64] = 1.0;
        let mut dense = vec![0.0; 128];
        for i in 0..128 {
            if i % 2 == 0 {
                dense[i] = 1.0;
            }
        }
        let ds = proximity_density(&sparse, Nm::new(12800.0), 1000.0);
        let dd = proximity_density(&dense, Nm::new(12800.0), 1000.0);
        assert!(dd[64] > 10.0 * ds[64], "dense {:.3} vs isolated {:.4}", dd[64], ds[64]);
        assert!(ds[64] < 0.05, "an isolated line has almost no background: {:.4}", ds[64]);
    }

    #[test]
    fn the_dose_falls_where_the_background_is_high() {
        let dense = vec![0.5; 8];
        let d = pec_dose(&dense, 0.5);
        assert!(d.iter().all(|v| *v < 1.0), "a dense area must be dosed down");
        assert!((d[0] - 1.0 / 1.25).abs() < 1e-12);
        // no eta, no correction
        assert!(pec_dose(&dense, 0.0).iter().all(|v| (v - 1.0).abs() < 1e-12));
    }

    #[test]
    fn a_low_density_gets_no_correction() {
        let iso = vec![0.001; 8];
        let d = pec_dose(&iso, 0.5);
        assert!(d.iter().all(|v| (v - 1.0).abs() < 0.01), "isolated features are left alone");
    }

    #[test]
    fn bad_input_is_rejected() {
        assert!(proximity_density(&[1.0, 0.0], Nm::new(0.0), 100.0).iter().all(|v| *v == 0.0));
        assert!(proximity_density(&[1.0, 0.0], Nm::new(100.0), 0.0).iter().all(|v| *v == 0.0));
        assert!(proximity_density(&[], Nm::new(100.0), 10.0).is_empty());
    }
}

#[cfg(test)]
mod end_to_end_tests {
    //! GDS layout through optics, resist and plasma to a cross-section.
    use atto_units::{MTorr, Mm, Watts};
    use super::*;
    use atto_geom::column::{etch, solid_top, ColumnGrid, EtchSpec, ProcessWindow};
    use atto_geom::gds::{Layout, Polygon, Structure};
    use atto_geom::svo::{Material, Wafer};

    fn lines() -> Layout {
        let mut polygons = Vec::new();
        for i in 0..6 {
            let x = i as f64 * 0.5;
            polygons.push(Polygon {
                layer: 1,
                datatype: 0,
                points: vec![(x, 0.0), (x + 0.25, 0.0), (x + 0.25, 1.0), (x, 1.0), (x, 0.0)],
            });
        }
        Layout { units_per_micron: 1000.0, structures: vec![Structure { name: "TOP".into(), polygons }] }
    }

    #[test]
    fn a_layout_becomes_a_cross_section() {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let mut g = ColumnGrid::sample(&w, &win).unwrap();

        // oxide on top, then resist from the layout
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            atto_geom::column::set_interval(&mut g.cols[k], top, top + 120.0, Material::SiliconDioxide);
        }
        let o = Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 };
        let pattern = expose_layout(&lines(), 1, &o, 3000.0, 500.0, 256, 0.35, true).unwrap();
        apply_resist(&mut g, &pattern, 3000.0, 80.0);

        let model = crate::etch::EtchModel::bundled_oxide();
        let rates = model.rates(Watts::new(200.0), MTorr::new(20.0));
        let spec = EtchSpec { rates, time_min: 1.0 };
        let (etched, stats) = etch(&g, &spec, 1.0).expect("etch runs");

        // some columns are protected by resist, some are open, and the open ones
        // lost oxide
        assert!(stats.feature_depth_max_nm > 0.0, "something must etch");
        let mut protected = 0;
        let mut opened = 0;
        for k in 0..g.cols.len() {
            let had_resist = matches!(solid_top(&g.cols[k]), Some((_, Material::Photoresist)));
            let now = solid_top(&etched.cols[k]).map(|(_, m)| m);
            if had_resist {
                if now == Some(Material::Photoresist) {
                    protected += 1;
                }
            } else if now == Some(Material::SiliconDioxide) || now == Some(Material::Silicon) {
                opened += 1;
            }
        }
        assert!(protected > 0, "resist-protected columns must survive");
        assert!(opened > 0, "open columns must be reached");
        assert!(stats.feature_depth_max_nm < 200.0, "the etch must not run away");
    }

    #[test]
    fn a_longer_etch_goes_deeper() {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let mut g = ColumnGrid::sample(&w, &win).unwrap();
        for k in 0..g.cols.len() {
            let top = solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            atto_geom::column::set_interval(&mut g.cols[k], top, top + 200.0, Material::SiliconDioxide);
        }
        let o = Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 };
        let pattern = expose_layout(&lines(), 1, &o, 3000.0, 500.0, 128, 0.35, true).unwrap();
        apply_resist(&mut g, &pattern, 3000.0, 80.0);
        let model = crate::etch::EtchModel::bundled_oxide();
        let rates = model.rates(Watts::new(200.0), MTorr::new(20.0));
        let (a, _) = etch(&g, &EtchSpec { rates, time_min: 0.5 }, 1.0).unwrap();
        let (b, _) = etch(&g, &EtchSpec { rates, time_min: 1.5 }, 1.0).unwrap();
        use atto_measure::etch_depth_nm;
        let da = etch_depth_nm(&g, &a).unwrap_or(0.0);
        let db = etch_depth_nm(&g, &b).unwrap_or(0.0);
        assert!(db > da, "three times the time must go deeper: {da:.1} -> {db:.1} nm");
    }
}

#[cfg(test)]
mod resist_tests {
    use atto_units::{MTorr, Mm, Watts};
    use super::*;
    use atto_geom::column::{solid_top, ColumnGrid, ProcessWindow};
    use atto_geom::svo::{Material, Wafer};
    use atto_geom::gds::{Layout, Polygon, Structure};

    fn lines() -> Layout {
        let mut polygons = Vec::new();
        for i in 0..6 {
            let x = i as f64 * 0.5;
            polygons.push(Polygon {
                layer: 1,
                datatype: 0,
                points: vec![(x, 0.0), (x + 0.25, 0.0), (x + 0.25, 1.0), (x, 1.0), (x, 0.0)],
            });
        }
        Layout { units_per_micron: 1000.0, structures: vec![Structure { name: "TOP".into(), polygons }] }
    }

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        ColumnGrid::sample(&w, &win).unwrap()
    }

    #[test]
    fn a_positive_resist_keeps_the_dark_half_of_the_image() {
        let o = Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 };
        let p = expose_layout(&lines(), 1, &o, 3000.0, 500.0, 256, 0.35, true).unwrap();
        let kept = p.iter().filter(|b| **b).count();
        // 6 lines of 250 nm in a 3000 nm field: half the field is dark
        let frac = kept as f64 / p.len() as f64;
        assert!((frac - 0.5).abs() < 0.12, "resist covers {frac:.2} of the field");
        assert!(p.iter().any(|b| *b) && p.iter().any(|b| !*b), "the pattern must vary");
    }

    #[test]
    fn a_negative_resist_keeps_the_complement() {
        let o = Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 };
        let pos = expose_layout(&lines(), 1, &o, 3000.0, 500.0, 256, 0.35, true).unwrap();
        let neg = expose_layout(&lines(), 1, &o, 3000.0, 500.0, 256, 0.35, false).unwrap();
        for (a, b) in pos.iter().zip(neg.iter()) {
            assert_ne!(a, b, "the two polarities must be complementary");
        }
    }

    #[test]
    fn the_resist_lands_on_the_grid_as_a_film() {
        let o = Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 };
        let p = expose_layout(&lines(), 1, &o, 3000.0, 500.0, 256, 0.35, true).unwrap();
        let mut g = grid();
        let before: Vec<Option<(f64, Material)>> = g.cols.iter().map(|c| solid_top(c)).collect();
        apply_resist(&mut g, &p, 3000.0, 80.0);
        let mut coated = 0;
        for (k, b) in before.iter().enumerate() {
            let top = solid_top(&g.cols[k]);
            if let (Some((tb, _)), Some((ta, m))) = (b, top) {
                if m == Material::Photoresist && ta > *tb + 1.0 {
                    coated += 1;
                }
            }
        }
        assert!(coated > 0, "some columns must carry resist");
        assert!(coated < g.cols.len(), "not all of them");
    }

    #[test]
    fn an_empty_pattern_changes_nothing() {
        let mut g = grid();
        let before: Vec<usize> = g.cols.iter().map(|c| c.len()).collect();
        apply_resist(&mut g, &[], 3000.0, 80.0);
        let after: Vec<usize> = g.cols.iter().map(|c| c.len()).collect();
        assert_eq!(before, after);
    }
}

#[cfg(test)]
mod layout_tests {
    use atto_units::{MTorr, Mm, Watts};
    use super::*;
    use atto_geom::gds::{Layout, Polygon, Structure};

    /// Ten 125 nm lines on a 250 nm pitch, in microns.
    fn lines() -> Layout {
        let mut polygons = Vec::new();
        for i in 0..10 {
            let x = i as f64 * 0.25;
            polygons.push(Polygon {
                layer: 1,
                datatype: 0,
                points: vec![(x, 0.0), (x + 0.125, 0.0), (x + 0.125, 1.0), (x, 1.0), (x, 0.0)],
            });
        }
        Layout { units_per_micron: 1000.0, structures: vec![Structure { name: "TOP".into(), polygons }] }
    }

    fn optics(z: f64) -> Optics {
        Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: z }
    }

    #[test]
    fn a_layout_images_to_a_pattern_with_contrast() {
        let r = image_layout(&lines(), 1, &optics(0.0), 2500.0, 500.0, 256, 0.35).unwrap();
        assert!(r.contrast > 0.7, "contrast {:.3}", r.contrast);
        assert!(r.intensity.iter().all(|v| v.is_finite() && *v >= -1e-9));
        // half the field is drawn, so about half the field prints
        let drawn = 10.0 * 125.0;
        assert!((r.printed_nm - drawn).abs() / drawn < 0.15, "printed {:.0} vs drawn {:.0} nm", r.printed_nm, drawn);
    }

    #[test]
    fn defocus_costs_contrast_over_the_process_window() {
        let a = image_layout(&lines(), 1, &optics(0.0), 2500.0, 500.0, 256, 0.35).unwrap();
        let b = image_layout(&lines(), 1, &optics(150.0), 2500.0, 500.0, 256, 0.35).unwrap();
        assert!(b.contrast < a.contrast, "defocus must cost contrast: {:.3} -> {:.3}", a.contrast, b.contrast);
    }

    #[test]
    fn an_empty_layer_images_to_dark() {
        let empty = Layout { units_per_micron: 1000.0, structures: vec![] };
        let r = image_layout(&empty, 1, &optics(0.0), 2500.0, 500.0, 128, 0.35).unwrap();
        assert!(r.intensity.iter().all(|v| v.abs() < 1e-9), "no polygons means no light");
        assert!(r.contrast.abs() < 1e-9);
        assert!(r.printed_nm.abs() < 1e-9);
    }

    #[test]
    fn a_clear_field_images_to_one() {
        let clear = Layout {
            units_per_micron: 1000.0,
            structures: vec![Structure {
                name: "CLEAR".into(),
                polygons: vec![Polygon {
                    layer: 1,
                    datatype: 0,
                    points: vec![(0.0, 0.0), (3.0, 0.0), (3.0, 2.0), (0.0, 2.0), (0.0, 0.0)],
                }],
            }],
        };
        let r = image_layout(&clear, 1, &optics(0.0), 2500.0, 500.0, 128, 0.35).unwrap();
        let mean = r.intensity.iter().sum::<f64>() / 128.0;
        assert!((mean - 1.0).abs() < 0.02, "a clear field must read 1.0, got {mean:.4}");
    }

    #[test]
    fn a_lower_threshold_prints_more() {
        let a = image_layout(&lines(), 1, &optics(0.0), 2500.0, 500.0, 256, 0.45).unwrap();
        let b = image_layout(&lines(), 1, &optics(0.0), 2500.0, 500.0, 256, 0.25).unwrap();
        assert!(b.printed_nm > a.printed_nm, "a lower threshold prints wider: {:.0} -> {:.0}", a.printed_nm, b.printed_nm);
    }

    #[test]
    fn a_layout_on_another_layer_images_to_dark() {
        let r = image_layout(&lines(), 5, &optics(0.0), 2500.0, 500.0, 128, 0.35).unwrap();
        assert!(r.contrast.abs() < 1e-9, "layer 5 is empty");
    }
}

#[cfg(test)]
mod tests {
    use atto_units::{MTorr, Mm, Watts};
    use super::*;

    fn optics(defocus: f64, sigma: f64) -> Optics {
        Optics { lambda_nm: 193.0, na: 0.93, sigma, defocus_nm: defocus }
    }

    /// 250 nm pitch: k1 = 0.60, comfortably inside the lambda/NA = 207 nm limit
    /// for this NA. A denser pitch has no first order in the pupil and images to
    /// a flat field, which is correct and is tested separately.
    fn grating() -> Grating {
        Grating { pitch_nm: 250.0, line_nm: 125.0 }
    }

    #[test]
    fn a_coherent_image_oscillates_with_defocus() {
        // Defocus is quadratic in pupil coordinate, so the +1 and -1 orders get
        // the *same* phase and the contrast oscillates rather than decaying: it
        // recovers at 400 nm after falling at 200 nm. A model that decays
        // monotonically here would be wrong.
        let g = grating();
        let a = optics(0.0, 0.0).contrast(&g).unwrap();
        let b = optics(200.0, 0.0).contrast(&g).unwrap();
        let c = optics(400.0, 0.0).contrast(&g).unwrap();
        assert!(a > 0.95, "best focus is near-perfect: {a:.3}");
        assert!(b < a, "contrast falls at 200 nm: {a:.3} -> {b:.3}");
        assert!(c > b, "and recovers at 400 nm: {b:.3} -> {c:.3}");
    }

    #[test]
    fn partial_coherence_collapses_monotonically() {
        // Averaging over source points washes the oscillation out, which is what
        // makes focus a process window at all.
        let g = grating();
        let mut prev = f64::INFINITY;
        for z in [0.0, 100.0, 200.0, 400.0] {
            let c = optics(z, 0.6).contrast(&g).unwrap();
            assert!(c < prev + 1e-9, "contrast rose at {z} nm: {prev:.4} -> {c:.4}");
            prev = c;
        }
        let a = optics(0.0, 0.6).contrast(&g).unwrap();
        let d = optics(400.0, 0.6).contrast(&g).unwrap();
        assert!(d < 0.5 * a, "a quarter of the contrast is left at 400 nm: {a:.3} -> {d:.3}");
    }

    #[test]
    fn the_bossung_curve_is_symmetric_about_best_focus() {
        let g = grating();
        for z in [50.0, 150.0, 300.0] {
            let a = optics(z, 0.6).contrast(&g).unwrap();
            let b = optics(-z, 0.6).contrast(&g).unwrap();
            assert!((a - b).abs() < 1e-9, "focus must be symmetric: {a:.5} vs {b:.5} at {z} nm");
        }
    }

    #[test]
    fn the_printed_cd_grows_with_dose() {
        let g = grating();
        let o = optics(0.0, 0.6);
        let img = o.image(&g, 256).unwrap();
        let a = cd_from_threshold(&img, Nm::new(g.pitch_nm), 0.45, true);
        let b = cd_from_threshold(&img, Nm::new(g.pitch_nm), 0.30, true);
        assert!(b > a, "a lower threshold prints a wider line: {a:.1} -> {b:.1} nm");
        assert!(a > 0.0 && a < g.pitch_nm);
    }

    #[test]
    fn the_image_is_normalised_and_physical() {
        let g = grating();
        let img = optics(0.0, 0.6).image(&g, 128).unwrap();
        assert_eq!(img.len(), 128);
        assert!(img.iter().all(|v| *v >= -1e-12 && v.is_finite()));
        let mean = img.iter().sum::<f64>() / img.len() as f64;
        assert!(mean > 0.28 && mean < 0.85, "mean intensity {mean:.3}");
    }

    #[test]
    fn a_pitch_below_the_resolution_limit_images_to_nothing() {
        // pitch < lambda/NA: no first order enters the pupil, so no modulation
        let g = Grating { pitch_nm: 100.0, line_nm: 50.0 };
        let c = optics(0.0, 0.0).contrast(&g).unwrap();
        assert!(c < 1e-9, "an unresolvable pitch must have no contrast: {c:.2e}");
    }

    #[test]
    fn bad_input_is_rejected() {
        let g = grating();
        assert!(Optics { lambda_nm: 0.0, ..optics(0.0, 0.6) }.image(&g, 64).is_err());
        assert!(Optics { na: 1.5, ..optics(0.0, 0.6) }.image(&g, 64).is_err());
        assert!(Optics { sigma: 2.0, ..optics(0.0, 0.6) }.image(&g, 64).is_err());
        assert!(optics(0.0, 0.6).image(&g, 4).is_err());
    }
}
