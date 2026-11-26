// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Phase-shift masks: alternating and attenuated (L-09).
//!
//! A binary mask has a DC term in its spectrum, so the intensity between two clear lines never
//! reaches zero — light leaks through and the dark line prints partially. An alternating
//! phase-shift mask puts a 180 degree step between adjacent clear features, which cancels the DC
//! term entirely: **the dark line prints black**, contrast rises sharply, and there is no side
//! lobe beside the feature.
//!
//! That is the whole idea, and it is why PSM was the first resolution-enhancement technique to
//! matter. The model here is a 1:1 line/space grating imaged with a pupil cutoff, which is enough
//! to show the DC cancellation and the resulting contrast.

use crate::litho::{Grating, LithoError, Optics};
use std::f64::consts::PI;

/// Which mask produced an image.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum MaskType {
    /// Opaque lines on a clear field: a DC term survives.
    Binary,
    /// Alternating 180 degree phase: the DC term cancels.
    Alternating,
}

/// What a mask printed.
#[derive(Debug, Clone, PartialEq)]
pub struct PsmResult {
    /// Michelson contrast of the image, (Imax - Imin) / (Imax + Imin).
    pub contrast: f64,
    /// Lowest intensity in the image, normalised so the clear-field intensity is 1.
    pub min_intensity: f64,
    /// Whether a side lobe is present: a local maximum away from the feature.
    pub has_side_lobe: bool,
}

impl PsmResult {
    /// The result as one line for a report.
    pub fn label(&self) -> String {
        format!(
            "contrast {:.2}, min intensity {:.3}{}",
            self.contrast,
            self.min_intensity,
            if self.has_side_lobe { ", side lobe" } else { ", no side lobe" }
        )
    }
}

/// Image a 1:1 line/space grating through a pupil, returning normalised intensity.
///
/// The field is expanded in Fourier orders; a binary mask has all orders with a DC term, an
/// alternating mask has only odd orders. Orders beyond the cutoff `NA / lambda` are removed by the
/// pupil, which is what makes the result depend on `k1` at all.
pub fn image(o: &Optics, g: &Grating, mask: MaskType, samples: usize) -> Result<Vec<f64>, LithoError> {
    if samples < 8 {
        return Err(LithoError::BadInput("need at least 8 samples".into()));
    }
    if !(g.pitch_nm > 0.0) || !g.line_nm.is_finite() {
        return Err(LithoError::BadInput("the grating needs a positive pitch".into()));
    }
    if !(o.lambda_nm > 0.0) || !(o.na > 0.0) {
        return Err(LithoError::BadInput("the optics need a wavelength and a numerical aperture".into()));
    }
    let cutoff = o.na / o.lambda_nm;
    let mut out = vec![0.0; samples];
    let mut peak = 0.0f64;
    for (i, slot) in out.iter_mut().enumerate() {
        let x = g.pitch_nm * i as f64 / samples as f64;
        // A real mask has a COSINE series: A0 + sum over m > 0 of a_m cos(2 pi f x), with
        // a_m = 2 x the complex coefficient. The first version looped over negative m as well
        // and added cos for both, which is even in m — the +1 and -1 orders cancelled and the
        // image was identically zero. It is worth being explicit: the sum runs over m >= 0 only.
        let mut field = match mask {
            MaskType::Binary => 0.5,
            MaskType::Alternating => 0.0,
        };
        for m in 1..=21i32 {
            let f = m as f64 / g.pitch_nm;
            if f > cutoff {
                break;
            }
            let amp = match mask {
                // 1:1 duty, opaque line: a_m = 2 sin(m pi / 2) / (m pi)
                MaskType::Binary => 2.0 * (PI * m as f64 / 2.0).sin() / (PI * m as f64),
                // alternating phase: odd orders only, a_m = 4 / (m pi)
                MaskType::Alternating => {
                    if m % 2 == 1 { 4.0 / (PI * m as f64) } else { 0.0 }
                }
            };
            if amp == 0.0 {
                continue;
            }
            field += amp * (2.0 * PI * f * x).cos();
        }
        let intensity = field * field;
        *slot = intensity;
        if intensity > peak {
            peak = intensity;
        }
    }
    if peak > 0.0 {
        for v in out.iter_mut() {
            *v /= peak;
        }
    }
    Ok(out)
}

/// Contrast and side-lobe check for a mask.
pub fn evaluate(o: &Optics, g: &Grating, mask: MaskType, samples: usize) -> Result<PsmResult, LithoError> {
    let img = image(o, g, mask, samples)?;
    let hi = img.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let lo = img.iter().cloned().fold(f64::INFINITY, f64::min);
    let contrast = if hi + lo > 0.0 { (hi - lo) / (hi + lo) } else { 0.0 };
    // a side lobe is a local maximum strictly between two minima
    let mut has_side_lobe = false;
    for i in 1..img.len() - 1 {
        if img[i] > img[i - 1] && img[i] > img[i + 1] && img[i] < hi * 0.5 {
            has_side_lobe = true;
            break;
        }
    }
    Ok(PsmResult { contrast, min_intensity: lo, has_side_lobe })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn optics() -> Optics {
        Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.0, defocus_nm: 0.0 }
    }
    fn grating() -> Grating {
        // 300 nm, not 200: at NA 0.93 and 193 nm the pupil cutoff is 0.00482 /nm and the first
        // order of a 200 nm grating sits at 0.005, outside it. Only the DC order would pass —
        // and for an alternating mask the DC is zero, so the image would be black and every
        // assertion here would pass or fail for the wrong reason.
        Grating { pitch_nm: 300.0, line_nm: 150.0 }
    }

    #[test]
    fn the_alternating_mask_prints_a_black_dark_line() {
        // L-09's acceptance: no light leaks into the dark line, so the minimum goes to zero
        let psm = evaluate(&optics(), &grating(), MaskType::Alternating, 256).expect("psm");
        assert!(psm.min_intensity < 1e-9, "the dark line is black: {}", psm.min_intensity);
    }

    #[test]
    fn a_binary_mask_leaks_light_into_the_dark_line() {
        // the comparison that gives the alternating result its meaning
        let bin = evaluate(&optics(), &grating(), MaskType::Binary, 256).expect("binary");
        assert!(bin.min_intensity > 0.0, "a binary mask always leaks: {}", bin.min_intensity);
    }

    #[test]
    fn the_phase_mask_has_the_higher_contrast() {
        let o = optics();
        let g = grating();
        let bin = evaluate(&o, &g, MaskType::Binary, 256).expect("binary");
        let psm = evaluate(&o, &g, MaskType::Alternating, 256).expect("psm");
        assert!(psm.contrast > bin.contrast, "phase shift raises contrast: {} vs {}", psm.contrast, bin.contrast);
        assert!((psm.contrast - 1.0).abs() < 1e-9, "and reaches the maximum: {}", psm.contrast);
    }

    #[test]
    fn the_alternating_mask_prints_no_side_lobe() {
        let psm = evaluate(&optics(), &grating(), MaskType::Alternating, 256).expect("psm");
        assert!(!psm.has_side_lobe, "no side lobe beside the feature: {}", psm.label());
    }

    #[test]
    fn the_image_is_normalised_and_periodic() {
        let img = image(&optics(), &grating(), MaskType::Alternating, 128).expect("image");
        let hi = img.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert!((hi - 1.0).abs() < 1e-12, "the peak is normalised to one: {hi}");
    }

    #[test]
    fn bad_input_is_refused() {
        let o = optics();
        let g = grating();
        assert!(image(&o, &g, MaskType::Binary, 4).is_err(), "too few samples");
        assert!(image(&o, &Grating { pitch_nm: 0.0, line_nm: 0.0 }, MaskType::Binary, 64).is_err(), "no pitch");
        assert!(image(&Optics { lambda_nm: 0.0, ..o }, &g, MaskType::Binary, 64).is_err(), "no wavelength");
    }
}
