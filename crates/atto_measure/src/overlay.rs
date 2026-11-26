// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Cross-section image comparison (M-10).
//!
//! A TEM or SEM cross-section is the ground truth a process simulator is judged against, and
//! comparing one with a simulation is harder than it looks: the two images have different
//! brightness, different contrast, and are registered to each other only approximately.
//!
//! This module compares **binarised boundaries**, not grey levels. A grey-level difference
//! between a TEM and a renderer is dominated by imaging, not by geometry, so the metric here
//! asks the only question a process engineer asks of a cross-section: are the material
//! boundaries in the same place?

/// A single-channel image, row-major, values in 0..1.
#[derive(Debug, Clone, PartialEq)]
pub struct Image {
    pub width: usize,
    pub height: usize,
    pub pixels: Vec<f64>,
}

/// Why an overlay could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub enum OverlayError {
    /// Dimensions disagree, or a buffer length does not match them.
    BadImage(&'static str),
    /// An image had no foreground at all, so there is nothing to overlay.
    Empty,
}

impl std::fmt::Display for OverlayError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OverlayError::BadImage(w) => write!(f, "unusable image: {w}"),
            OverlayError::Empty => write!(f, "an image has no foreground"),
        }
    }
}

/// How far apart two cross-sections are.
#[derive(Debug, Clone, PartialEq)]
pub struct Overlay {
    /// Intersection over union of the foreground masks, 0..1. One means identical.
    pub iou: f64,
    /// Mean distance from a simulated boundary pixel to the nearest measured one, in pixels.
    pub mean_boundary_error_px: f64,
    /// Largest such distance, pixels.
    pub max_boundary_error_px: f64,
    /// Best integer registration offset, (dx, dy) pixels.
    ///
    /// The measured image was sampled at `(x + dx, y + dy)`, so a positive `dx` means the
    /// measured content lay this far to the right of the simulation and was pulled back.
    pub offset_px: (i32, i32),
}

impl Image {
    fn check(&self) -> Result<(), OverlayError> {
        if self.width == 0 || self.height == 0 {
            return Err(OverlayError::BadImage("zero extent"));
        }
        if self.pixels.len() != self.width * self.height {
            return Err(OverlayError::BadImage("buffer length does not match the extent"));
        }
        if self.pixels.iter().any(|p| !p.is_finite()) {
            return Err(OverlayError::BadImage("a pixel is not finite"));
        }
        Ok(())
    }

    /// Foreground mask at a threshold: 1 where the pixel is at or above it.
    fn mask(&self, threshold: f64) -> Vec<bool> {
        self.pixels.iter().map(|p| *p >= threshold).collect()
    }

    fn shifted_mask(&self, threshold: f64, dx: i32, dy: i32) -> Vec<bool> {
        let m = self.mask(threshold);
        let mut out = vec![false; m.len()];
        for y in 0..self.height {
            for x in 0..self.width {
                let sx = x as i32 + dx;
                let sy = y as i32 + dy;
                if sx < 0 || sy < 0 || sx as usize >= self.width || sy as usize >= self.height {
                    continue;
                }
                out[y * self.width + x] = m[sy as usize * self.width + sx as usize];
            }
        }
        out
    }

    /// Pixels on the boundary of the mask: foreground with at least one background neighbour.
    fn boundary(&self, mask: &[bool]) -> Vec<(usize, usize)> {
        let mut out = Vec::new();
        for y in 0..self.height {
            for x in 0..self.width {
                if !mask[y * self.width + x] {
                    continue;
                }
                let mut edge = false;
                for (dx, dy) in [(-1i32, 0i32), (1, 0), (0, -1), (0, 1)] {
                    let (nx, ny) = (x as i32 + dx, y as i32 + dy);
                    if nx < 0 || ny < 0 || nx as usize >= self.width || ny as usize >= self.height {
                        edge = true;
                        continue;
                    }
                    if !mask[ny as usize * self.width + nx as usize] {
                        edge = true;
                    }
                }
                if edge {
                    out.push((x, y));
                }
            }
        }
        out
    }
}

fn iou_of(a: &[bool], b: &[bool]) -> f64 {
    let mut inter = 0usize;
    let mut union = 0usize;
    for i in 0..a.len() {
        match (a[i], b[i]) {
            (true, true) => {
                inter += 1;
                union += 1;
            }
            (false, false) => {}
            _ => union += 1,
        }
    }
    if union == 0 {
        return 1.0;
    }
    inter as f64 / union as f64
}

/// Compare a simulated cross-section with a measured one.
///
/// The registration search covers `±search_px` in both axes, which is the range a person would
/// try by eye. Everything is measured on binarised boundaries, so brightness and contrast
/// differences between a TEM and a renderer do not enter the number.
pub fn compare(
    simulated: &Image,
    measured: &Image,
    threshold: f64,
    search_px: i32,
) -> Result<Overlay, OverlayError> {
    simulated.check()?;
    measured.check()?;
    if simulated.width != measured.width || simulated.height != measured.height {
        return Err(OverlayError::BadImage("the two images differ in extent"));
    }
    if !threshold.is_finite() {
        return Err(OverlayError::BadImage("threshold is not a number"));
    }
    let search = search_px.max(0);

    let sim = simulated.mask(threshold);
    if !sim.iter().any(|b| *b) {
        return Err(OverlayError::Empty);
    }

    let mut best: Option<(f64, i32, i32)> = None;
    for dy in -search..=search {
        for dx in -search..=search {
            let meas = measured.shifted_mask(threshold, dx, dy);
            let v = iou_of(&sim, &meas);
            let better = match best {
                None => true,
                Some((b, _, _)) => v > b + 1e-12,
            };
            if better {
                best = Some((v, dx, dy));
            }
        }
    }
    let (iou, dx, dy) = best.ok_or(OverlayError::Empty)?;

    let meas = measured.shifted_mask(threshold, dx, dy);
    let sim_edge = simulated.boundary(&sim);
    let meas_edge = measured.boundary(&meas);
    if sim_edge.is_empty() || meas_edge.is_empty() {
        return Err(OverlayError::Empty);
    }

    // nearest measured boundary pixel for each simulated boundary pixel, brute force
    let mut sum = 0.0;
    let mut worst: f64 = 0.0;
    for (x, y) in &sim_edge {
        let mut nearest = f64::INFINITY;
        for (mx, my) in &meas_edge {
            let d = (((*x as f64) - (*mx as f64)).powi(2) + ((*y as f64) - (*my as f64)).powi(2)).sqrt();
            if d < nearest {
                nearest = d;
            }
        }
        sum += nearest;
        if nearest > worst {
            worst = nearest;
        }
    }
    let mean = sum / sim_edge.len() as f64;

    Ok(Overlay { iou, mean_boundary_error_px: mean, max_boundary_error_px: worst, offset_px: (dx, dy) })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A filled rectangle: a crude but unambiguous cross-section.
    fn rect(w: usize, h: usize, x0: usize, y0: usize, rw: usize, rh: usize) -> Image {
        let mut px = vec![0.0; w * h];
        for y in y0..(y0 + rh).min(h) {
            for x in x0..(x0 + rw).min(w) {
                px[y * w + x] = 1.0;
            }
        }
        Image { width: w, height: h, pixels: px }
    }

    #[test]
    fn an_identical_pair_has_zero_overlay_error() {
        let a = rect(40, 30, 10, 8, 16, 12);
        let o = compare(&a, &a, 0.5, 3).expect("an overlay");
        assert!((o.iou - 1.0).abs() < 1e-12, "identical masks overlap exactly: {}", o.iou);
        assert!(o.mean_boundary_error_px.abs() < 1e-12, "and no boundary is displaced: {}", o.mean_boundary_error_px);
        assert_eq!(o.offset_px, (0, 0), "and nothing needed registering");
    }

    #[test]
    fn a_shifted_pair_is_registered_and_reports_the_offset() {
        let a = rect(40, 30, 10, 8, 16, 12);
        let b = rect(40, 30, 13, 8, 16, 12); // 3 px to the right
        let o = compare(&a, &b, 0.5, 5).expect("an overlay");
        assert_eq!(o.offset_px, (3, 0), "the measurement lay 3 px right and was pulled back");
        assert!((o.iou - 1.0).abs() < 1e-12, "once registered they coincide: {}", o.iou);
    }

    #[test]
    fn a_genuinely_different_pair_reports_an_error() {
        // same extent, different shape: this must not register away
        let a = rect(40, 30, 10, 8, 16, 12);
        let b = rect(40, 30, 10, 8, 16, 6); // half the height
        let o = compare(&a, &b, 0.5, 2).expect("an overlay");
        assert!(o.iou < 0.7, "half the area is missing: {}", o.iou);
        assert!(o.mean_boundary_error_px > 1.0, "and the boundary is displaced: {}", o.mean_boundary_error_px);
    }

    #[test]
    fn brightness_and_contrast_do_not_change_the_answer() {
        // a dimmer measured image with the same geometry must compare identically
        let a = rect(40, 30, 10, 8, 16, 12);
        let dim = Image {
            width: a.width,
            height: a.height,
            pixels: a.pixels.iter().map(|p| p * 0.6).collect(),
        };
        let o = compare(&a, &dim, 0.3, 2).expect("an overlay");
        assert!((o.iou - 1.0).abs() < 1e-12, "a dimmer image of the same shape overlays exactly: {}", o.iou);
    }

    #[test]
    fn unusable_images_are_refused() {
        let a = rect(20, 20, 5, 5, 8, 8);
        let mut short = a.clone();
        short.pixels.truncate(10);
        assert!(compare(&short, &a, 0.5, 2).is_err(), "buffer too short");
        let empty = Image { width: 20, height: 20, pixels: vec![0.0; 400] };
        assert!(compare(&empty, &a, 0.5, 2).is_err(), "nothing to overlay");
        let b = rect(10, 10, 2, 2, 4, 4);
        assert!(compare(&a, &b, 0.5, 2).is_err(), "different extents");
    }
}
