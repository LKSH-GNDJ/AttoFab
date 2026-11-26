// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Line-of-sight flux: what a surface can actually see of the arriving species.
//!
//! Step coverage is not a single number. In a physical vapour deposition the source is far
//! away and arrives at an angle, so a surface in a trench bottom sees almost none of it
//! while the field beside the trench sees all of it. That is why a PVD film is thin at the
//! bottom of a hole and thick on the flat — and why the mouth pinches shut.
//!
//! The existing `step_coverage` parameter states the answer as one number for the whole
//! grid. This computes it from the geometry instead: for each column, trace the ray toward
//! the source and ask whether any nearer column blocks it.
//!
//! It is a 2-D cross-section, not a 3-D field solve: no re-emission, no surface diffusion,
//! no angular distribution beyond the spread given. What it does get right is shadowing,
//! which is the dominant term.
//!
//! For the three-dimensional solve with sticking and re-emission, see `field3d`.

/// The direction the arriving species comes from, as an angle from the surface normal.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Source {
    /// Tilt from vertical, degrees. 0 is straight down, 90 is fully grazing.
    pub angle_deg: f64,
    /// Angular spread of the source, degrees. 0 is a single direction.
    pub spread_deg: f64,
}

impl Default for Source {
    fn default() -> Self {
        Source { angle_deg: 0.0, spread_deg: 10.0 }
    }
}

impl Source {
    /// Refuse a source whose angle or spread is outside its domain.
    pub fn validate(&self) -> Result<(), String> {
        if !self.angle_deg.is_finite() || self.angle_deg.abs() >= 90.0 {
            return Err("source angle must be finite and within (-90, 90) degrees".into());
        }
        if !self.spread_deg.is_finite() || !(0.0..=90.0).contains(&self.spread_deg) {
            return Err("source spread must be finite and within [0, 90] degrees".into());
        }
        Ok(())
    }
}

/// Is the ray from column `i` toward `angle_deg` blocked by a nearer column?
///
/// The ray leaves the surface at `(x_i, h_i)` heading toward the source. A column `j` blocks
/// it when the surface at `j` rises above the ray's height at that x. This is the whole of
/// shadowing: no other term is needed to explain why a trench bottom stays bare.
pub fn is_shadowed(heights: &[f64], pitch: f64, i: usize, angle_deg: f64) -> bool {
    if i >= heights.len() || !(pitch > 0.0) || !angle_deg.is_finite() {
        return false;
    }
    let theta = angle_deg.to_radians();
    let (sx, cz) = (theta.sin(), theta.cos());
    if cz <= 1e-9 {
        return true; // fully grazing: nothing arrives
    }
    let h0 = heights[i];
    // the source is above; walk away from the column in the direction the light comes from
    let dir: isize = if sx >= 0.0 { -1 } else { 1 };
    let steps = heights.len() as isize;
    let mut k = i as isize + dir;
    let mut n = 1isize;
    while n <= steps && k >= 0 && k < heights.len() as isize {
        let dx = n as f64 * pitch;
        let ray_z = h0 + dx * cz / sx.abs().max(1e-9);
        if heights[k as usize] > ray_z {
            return true;
        }
        k += dir;
        n += 1;
    }
    false
}

/// Fraction of the source each column can see, 0..1, and the cosine of incidence.
///
/// Two effects multiply. Shadowing is binary per direction: blocked or not. Incidence is
/// `cos(angle)`: a tilted source spreads the same flux over a larger area. Averaging over
/// the source spread gives a smooth profile instead of a hard edge.
pub fn visible_fraction(heights: &[f64], pitch: f64, src: &Source) -> Vec<f64> {
    let n = heights.len();
    let mut out = vec![0.0f64; n];
    if n == 0 || !(pitch > 0.0) {
        return out;
    }
    let samples = 9usize;
    let spread = src.spread_deg.max(0.0);
    for i in 0..n {
        let mut acc = 0.0;
        let mut wsum = 0.0;
        for s in 0..samples {
            let frac = if samples == 1 { 0.5 } else { s as f64 / (samples - 1) as f64 - 0.5 };
            let a = src.angle_deg + frac * spread;
            if a.abs() >= 90.0 {
                continue;
            }
            let w = a.to_radians().cos().max(0.0);
            wsum += w;
            if !is_shadowed(heights, pitch, i, a) {
                acc += w;
            }
        }
        out[i] = if wsum > 0.0 { (acc / wsum).clamp(0.0, 1.0) } else { 0.0 };
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_flat_surface_sees_the_whole_source() {
        let h = vec![0.0f64; 32];
        let v = visible_fraction(&h, 10.0, &Source { angle_deg: 0.0, spread_deg: 0.0 });
        assert!(v.iter().all(|x| (*x - 1.0).abs() < 1e-9), "nothing blocks a flat field");
    }

    #[test]
    fn a_trench_bottom_is_shadowed_from_a_tilted_source() {
        // a 200 nm wide, 400 nm deep trench in a 10 nm grid
        let mut h = vec![0.0f64; 64];
        for i in 28..36 { h[i] = -400.0; }
        let v = visible_fraction(&h, 10.0, &Source { angle_deg: 60.0, spread_deg: 0.0 });
        let bottom = v[32];
        let field = v[5];
        assert!(field > 0.99, "the open field is fully exposed: {field}");
        assert!(bottom < 0.05, "the trench bottom sees almost nothing at 60 degrees: {bottom}");
    }

    #[test]
    fn a_vertical_source_reaches_the_trench_bottom() {
        // straight down: a vertical-walled trench is not shadowed by its own walls
        let mut h = vec![0.0f64; 64];
        for i in 28..36 { h[i] = -400.0; }
        let v = visible_fraction(&h, 10.0, &Source { angle_deg: 0.0, spread_deg: 0.0 });
        assert!(v[32] > 0.99, "a vertical source is not blocked: {}", v[32]);
    }

    #[test]
    fn a_grazing_source_is_blocked_by_any_rise() {
        let mut h = vec![0.0f64; 32];
        h[20] = 100.0; // a single step
        // the source comes from the right, so the step at column 20 casts shadow leftward
        let v = visible_fraction(&h, 10.0, &Source { angle_deg: -80.0, spread_deg: 0.0 });
        assert!(v[10] < 0.05, "the far side of a step is in shadow at 80 degrees: {}", v[10]);
    }

    #[test]
    fn a_source_that_cannot_exist_is_refused() {
        assert!(Source { angle_deg: 90.0, spread_deg: 0.0 }.validate().is_err(), "grazing is not a source");
        assert!(Source { angle_deg: 0.0, spread_deg: -1.0 }.validate().is_err(), "negative spread");
        assert!(Source { angle_deg: f64::NAN, spread_deg: 0.0 }.validate().is_err(), "not a number");
        assert!(Source::default().validate().is_ok());
    }
}
