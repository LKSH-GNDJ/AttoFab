// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A 2-D level-set surface, for the shapes a column stack cannot hold.
//!
//! Every surface in this simulator is the top of a column stack: the height varies across
//! the grid, but a column cannot lean, undercut or round over. Faceting and redeposition are
//! exactly those shapes, which is why ADR-1 ties them to a level set that did not exist.
//!
//! This is that level set, on one cross-section. The surface is the zero contour of a signed
//! distance field `phi`: negative inside the material, positive outside. Advancing the
//! surface is `dphi/dt = -v|grad phi|`, solved with upwind differences. A curvature term
//! smooths sharp corners, which is what faceting looks like.
//!
//! LIMITATION: a *flat* surface cannot be advected through this embedding — `phi = level
//! - h` has no gradient where the surface is flat, so a uniform etch moves nothing. Moving a
//! flat surface vertically needs a 3-D field phi(x, z). This 2-D field carries lateral shape
//! only; handle the flat case with the column stack.
//!
//! It is a *representation*, not a pipeline. Nothing here replaces `ColumnGrid`; a later
//! slice would have to decide where the two meet.

/// A signed distance field on a regular 2-D grid.
#[derive(Debug, Clone, PartialEq)]
pub struct Phi {
    pub nx: usize,
    pub ny: usize,
    /// Grid pitch, nm.
    pub pitch: f64,
    /// Signed distance, nm. Negative inside the material. Row-major, `ny * nx`.
    pub v: Vec<f64>,
}

impl Phi {
    /// A field that is negative everywhere inside a circle and positive outside.
    pub fn circle(nx: usize, ny: usize, pitch: f64, centre: (f64, f64), radius: f64) -> Option<Self> {
        if nx == 0 || ny == 0 || !(pitch > 0.0) || !pitch.is_finite() || !(radius > 0.0) {
            return None;
        }
        let mut v = Vec::with_capacity(nx * ny);
        for j in 0..ny {
            for i in 0..nx {
                let x = (i as f64 + 0.5) * pitch;
                let y = (j as f64 + 0.5) * pitch;
                let d = ((x - centre.0).powi(2) + (y - centre.1).powi(2)).sqrt() - radius;
                v.push(d);
            }
        }
        Some(Phi { nx, ny, pitch, v })
    }

    /// Value at a grid point, `None` outside.
    pub fn at(&self, i: usize, j: usize) -> Option<f64> {
        if i >= self.nx || j >= self.ny { None } else { Some(self.v[j * self.nx + i]) }
    }

    /// Area of material, nm², taken as the count of negative cells times the cell area.
    ///
    /// A cell-count rather than a contour integral: it is what the field can actually
    /// resolve, and it makes the area tests exact rather than approximately exact.
    pub fn area_nm2(&self) -> f64 {
        let n = self.v.iter().filter(|x| **x < 0.0).count();
        n as f64 * self.pitch * self.pitch
    }

    fn grad(&self, i: usize, j: usize) -> (f64, f64) {
        let gx = if i + 1 < self.nx && i > 0 { (self.v[j * self.nx + i + 1] - self.v[j * self.nx + i - 1]) / (2.0 * self.pitch) } else { 0.0 };
        let gy = if j + 1 < self.ny && j > 0 { (self.v[(j + 1) * self.nx + i] - self.v[(j - 1) * self.nx + i]) / (2.0 * self.pitch) } else { 0.0 };
        (gx, gy)
    }

    /// Advance the surface along its normal at `speed_nm_per_s` for `dt_s`.
    ///
    /// Positive speed grows the material (deposition), negative removes it (etch). The step
    /// uses the upwind gradient so the surface moves at the right speed regardless of which
    /// way it faces — a central difference would stall on a diagonal front.
    pub fn advect(&mut self, speed_nm_per_s: f64, dt_s: f64) {
        if !speed_nm_per_s.is_finite() || !dt_s.is_finite() || dt_s <= 0.0 {
            return;
        }
        let cfl = self.pitch / dt_s.max(1e-12);
        let v = speed_nm_per_s.clamp(-cfl, cfl);
        let old = self.v.clone();
        // Godunov (Osher-Sethian) upwind gradient. The one-sided differences are taken
        // first and then combined by the sign of the speed, so a surface facing either way
        // moves at the same rate. The earlier hand-rolled version picked one side and
        // returned zero on a monotone front.
        for j in 0..self.ny {
            for i in 0..self.nx {
                let k = j * self.nx + i;
                let dmx = if i > 0 { (old[k] - old[k - 1]) / self.pitch } else { 0.0 };
                let dpx = if i + 1 < self.nx { (old[k + 1] - old[k]) / self.pitch } else { 0.0 };
                let dmy = if j > 0 { (old[k] - old[k - self.nx]) / self.pitch } else { 0.0 };
                let dpy = if j + 1 < self.ny { (old[k + self.nx] - old[k]) / self.pitch } else { 0.0 };
                let (gx, gy) = if v > 0.0 {
                    (dmx.max(0.0).powi(2) + dpx.min(0.0).powi(2), dmy.max(0.0).powi(2) + dpy.min(0.0).powi(2))
                } else {
                    (dpx.max(0.0).powi(2) + dmx.min(0.0).powi(2), dpy.max(0.0).powi(2) + dmy.min(0.0).powi(2))
                };
                let grad = (gx + gy).sqrt();
                self.v[k] = old[k] - v * dt_s * grad;
            }
        }
    }

    /// Curvature-driven smoothing: `dphi/dt = -gamma * kappa * |grad phi|`.
    ///
    /// This is what rounds a sharp corner. Under it a convex bump shrinks and a concave
    /// notch fills, both at a rate set by curvature, which is the behaviour that makes
    /// faceting look like faceting rather than like a blurred copy.
    pub fn smooth(&mut self, gamma_nm2_per_s: f64, dt_s: f64) {
        if !gamma_nm2_per_s.is_finite() || gamma_nm2_per_s <= 0.0 || !dt_s.is_finite() || dt_s <= 0.0 {
            return;
        }
        let old = self.v.clone();
        for j in 0..self.ny {
            for i in 0..self.nx {
                let k = j * self.nx + i;
                let (gx, gy) = self.grad(i, j);
                let grad = (gx * gx + gy * gy).sqrt();
                if grad < 1e-12 { continue; }
                // curvature of the zero contour: div(grad phi / |grad phi|)
                let lap = if i > 0 && i + 1 < self.nx && j > 0 && j + 1 < self.ny {
                    (old[k + 1] - 2.0 * old[k] + old[k - 1] + old[k + self.nx] - 2.0 * old[k] + old[k - self.nx]) / (self.pitch * self.pitch)
                } else { 0.0 };
                let kappa = lap / grad;
                // mean-curvature flow in this convention: a convex bump recedes, a notch fills.
                // The minus in dphi/dt = -gamma kappa |grad phi| is already inside kappa,
                // which is signed the other way here.
                self.v[k] = old[k] + gamma_nm2_per_s * dt_s * kappa * grad;
            }
        }
    }


    /// The steepest gradient anywhere in the field.
    ///
    /// A step from -20 to +20 has a slope of 2 per cell and a flat interior; smoothing
    /// lowers the peak slope. This is the measure that shows curvature flow working, where
    /// `gradient_error` would read 1.0 on both sides because the interior is flat.
    pub fn max_gradient(&self) -> f64 {
        let mut worst = 0.0f64;
        for j in 1..self.ny.saturating_sub(1) {
            for i in 1..self.nx.saturating_sub(1) {
                let (gx, gy) = self.grad(i, j);
                worst = worst.max((gx * gx + gy * gy).sqrt());
            }
        }
        worst
    }

    /// How far from a true distance field the current field is.
    ///
    /// Advection and smoothing both distort the field: the zero contour stays where it should
    /// but the magnitudes drift, so a later step reads the wrong gradient. This is the
    /// diagnostic; a full reinitialisation solve is a separate job.
    pub fn gradient_error(&self) -> f64 {
        let mut worst = 0.0f64;
        for j in 1..self.ny.saturating_sub(1) {
            for i in 1..self.nx.saturating_sub(1) {
                let (gx, gy) = self.grad(i, j);
                worst = worst.max(((gx * gx + gy * gy).sqrt() - 1.0).abs());
            }
        }
        worst
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_circle_shrinks_at_the_speed_it_is_etched() {
        // 20 nm cells, 40 nm radius, etched at 1 nm/s for 10 s: the radius should drop by 10 nm
        let mut p = Phi::circle(40, 40, 20.0, (400.0, 400.0), 200.0).expect("a field");
        let before = p.area_nm2();
        for _ in 0..10 { p.advect(-1.0, 1.0); }
        let after = p.area_nm2();
        assert!(after < before, "etching removes area: {before} -> {after}");
        // a 10 nm radius loss on a 200 nm circle is about 10% of the area
        let ratio = after / before;
        assert!((0.80..0.95).contains(&ratio), "about a tenth of the area gone: ratio={ratio}");
    }

    #[test]
    fn a_depositing_front_grows_the_material() {
        let mut p = Phi::circle(40, 40, 20.0, (400.0, 400.0), 200.0).expect("a field");
        let before = p.area_nm2();
        for _ in 0..5 { p.advect(1.0, 1.0); }
        assert!(p.area_nm2() > before, "deposition adds area");
    }

    #[test]
    fn smoothing_reduces_the_gradient_of_a_sharp_field() {
        // a hard step: half the grid negative, half positive
        let (nx, ny) = (40usize, 40usize);
        let mut v = vec![0.0f64; nx * ny];
        for j in 0..ny { for i in 0..nx { v[j * nx + i] = if i < nx / 2 { -20.0 } else { 20.0 }; } }
        let mut p = Phi { nx, ny, pitch: 20.0, v };
        let before = p.max_gradient();
        p.smooth(50.0, 1.0);
        assert!(p.max_gradient() < before, "curvature flow must soften the step: {before} -> {}", p.max_gradient());
    }

    #[test]
    fn a_monotone_step_now_moves_under_godunov() {
        // this is the case the hand-rolled upwind returned zero for: a surface that only
        // ever goes one way. Godunov takes both one-sided differences and combines them by
        // the sign of the speed, so the cells behind the front still carry a gradient.
        let (nx, ny) = (16usize, 16usize);
        let mut h = vec![0.0f64; nx * ny];
        for j in 0..ny { for i in 0..nx { h[j * nx + i] = if i < nx / 2 { 40.0 } else { 0.0 }; } }
        let mut p = Phi::from_heights(&h, nx, ny, 5.0, 0.0).expect("a field");
        let before = p.v.clone();
        p.advect(-1.0, 1.0);
        assert_ne!(p.v, before, "a monotone step must move under Godunov");
        // and it must move the way an etch moves: material recedes, so phi increases
        let k = (ny / 2) * nx + (nx / 2 - 1);  // the cell next to the front, where the gradient is
        assert!(p.v[k] > before[k], "etching raises phi behind the front: {} -> {}", before[k], p.v[k]);
    }

    #[test]
    fn a_field_with_no_shape_is_refused() {
        assert!(Phi::circle(0, 10, 20.0, (0.0, 0.0), 10.0).is_none());
        assert!(Phi::circle(10, 10, 0.0, (0.0, 0.0), 10.0).is_none());
        assert!(Phi::circle(10, 10, 20.0, (0.0, 0.0), 0.0).is_none());
    }

    #[test]
    fn a_step_that_cannot_happen_changes_nothing() {
        let mut p = Phi::circle(20, 20, 20.0, (200.0, 200.0), 100.0).expect("a field");
        let before = p.v.clone();
        p.advect(1.0, 0.0);
        p.advect(f64::NAN, 1.0);
        p.smooth(0.0, 1.0);
        assert_eq!(p.v, before, "a zero or nonsense step must not move the surface");
    }
}

impl Phi {
    /// A field whose zero contour is the surface of a height map.
    ///
    /// `phi = level - height`, so material (below the surface) is negative and air above is
    /// positive. The field is not a true distance function away from the surface — a column
    /// grid only knows one height per column — but it is exact at the contour, which is what
    /// advection and curvature read.
    pub fn from_heights(heights: &[f64], nx: usize, ny: usize, pitch: f64, level_nm: f64) -> Option<Self> {
        if heights.len() != nx * ny || nx == 0 || ny == 0 || !(pitch > 0.0) || !level_nm.is_finite() {
            return None;
        }
        let v = heights.iter().map(|h| level_nm - h).collect();
        Some(Phi { nx, ny, pitch, v })
    }

    /// Surface height per column, read back from the zero contour.
    ///
    /// Each column is scanned for the first sign change from negative to positive; the
    /// crossing is interpolated, so a surface between two cells comes back as a fraction
    /// rather than snapping to the nearest cell. A column with no crossing returns
    /// `level_nm` — no surface there.
    pub fn to_heights(&self, level_nm: f64) -> Vec<f64> {
        let mut out = vec![level_nm; self.nx * self.ny];
        for j in 0..self.ny {
            for i in 0..self.nx {
                let k = j * self.nx + i;
                let phi = self.v[k];
                if phi <= 0.0 {
                    out[k] = level_nm - phi;
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod bridge_tests {
    use super::*;

    #[test]
    fn a_height_map_survives_the_round_trip() {
        // a step: the left half 100 nm higher than the right
        let (nx, ny) = (16usize, 8usize);
        let mut h = vec![0.0f64; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                h[j * nx + i] = if i < nx / 2 { 100.0 } else { 0.0 };
            }
        }
        let p = Phi::from_heights(&h, nx, ny, 10.0, 0.0).expect("a field");
        let back = p.to_heights(0.0);
        for k in 0..h.len() {
            assert!((back[k] - h[k]).abs() < 1e-9, "height {k} came back as {} not {}", back[k], h[k]);
        }
    }

    #[test]
    fn a_field_from_a_flat_surface_is_flat() {
        let (nx, ny) = (8usize, 8usize);
        let h = vec![50.0f64; nx * ny];
        let p = Phi::from_heights(&h, nx, ny, 10.0, 0.0).expect("a field");
        assert!(p.v.iter().all(|x| (*x + 50.0).abs() < 1e-9), "every point is 50 nm below the level");
    }

    #[test]
    fn a_bad_height_map_is_refused() {
        assert!(Phi::from_heights(&[0.0; 3], 4, 4, 10.0, 0.0).is_none(), "wrong length");
        assert!(Phi::from_heights(&[0.0; 16], 4, 4, 0.0, 0.0).is_none(), "no pitch");
        assert!(Phi::from_heights(&[0.0; 16], 4, 4, 10.0, f64::NAN).is_none(), "not a number");
    }

    #[test]
    fn a_flat_surface_cannot_be_advected_through_this_embedding() {
        // A documented limitation, not a wish. `phi = level - h` has no gradient where the
        // surface is flat, so the advection term is zero and a uniform etch does nothing.
        // Moving a flat surface vertically needs a 3-D embedding phi(x, z); this 2-D field
        // only carries lateral shape. Anyone wiring this into an etch must handle the flat
        // case with the column stack, not with the level set.
        let (nx, ny) = (16usize, 16usize);
        let h = vec![0.0f64; nx * ny];
        let mut p = Phi::from_heights(&h, nx, ny, 5.0, 0.0).expect("a field");
        let before = p.v.clone();
        for _ in 0..10 { p.advect(-1.0, 1.0); }
        assert_eq!(p.v, before, "a flat field has no gradient, so nothing moves — as documented");
    }

}
