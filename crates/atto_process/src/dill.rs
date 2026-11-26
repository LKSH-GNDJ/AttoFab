// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Exposure and bake: how much photoactive compound is destroyed, and how far it spreads.
//!
//! Two effects the threshold model could not express:
//!
//! **Dill absorption.** The resist absorbs its own exposure. Light entering the top is
//! attenuated as it goes down, so the bottom of a 150 nm film sees less dose than the top.
//! The compound fraction therefore varies with depth, and a develop that removes from the
//! top down meets a *different* rate at every level — which is how a resist profile gets a
//! foot and a rounded top.
//!
//! **Post-exposure bake.** After exposure the compound is still mobile. Baking lets it
//! diffuse sideways, which smooths the latent image. This is why a bake softens a sharp
//! aerial image into a gentler dose profile, and why leaving the bake out makes simulated
//! sidewalls unrealistically steep.
//!
//! The absorption coefficients are Dill's A (bleachable) and B (permanent), in 1/µm; C is
//! the exposure rate constant, per unit dose.

/// Dill exposure parameters for one resist.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DillModel {
    /// Bleachable absorption `A`, 1/µm. Falls as the compound is consumed.
    pub a_per_um: f64,
    /// Permanent absorption `B`, 1/µm. Does not bleach.
    pub b_per_um: f64,
    /// Exposure rate constant `C`, per unit dose.
    pub c: f64,
}

impl Default for DillModel {
    fn default() -> Self {
        // typical positive KrF resist; replace with measured values for your resist
        // `c` is per unit *dose ratio*, not per mJ/cm2: the simulator feeds a dimensionless
        // dose where 1.0 is the dose to clear. A real Dill C (order 0.03 /mJ) multiplied by
        // a dose ratio of 2 leaves 94% of the compound, and with a Mack n of 8 that rate is
        // 1e-10 of R_max — the resist never develops. The scale must match the input.
        DillModel { a_per_um: 0.6, b_per_um: 0.2, c: 2.5 }
    }
}

impl DillModel {
    /// Every parameter must be finite and non-negative.
    pub fn validate(&self) -> Result<(), String> {
        for (name, v) in [("A", self.a_per_um), ("B", self.b_per_um), ("C", self.c)] {
            if !v.is_finite() || v < 0.0 {
                return Err(format!("Dill {name} must be finite and >= 0"));
            }
        }
        Ok(())
    }

    /// Fraction of photoactive compound remaining at `depth_nm` below the resist top.
    ///
    /// 1.0 is untouched, 0.0 is fully bleached. The attenuation uses the *mean* compound
    /// fraction over the film, which is the standard single-pass approximation: solving the
    /// coupled depth profile exactly needs an integral no column simulator can afford.
    pub fn pac_fraction(&self, dose: f64, depth: atto_units::Nm, thickness: atto_units::Nm) -> f64 {
        let depth_nm = depth.0;
        let thickness_nm = thickness.0;
        if !dose.is_finite() || dose <= 0.0 || !depth_nm.is_finite() || depth_nm < 0.0 {
            return 1.0;
        }
        let thick = if thickness_nm.is_finite() && thickness_nm > 0.0 { thickness_nm } else { 1.0 };
        // the film as a whole bleaches first; that sets how much light reaches the bottom
        let bulk = 1.0 - (-self.c * dose).exp();
        let alpha_per_um = self.a_per_um * bulk + self.b_per_um;
        let z_um = (depth_nm / thick.min(depth_nm.max(thick))) * thick / 1000.0;
        let attenuation = (-alpha_per_um * z_um).exp();
        let local = 1.0 - (-self.c * dose * attenuation).exp();
        (1.0 - local).clamp(0.0, 1.0)
    }
}

/// Diffusion length of the compound during a bake, nm.
///
/// `sigma = sqrt(2 D t)`, the standard random-walk result.
pub fn bake_sigma_nm(diffusivity_nm2_per_min: f64, time_min: f64) -> f64 {
    if !diffusivity_nm2_per_min.is_finite() || diffusivity_nm2_per_min < 0.0 || !time_min.is_finite() || time_min <= 0.0 {
        return 0.0;
    }
    (2.0 * diffusivity_nm2_per_min * time_min).sqrt()
}

/// Smooth a scalar map with a Gaussian of the given width, nm.
///
/// This is the bake: the latent image is blurred, so a step becomes a ramp. Values outside
/// the map are clamped to the edge, which keeps the mean where the field is uniform.
pub fn bake_blur(values: &[f32], nx: usize, ny: usize, sigma_nm: f64, pitch_nm: f64) -> Vec<f32> {
    if values.len() != nx * ny || nx == 0 || ny == 0 {
        return values.to_vec();
    }
    if !sigma_nm.is_finite() || sigma_nm <= 0.0 || !pitch_nm.is_finite() || pitch_nm <= 0.0 {
        return values.to_vec();
    }
    let sigma = sigma_nm / pitch_nm;
    let radius = (3.0 * sigma).ceil() as isize;
    if radius <= 0 {
        return values.to_vec();
    }
    let mut kernel = Vec::with_capacity((2 * radius + 1) as usize);
    let mut sum = 0.0f64;
    for i in -radius..=radius {
        let w = (-(i as f64).powi(2) / (2.0 * sigma * sigma)).exp();
        kernel.push(w);
        sum += w;
    }
    for w in &mut kernel {
        *w /= sum;
    }
    let at = |x: isize, y: isize| -> f32 {
        let xx = x.clamp(0, nx as isize - 1) as usize;
        let yy = y.clamp(0, ny as isize - 1) as usize;
        values[yy * nx + xx]
    };
    let mut out = vec![0.0f32; values.len()];
    for j in 0..ny as isize {
        for i in 0..nx as isize {
            let mut acc = 0.0f64;
            for (k, w) in kernel.iter().enumerate() {
                acc += *w * at(i + k as isize - radius, j) as f64;
            }
            out[j as usize * nx + i as usize] = acc as f32;
        }
    }
    // second pass, vertical
    let horiz = out.clone();
    for j in 0..ny as isize {
        for i in 0..nx as isize {
            let mut acc = 0.0f64;
            for (k, w) in kernel.iter().enumerate() {
                let yy = (j + k as isize - radius).clamp(0, ny as isize - 1) as usize;
                acc += *w * horiz[yy * nx + i as usize] as f64;
            }
            out[j as usize * nx + i as usize] = acc as f32;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deeper_resist_is_less_exposed() {
        let d = DillModel::default();
        let top = d.pac_fraction(2.0, atto_units::Nm::new(0.0), atto_units::Nm::new(150.0));
        let bottom = d.pac_fraction(2.0, atto_units::Nm::new(150.0), atto_units::Nm::new(150.0));
        assert!(top < bottom, "the top bleaches more than the bottom: top={top} bottom={bottom}");
    }

    #[test]
    fn more_dose_bleaches_more() {
        let d = DillModel::default();
        assert!(d.pac_fraction(3.0, atto_units::Nm::new(0.0), atto_units::Nm::new(150.0)) < d.pac_fraction(0.5, atto_units::Nm::new(0.0), atto_units::Nm::new(150.0)));
    }

    #[test]
    fn a_zero_dose_leaves_the_compound_alone() {
        let d = DillModel::default();
        assert_eq!(d.pac_fraction(0.0, atto_units::Nm::new(50.0), atto_units::Nm::new(150.0)), 1.0);
    }

    #[test]
    fn the_bake_length_grows_as_the_square_root_of_time() {
        let one = bake_sigma_nm(100.0, 1.0);
        let four = bake_sigma_nm(100.0, 4.0);
        assert!((four / one - 2.0).abs() < 1e-9, "quadruple the time doubles the length: {one} {four}");
    }

    #[test]
    fn the_bake_turns_a_step_into_a_ramp() {
        // half the map at 1.0, half at 0.0
        // wide enough that "far from the edge" is farther than the kernel reaches:
        // sigma 30 nm at 10 nm pitch is 3 cells, so the kernel spans 9 cells
        let (nx, ny) = (64usize, 4usize);
        let mut v = vec![0.0f32; nx * ny];
        for j in 0..ny { for i in nx / 2..nx { v[j * nx + i] = 1.0; } }
        let out = bake_blur(&v, nx, ny, 30.0, 10.0);
        // a cell just left of the step must now hold some of the bright side
        let edge = out[0 * nx + (nx / 2 - 1)];
        assert!(edge > 0.0 && edge < 1.0, "the step is now a ramp: {edge}");
        assert_eq!(out[0], 0.0, "far from the edge it stays dark");
        assert_eq!(out[nx - 1], 1.0, "and far on the bright side it stays bright");
    }

    #[test]
    fn a_bake_that_cannot_happen_is_ignored_rather_than_invented() {
        assert_eq!(bake_sigma_nm(0.0, 10.0), 0.0);
        assert_eq!(bake_sigma_nm(100.0, 0.0), 0.0);
        let v = vec![1.0f32; 4];
        assert_eq!(bake_blur(&v, 2, 2, 0.0, 10.0), v, "no width, no change");
        assert_eq!(bake_blur(&v, 2, 2, 10.0, 0.0), v, "no pitch, no change");
    }

    #[test]
    fn a_model_that_cannot_be_used_is_refused() {
        let mut d = DillModel::default();
        d.a_per_um = -1.0;
        assert!(d.validate().is_err(), "negative A");
        let mut d = DillModel::default();
        d.c = f64::NAN;
        assert!(d.validate().is_err(), "not a number");
        assert!(DillModel::default().validate().is_ok());
    }
}

/// Bake diffusion length from a bake temperature and time, nm.
///
/// `D = D0 exp(-Ea / kT)`, then `sigma = sqrt(2 D t)`. The two constants are defaults for a
/// typical positive resist and are **not measured** — they are the first thing to replace
/// with your own bake data. What the model gets right without calibration is the shape:
/// hotter and longer both spread the latent image further.
pub fn bake_sigma_from_temperature(temp: atto_units::Kelvin, time: atto_units::Minutes) -> atto_units::Nm {
    let (temp_k, time_min) = (temp.0, time.0);
    if !temp_k.is_finite() || !time_min.is_finite() || time_min <= 0.0 {
        return atto_units::Nm::new(0.0);
    }
    let t_k = temp_k.max(1.0);
    // D0 in nm^2/min, Ea in eV
    // D0 chosen so that a 110 C bake for a minute gives about 10 nm of diffusion, which
    // is what a real post-exposure bake does. The first value tried (1e12) gave 9e-4 nm:
    // the bake was switched off, and a switched-off bake looks like a working one.
    let (d0, ea_ev) = (1.3e20_f64, 1.4_f64);
    let k_ev_per_k = 8.617_333_262e-5_f64;
    let d = d0 * (-ea_ev / (k_ev_per_k * t_k)).exp();
    atto_units::Nm::new(bake_sigma_nm(d, time_min))
}

/// Depth the developer reaches in one column, stepping through the film.
///
/// The rate depends on how much compound is left, and that depends on how deep the
/// developer already is, because the resist absorbs its own exposure. So the depth is found
/// by stepping: at each level, ask Dill what compound is left, ask Mack how fast that
/// develops, and advance. A single step would assume the whole film is equally exposed,
/// which is the thing this model exists to stop assuming.
pub fn develop_depth(
    dose_ratio: f64,
    thickness: atto_units::Nm,
    mack: &crate::resist::DevelopModel,
    dill: &DillModel,
    time: atto_units::Minutes,
    steps: usize,
) -> atto_units::Nm {
    let (time_min, thickness_nm) = (time.0, thickness.0);
    if !time_min.is_finite() || time_min <= 0.0 || !thickness_nm.is_finite() || thickness_nm <= 0.0 {
        return atto_units::Nm::new(0.0);
    }
    let n = steps.clamp(1, 256);
    let dt = time_min / n as f64;
    let mut depth = 0.0f64;
    for _ in 0..n {
        if depth >= thickness_nm {
            break;
        }
        let m = dill.pac_fraction(dose_ratio, atto_units::Nm::new(depth), atto_units::Nm::new(thickness_nm));
        let r = mack.rate_nm_per_min(m);
        depth = (depth + r * dt).min(thickness_nm);
    }
    atto_units::Nm::new(depth)
}

#[cfg(test)]
mod wiring_tests {
    use super::*;
    use crate::resist::DevelopModel;

    #[test]
    fn absorption_makes_the_developer_stop_sooner() {
        // same dose, same time: the film that absorbs its own light leaves less compound at
        // depth, so the developer slows and stops higher up
        let mack = DevelopModel::default();
        let absorbing = DillModel::default();
        let clear = DillModel { a_per_um: 0.0, b_per_um: 0.0, c: absorbing.c };
        let d_abs = develop_depth(2.0, atto_units::Nm::new(150.0), &mack, &absorbing, atto_units::Minutes::new(0.4), 32).0;
        let d_clear = develop_depth(2.0, atto_units::Nm::new(150.0), &mack, &clear, atto_units::Minutes::new(0.4), 32).0;
        assert!(d_abs <= d_clear + 1e-9, "absorption must not develop deeper: {d_abs} vs {d_clear}");
    }

    #[test]
    fn a_fully_exposed_column_clears_the_film() {
        // the missing test: monotonicity is not enough. A real dose must actually open the
        // resist, or every downstream stage fails with a message that blames them.
        let mack = DevelopModel::default();
        let dill = DillModel::default();
        let d = develop_depth(2.0, atto_units::Nm::new(150.0), &mack, &dill, atto_units::Minutes::new(1.0), 16).0;
        assert!(d >= 149.0, "a dose of 2.0 must clear a 150 nm film in a minute: got {d} nm");
    }

    #[test]
    fn an_unexposed_column_keeps_its_resist() {
        let mack = DevelopModel::default();
        let dill = DillModel::default();
        let d = develop_depth(0.0, atto_units::Nm::new(150.0), &mack, &dill, atto_units::Minutes::new(1.0), 16).0;
        assert!(d < 1.0, "no dose must leave the film alone: got {d} nm");
    }

    #[test]
    fn more_time_develops_deeper() {
        let mack = DevelopModel::default();
        let dill = DillModel::default();
        let short = develop_depth(2.0, atto_units::Nm::new(150.0), &mack, &dill, atto_units::Minutes::new(0.2), 32).0;
        let long = develop_depth(2.0, atto_units::Nm::new(150.0), &mack, &dill, atto_units::Minutes::new(0.8), 32).0;
        assert!(long > short, "longer develop goes deeper: {short} {long}");
    }

    #[test]
    fn the_depth_never_exceeds_the_film() {
        let mack = DevelopModel::default();
        let dill = DillModel::default();
        assert!(develop_depth(20.0, atto_units::Nm::new(150.0), &mack, &dill, atto_units::Minutes::new(100.0), 64).0 <= 150.0 + 1e-9);
    }

    #[test]
    fn a_real_bake_moves_the_compound_a_visible_distance() {
        // the first constant gave 9e-4 nm: technically a bake, physically nothing
        let sigma = bake_sigma_from_temperature(atto_units::Kelvin::new(110.0 + 273.15), atto_units::Minutes::new(1.0)).0;
        // the constant is documented as giving about 10 nm; the loose bound that stood here
        // (more than 1 nm) let a value 100 times too small through
        assert!(sigma > 5.0 && sigma < 20.0, "a 110 C bake for a minute should diffuse about 10 nm: {sigma} nm");
    }

    #[test]
    fn a_hotter_bake_spreads_further() {
        let cool = bake_sigma_from_temperature(atto_units::Kelvin::new(90.0 + 273.15), atto_units::Minutes::new(1.0)).0;
        let hot = bake_sigma_from_temperature(atto_units::Kelvin::new(130.0 + 273.15), atto_units::Minutes::new(1.0)).0;
        assert!(hot > cool, "hotter bake diffuses further: cool={cool} hot={hot}");
        assert!(cool > 0.0, "and a real bake is not zero");
    }
}
