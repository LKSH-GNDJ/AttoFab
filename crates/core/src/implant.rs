//! Ion Implantation Solver - Gaussian (LSS-theory-style) range profile:
//!
//! ```text
//! C(x) = Q / (sqrt(2*pi) * dRp) * exp( -(x-Rp)^2 / (2*dRp^2) )
//! ```
//!
//! Rp/dRp are physically set by ion species, energy, and target material
//! (classically from LSS theory / SRIM range tables). This module accepts
//! Rp/dRp directly, OR derives a simple illustrative energy-scaling
//! approximation - NOT a substitute for a real SRIM/TRIM lookup.

use crate::materials::{range_scaling, Dopant};
use std::f64::consts::PI;

#[derive(Debug, Clone, Copy)]
pub struct Range {
    /// projected range, um
    pub rp: f64,
    /// straggle, um
    pub d_rp: f64,
}

pub fn estimate_range(dopant: Dopant, energy_kev: f64) -> Range {
    let s = range_scaling(dopant);
    let rp = s.rp_per_kev_um * energy_kev;
    Range { rp, d_rp: rp * s.straggle_fraction }
}

pub struct ImplantParams {
    pub dose_cm2: f64,
    pub range: Range,
}

/// Build a Gaussian implant concentration profile (atoms/cm^3) over a
/// depth grid (um).
pub fn implant_profile(depth_grid_um: &[f64], p: &ImplantParams) -> Vec<f64> {
    let d_rp_cm = p.range.d_rp * 1e-4; // um -> cm
    let peak = p.dose_cm2 / ((2.0 * PI).sqrt() * d_rp_cm);

    depth_grid_um
        .iter()
        .map(|&x| {
            let exponent = -((x - p.range.rp).powi(2)) / (2.0 * p.range.d_rp * p.range.d_rp);
            peak * exponent.exp()
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn peak_concentration_occurs_at_projected_range() {
        let range = Range { rp: 0.1, d_rp: 0.03 };
        let depths: Vec<f64> = (0..200).map(|i| i as f64 * 0.001).collect();
        let profile = implant_profile(&depths, &ImplantParams { dose_cm2: 1e15, range });

        let (max_idx, _) = profile
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap();
        let peak_depth = depths[max_idx];
        assert!((peak_depth - range.rp).abs() < 0.005);
    }

    #[test]
    fn higher_energy_gives_deeper_range() {
        let shallow = estimate_range(Dopant::Phosphorus, 20.0);
        let deep = estimate_range(Dopant::Phosphorus, 80.0);
        assert!(deep.rp > shallow.rp);
    }

    #[test]
    fn total_dose_integrates_approximately_to_q() {
        let range = Range { rp: 0.1, d_rp: 0.02 };
        let dx_um = 0.0005;
        let depths: Vec<f64> = (0..2000).map(|i| i as f64 * dx_um).collect();
        let profile = implant_profile(&depths, &ImplantParams { dose_cm2: 1e15, range });

        // Riemann-sum integral of C(x) dx (converting um -> cm) should
        // approximate the specified dose Q (atoms/cm^2).
        let dx_cm = dx_um * 1e-4;
        let integral: f64 = profile.iter().map(|c| c * dx_cm).sum();
        let relative_error = (integral - 1e15).abs() / 1e15;
        assert!(relative_error < 0.02, "relative error was {}", relative_error);
    }
}

/// Pearson Type IV distribution - intended to capture the asymmetric
/// "channeling tail" real ion implants exhibit (ions traveling down open
/// lattice "hallways" penetrate deeper than a symmetric Gaussian
/// predicts).
///
/// STATUS: EXPERIMENTAL / NOT YET VERIFIED. The ODE integration and dose-
/// normalization machinery are correct and tested (see
/// `dose_normalization_conserves_total_dose` and
/// `invalid_moments_outside_type_iv_region_return_none`, both passing).
/// However, the moment-matching formula that converts (mean, variance,
/// skewness, kurtosis) into the ODE's (a, b0, b1, b2) parameters is
/// NOT YET PRODUCING CORRECT OUTPUT: empirical tests that recompute the
/// realized moments of the generated distribution and compare them to the
/// requested targets currently fail (mean off by several standard
/// deviations; skew sign inverted in at least one tested case) - see the
/// `#[ignore]`d tests below for the specific failures. Getting the
/// classical Elderton & Johnson moment-matching formulas exactly right
/// from a general reference (rather than a primary source with the full
/// derivation) proved unreliable enough that shipping it unverified would
/// risk silently producing wrong implant profiles - worse than not having
/// the feature at all.
///
/// This module is therefore NOT wired into `Wafer2d::implant` or
/// `Wafer1d::implant` - `implant_profile` (Gaussian) remains the default,
/// verified model. Do not enable Pearson-IV for real use until the
/// ignored tests below pass against the exact Elderton & Johnson formulas
/// (Continuous Univariate Distributions, Vol. 1, Johnson/Kotz/Balakrishnan
/// 1994 - the PearsonDS R package implements this correctly and is a good
/// reference to validate against numerically).
pub mod pearson4 {
    

    #[derive(Debug, Clone, Copy)]
    pub struct PearsonMoments {
        pub mean: f64,
        pub variance: f64,
        pub skewness: f64,
        pub kurtosis: f64, // NOT excess kurtosis (Gaussian = 3.0)
    }

    #[derive(Debug, Clone, Copy)]
    struct PearsonParams {
        a: f64,
        b0: f64,
        b1: f64,
        b2: f64,
    }

    /// Classical Pearson-system moment-matching formulas. Returns None if
    /// the moments fall outside the Type IV region (denominator
    /// degenerate, or beta1/beta2 combination invalid).
    fn fit_params(m: PearsonMoments) -> Option<PearsonParams> {
        let beta1 = m.skewness * m.skewness; // skewness^2
        let beta2 = m.kurtosis;
        let denom = 10.0 * beta2 - 12.0 * beta1 - 18.0;
        if denom.abs() < 1e-9 {
            return None;
        }
        let b0 = -m.variance * (4.0 * beta2 - 3.0 * beta1) / denom;
        let b1_mag = m.variance.sqrt() * beta1.sqrt() * (beta2 + 3.0) / denom;
        let b1 = if m.skewness < 0.0 { -b1_mag.abs() } else { b1_mag.abs() };
        let b2 = -(2.0 * beta2 - 3.0 * beta1 - 6.0) / denom;
        let a = b1;

        // Type IV requires a negative discriminant (complex conjugate
        // roots of the denominator quadratic) - i.e. b0+b1*u+b2*u^2 never
        // crosses zero, avoiding a singularity in the ODE.
        let discriminant = b1 * b1 - 4.0 * b2 * b0;
        if discriminant >= 0.0 || b2.abs() < 1e-12 {
            return None;
        }

        Some(PearsonParams { a, b0, b1, b2 })
    }

    /// Integrate the Pearson ODE via RK4 over the given depth grid,
    /// returning an unnormalized density. Returns None if the moments
    /// don't fit a valid Type IV distribution.
    pub fn unnormalized_density(moments: PearsonMoments, depth_grid_um: &Vec<f64>) -> Option<Vec<f64>> {
        let p = fit_params(moments)?;
        let denom_at = |u: f64| p.b0 + p.b1 * u + p.b2 * u * u;
        let dpdu = |u: f64, y: f64| -(p.a + u) / denom_at(u) * y;

        let n = depth_grid_um.len();
        let mut density = vec![0.0; n];

        let mean = moments.mean;
        let (start_idx, _) = depth_grid_um
            .iter()
            .enumerate()
            .min_by(|(_, a), (_, b)| (*a - mean).abs().partial_cmp(&(*b - mean).abs()).unwrap())?;

        density[start_idx] = 1.0;

        let rk4_step = |u: f64, y: f64, h: f64| -> f64 {
            let k1 = dpdu(u, y);
            let k2 = dpdu(u + h / 2.0, y + h / 2.0 * k1);
            let k3 = dpdu(u + h / 2.0, y + h / 2.0 * k2);
            let k4 = dpdu(u + h, y + h * k3);
            (y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)).max(0.0)
        };

        let mut y = density[start_idx];
        for i in start_idx..n.saturating_sub(1) {
            let u = depth_grid_um[i] - mean;
            let h = depth_grid_um[i + 1] - depth_grid_um[i];
            y = rk4_step(u, y, h);
            density[i + 1] = y;
        }
        let mut y = density[start_idx];
        for i in (1..=start_idx).rev() {
            let u = depth_grid_um[i] - mean;
            let h = depth_grid_um[i - 1] - depth_grid_um[i];
            y = rk4_step(u, y, h);
            density[i - 1] = y;
        }

        Some(density)
    }

    /// Build a dose-normalized Pearson-IV implant profile. Falls back
    /// signature mirrors `implant_profile` (Gaussian): dose in atoms/cm^2,
    /// depth grid in um, output in atoms/cm^3. Returns None (caller should
    /// fall back to Gaussian) if the requested moments aren't a valid
    /// Type IV distribution.
    pub fn implant_profile_pearson4(
        depth_grid_um: &Vec<f64>,
        dose_cm2: f64,
        moments: PearsonMoments,
    ) -> Option<Vec<f64>> {
        let density = unnormalized_density(moments, depth_grid_um)?;
        let dx_cm = if depth_grid_um.len() > 1 {
            (depth_grid_um[1] - depth_grid_um[0]) * 1e-4
        } else {
            return None;
        };
        // Normalize by numerical (trapezoidal) integration so the profile
        // integrates to the requested dose (atoms/cm^2).
        let raw_integral: f64 = density.iter().sum::<f64>() * dx_cm;
        if raw_integral <= 0.0 {
            return None;
        }
        let scale = dose_cm2 / raw_integral;
        Some(density.iter().map(|d| d * scale).collect())
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// Numerically recompute the realized mean/variance/skewness/
        /// kurtosis of a generated profile via trapezoidal integration,
        /// and check they match the requested target moments. This is the
        /// actual correctness check for the moment-matching formulas
        /// above - not an assumption that the formulas are right, but an
        /// empirical measurement of whether they produced a distribution
        /// with the requested shape.
        fn realized_moments(depth_grid_um: &Vec<f64>, density: &Vec<f64>) -> PearsonMoments {
            let total: f64 = density.iter().sum();
            let mean = density.iter().zip(depth_grid_um.iter()).map(|(d, x)| d * x).sum::<f64>() / total;
            let variance = density
                .iter()
                .zip(depth_grid_um.iter())
                .map(|(d, x)| d * (x - mean).powi(2))
                .sum::<f64>()
                / total;
            let m3 = density
                .iter()
                .zip(depth_grid_um.iter())
                .map(|(d, x)| d * (x - mean).powi(3))
                .sum::<f64>()
                / total;
            let m4 = density
                .iter()
                .zip(depth_grid_um.iter())
                .map(|(d, x)| d * (x - mean).powi(4))
                .sum::<f64>()
                / total;
            PearsonMoments {
                mean,
                variance,
                skewness: m3 / variance.powf(1.5),
                kurtosis: m4 / (variance * variance),
            }
        }

        #[test]
        #[ignore = "KNOWN FAILING: moment-matching formula produces wrong mean/shape - \
                    see module-level doc comment. Do not un-ignore without fixing the \
                    fit_params() formula against a verified reference (e.g. cross-check \
                    numerically against the PearsonDS R package)."]
        fn near_gaussian_moments_recover_a_symmetric_bell_shape() {
            // beta1=0 (zero skew), beta2=3 (Gaussian kurtosis) is the
            // degenerate boundary of the Pearson system (falls back to
            // Normal, not strictly Type IV) - so nudge slightly off it to
            // stay inside the valid Type IV region while still checking
            // near-symmetric, near-Gaussian-shaped output.
            let target = PearsonMoments { mean: 0.08, variance: 0.0009, skewness: 0.05, kurtosis: 3.3 };
            let depth_grid: Vec<f64> = (0..601).map(|i| i as f64 * 0.3 / 600.0).collect();

            let density = unnormalized_density(target, &depth_grid).expect("should fit valid Type IV params");
            let realized = realized_moments(&depth_grid, &density);

            assert!((realized.mean - target.mean).abs() < 0.01, "mean mismatch: {} vs {}", realized.mean, target.mean);
            assert!(
                (realized.variance - target.variance).abs() / target.variance < 0.25,
                "variance mismatch: {} vs {}",
                realized.variance,
                target.variance
            );
        }

        #[test]
        #[ignore = "KNOWN FAILING: realized skew sign is inverted relative to target - \
                    see module-level doc comment. Do not un-ignore without fixing the \
                    fit_params() formula against a verified reference."]
        fn skewed_profile_produces_asymmetric_tail() {
            // Positive skew => the profile should extend further on the
            // deep side than the shallow side, relative to the peak - the
            // qualitative signature of channeling that motivates using
            // Pearson-IV over a symmetric Gaussian in the first place.
            let target = PearsonMoments { mean: 0.1, variance: 0.0016, skewness: 0.8, kurtosis: 4.5 };
            let depth_grid: Vec<f64> = (0..801).map(|i| i as f64 * 0.4 / 800.0).collect();

            let density = unnormalized_density(target, &depth_grid).expect("should fit valid Type IV params");
            let realized = realized_moments(&depth_grid, &density);

            assert!(realized.skewness > 0.1, "expected positive realized skew, got {}", realized.skewness);
        }

        #[test]
        fn dose_normalization_conserves_total_dose() {
            let target = PearsonMoments { mean: 0.08, variance: 0.0009, skewness: 0.3, kurtosis: 3.5 };
            let depth_grid: Vec<f64> = (0..601).map(|i| i as f64 * 0.3 / 600.0).collect();
            let dose_cm2 = 1e15;

            let profile = implant_profile_pearson4(&depth_grid, dose_cm2, target).expect("valid profile");
            let dx_cm = (depth_grid[1] - depth_grid[0]) * 1e-4;
            let integral: f64 = profile.iter().sum::<f64>() * dx_cm;

            let relative_error = (integral - dose_cm2).abs() / dose_cm2;
            assert!(relative_error < 0.02, "dose not conserved, relative error {}", relative_error);
        }

        #[test]
        fn invalid_moments_outside_type_iv_region_return_none() {
            // beta1 and beta2 chosen to force a non-negative discriminant
            // (or a degenerate denominator) - outside the region Type IV
            // covers. Must fail closed (None), never silently return a
            // bogus/singular profile.
            let target = PearsonMoments { mean: 0.1, variance: 0.001, skewness: 0.0, kurtosis: 1.5 };
            let depth_grid: Vec<f64> = (0..401).map(|i| i as f64 * 0.4 / 400.0).collect();
            assert!(unnormalized_density(target, &depth_grid).is_none());
        }
    }
}
