//! Thermal Oxidation Solver - Deal & Grove (1965) analytical model.
//!
//! ```text
//! x_o^2 + A * x_o = B * (t + tau)
//! ```
//!
//! Closed-form solution:
//!
//! ```text
//! x_o(t) = (A/2) * ( sqrt(1 + (t+tau)/(A^2/4B)) - 1 )
//! ```
//!
//! where tau offsets the parabola so it passes through any oxide already
//! present at t=0. Purely analytical - matches how Deal-Grove is taught:
//! a first-principles physical model, not an iterative/fitted curve.

use crate::materials::{
    arrhenius, oxidation_kinetics, Ambient, CrystalOrientation, NATIVE_OXIDE_UM,
    ORIENTATION_111_TO_100_LINEAR_RATIO,
};

#[derive(Debug, Clone, Copy)]
pub struct RateConstants {
    /// um
    pub a: f64,
    /// um^2/hr
    pub b: f64,
}

pub fn rate_constants(t_kelvin: f64, ambient: Ambient, orientation: CrystalOrientation) -> RateConstants {
    let k = oxidation_kinetics(ambient);
    let b = arrhenius(k.c_b, k.ea_b, t_kelvin);
    let mut b_over_a = arrhenius(k.c_ba, k.ea_ba, t_kelvin);
    // (111) is reaction-limited faster than (100); B (diffusion-limited)
    // is orientation-independent, so only B/A is scaled.
    if orientation == CrystalOrientation::Si111 {
        b_over_a *= ORIENTATION_111_TO_100_LINEAR_RATIO;
    }
    RateConstants { a: b / b_over_a, b }
}

#[derive(Debug, Clone, Copy)]
pub struct OxidationResult {
    pub final_oxide_um: f64,
    pub rate_constants: RateConstants,
    pub tau_hr: f64,
}

pub struct GrowOxideParams {
    pub temperature_c: f64,
    pub time_hours: f64,
    pub ambient: Ambient,
    pub initial_oxide_um: f64,
    pub orientation: CrystalOrientation,
}

impl Default for GrowOxideParams {
    fn default() -> Self {
        Self {
            temperature_c: 1000.0,
            time_hours: 1.0,
            ambient: Ambient::Dry,
            initial_oxide_um: NATIVE_OXIDE_UM,
            orientation: CrystalOrientation::Si100,
        }
    }
}

/// Grow oxide for a given duration, returning the new total thickness (um).
pub fn grow_oxide(p: &GrowOxideParams) -> OxidationResult {
    let t_kelvin = p.temperature_c + 273.15;
    let rc = rate_constants(t_kelvin, p.ambient, p.orientation);

    let tau_hr = (p.initial_oxide_um.powi(2) + rc.a * p.initial_oxide_um) / rc.b;

    let t = p.time_hours;
    let final_oxide_um =
        (rc.a / 2.0) * (((1.0 + (t + tau_hr) / (rc.a * rc.a / (4.0 * rc.b))).sqrt()) - 1.0);

    OxidationResult { final_oxide_um, rate_constants: rc, tau_hr }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Regression benchmark, pinned to the *current* (PDF-sourced,
    /// corrected) Deal-Grove constants - NOT the original JS Phase-1
    /// engine's output, which used an earlier, less carefully sourced set
    /// of pre-exponential constants and will legitimately disagree by a
    /// few nm. See docs/math_references.md for the constants' provenance
    /// and the note on why this value changed from an earlier draft.
    #[test]
    fn dry_oxidation_1000c_45min_matches_benchmark() {
        let result = grow_oxide(&GrowOxideParams {
            temperature_c: 1000.0,
            time_hours: 0.75,
            ambient: Ambient::Dry,
            initial_oxide_um: NATIVE_OXIDE_UM,
            orientation: crate::materials::CrystalOrientation::Si100,
        });
        let nm = result.final_oxide_um * 1000.0;
        assert_relative_eq!(nm, 31.4, epsilon = 0.5);
    }

    #[test]
    fn wet_oxidation_grows_faster_than_dry_at_same_conditions() {
        let dry = grow_oxide(&GrowOxideParams {
            temperature_c: 1000.0,
            time_hours: 1.0,
            ambient: Ambient::Dry,
            ..Default::default()
        });
        let wet = grow_oxide(&GrowOxideParams {
            temperature_c: 1000.0,
            time_hours: 1.0,
            ambient: Ambient::Wet,
            ..Default::default()
        });
        assert!(wet.final_oxide_um > dry.final_oxide_um);
    }

    #[test]
    fn zero_time_returns_initial_oxide() {
        let result = grow_oxide(&GrowOxideParams {
            time_hours: 0.0,
            initial_oxide_um: 0.01,
            ..Default::default()
        });
        assert_relative_eq!(result.final_oxide_um, 0.01, epsilon = 1e-9);
    }
}
