// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! EUV stochastics: photon shot noise and line-width roughness (L-10).
//!
//! At 13.5 nm a photon carries 91.8 eV, so a 20 mJ/cm^2 dose is a countable number of photons
//! per square nanometre — around 14. Counting statistics then set a floor on how straight a
//! printed edge can be, and the floor moves as `1 / sqrt(N)`.
//!
//! That square-root law is the whole content of this module, and it is worth stating plainly:
//! **doubling the dose reduces the roughness by 1.41, not by 2.** A model that reports roughness
//! as a constant, or that scales it linearly with dose, is wrong in a way that matters — it
//! predicts a resist that improves twice as fast as it can.
use atto_units::{MjPerCm2, Nm};

/// Photon energy at 13.5 nm, joules. `h c / lambda`.
pub const EUV_PHOTON_J: f64 = 1.463_3e-17;

/// Why a stochastic estimate could not run.
#[derive(Debug, Clone, PartialEq)]
pub struct StochasticError {
    pub what: String,
}

impl std::fmt::Display for StochasticError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "stochastic: {}", self.what)
    }
}

/// Photons arriving per square nanometre at a given dose.
///
/// `N = dose [J/m^2] x area [m^2] / photon energy`. The area is the pixel the roughness is
/// measured over, which is the resist correlation length rather than the scanner grid.
pub fn photons_per_pixel(dose_mj_cm2: MjPerCm2, pixel_nm: Nm) -> Result<f64, StochasticError> {
    let dose_mj_cm2 = dose_mj_cm2.0;
    let pixel_nm = pixel_nm.0;
    if !dose_mj_cm2.is_finite() || dose_mj_cm2 <= 0.0 {
        return Err(StochasticError { what: "dose must be finite and > 0 mJ/cm^2".into() });
    }
    if !pixel_nm.is_finite() || pixel_nm <= 0.0 {
        return Err(StochasticError { what: "pixel must be finite and > 0 nm".into() });
    }
    // mJ/cm^2 -> J/m^2 is a factor of 10
    let dose_j_m2 = dose_mj_cm2 * 10.0;
    let area_m2 = (pixel_nm * 1.0e-9) * (pixel_nm * 1.0e-9);
    let n = dose_j_m2 * area_m2 / EUV_PHOTON_J;
    if !n.is_finite() || n <= 0.0 {
        return Err(StochasticError { what: "the photon count is not a number".into() });
    }
    Ok(n)
}

/// Line-width roughness from shot noise, nm.
///
/// The roughness is the edge displacement that a one-sigma fluctuation in the photon count
/// produces, scaled by the resist's response. `k` is that response in nm per unit fractional
/// dose fluctuation; it is the one assumed constant here.
pub fn lwr_nm(dose_mj_cm2: MjPerCm2, pixel_nm: Nm, k_nm: Nm) -> Result<f64, StochasticError> {
    let dose_mj_cm2 = dose_mj_cm2.0;
    let k_nm = k_nm.0;
    let pixel_nm = pixel_nm.0;
    if !k_nm.is_finite() || k_nm <= 0.0 {
        return Err(StochasticError { what: "the resist response must be finite and > 0 nm".into() });
    }
    let n = photons_per_pixel(atto_units::MjPerCm2(dose_mj_cm2), atto_units::Nm(pixel_nm))?;
    Ok(k_nm / n.sqrt())
}

/// The assumed resist response, nm.
///
/// **Assumed, not fitted.** L-10 asks that LWR fall as `1/sqrt(dose)` and that the reported LER
/// match a published trend; the square-root law is exact in this model and the magnitude is set
/// by this one number, which is not calibrated.
pub const LWR_RESPONSE_NM: f64 = 3.0;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roughness_falls_as_the_inverse_square_root_of_dose() {
        // L-10's acceptance in its simplest form: quadruple the dose, halve the roughness
        let a = lwr_nm(atto_units::MjPerCm2(20.0), atto_units::Nm(10.0), atto_units::Nm(LWR_RESPONSE_NM)).expect("a");
        let b = lwr_nm(atto_units::MjPerCm2(80.0), atto_units::Nm(10.0), atto_units::Nm(LWR_RESPONSE_NM)).expect("b");
        assert!((a / b - 2.0).abs() < 1e-9, "4x dose gives 1/2 the LWR: {} vs {}", a, b);
    }

    #[test]
    fn doubling_the_dose_improves_by_root_two_not_two() {
        // the mistake this module exists to prevent
        let a = lwr_nm(atto_units::MjPerCm2(20.0), atto_units::Nm(10.0), atto_units::Nm(LWR_RESPONSE_NM)).expect("a");
        let b = lwr_nm(atto_units::MjPerCm2(40.0), atto_units::Nm(10.0), atto_units::Nm(LWR_RESPONSE_NM)).expect("b");
        assert!((a / b - 2.0_f64.sqrt()).abs() < 1e-9, "1/sqrt(2), not 2: {}", a / b);
        assert!(a / b < 2.0, "and strictly less than a factor of two");
    }

    #[test]
    fn the_photon_count_matches_the_photon_density() {
        // 20 mJ/cm^2 at 91.8 eV is 1.37e15 photons/cm^2, which is 13.7 per nm^2 — so a 10x10 nm
        // pixel holds about 1370. The first version of this test said "tens" and failed: the
        // module doc's "around 14 per square nanometre" is right and my pixel arithmetic was not.
        let n = photons_per_pixel(atto_units::MjPerCm2(20.0), atto_units::Nm(10.0)).expect("a count");
        assert!((n - 1370.0).abs() / 1370.0 < 0.02, "about 1370 photons in a 10 nm pixel: {n}");
        let per_nm2 = photons_per_pixel(atto_units::MjPerCm2(20.0), atto_units::Nm(1.0)).expect("a count");
        assert!((per_nm2 - 13.7).abs() / 13.7 < 0.02, "and 13.7 per square nanometre: {per_nm2}");
    }

    #[test]
    fn a_larger_pixel_has_more_photons_and_less_noise() {
        let small = lwr_nm(atto_units::MjPerCm2(20.0), atto_units::Nm(5.0), atto_units::Nm(LWR_RESPONSE_NM)).expect("small");
        let large = lwr_nm(atto_units::MjPerCm2(20.0), atto_units::Nm(20.0), atto_units::Nm(LWR_RESPONSE_NM)).expect("large");
        assert!(large < small, "averaging over more area is smoother: {large} vs {small}");
    }

    #[test]
    fn an_impossible_dose_is_refused() {
        assert!(photons_per_pixel(atto_units::MjPerCm2(0.0), atto_units::Nm(10.0)).is_err(), "no dose");
        assert!(photons_per_pixel(atto_units::MjPerCm2(-1.0), atto_units::Nm(10.0)).is_err(), "negative dose");
        assert!(photons_per_pixel(atto_units::MjPerCm2(f64::NAN), atto_units::Nm(10.0)).is_err(), "not a number");
        assert!(lwr_nm(atto_units::MjPerCm2(20.0), atto_units::Nm(10.0), atto_units::Nm(0.0)).is_err(), "no response");
    }
}
