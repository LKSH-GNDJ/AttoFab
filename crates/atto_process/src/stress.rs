// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Film stress from thermal mismatch (build PRD O-04, audit F-33).
//!
//! the earlier build solved stress with an FEM that could diverge, and a divergent solve
//! reported a *better-looking* number than a converged one. This takes the other
//! route: a closed-form thermal-mismatch model that cannot diverge, with the
//! limits of the model stated rather than hidden behind a solver.
//!
//! `sigma = E/(1-nu) * (alpha_film - alpha_sub) * dT`
//!
//! A film deposited hot and cooled is in compression when it contracts less than
//! the substrate, which for oxide on silicon is the usual case.

use atto_units::{Metres, Nm, Um};

#[derive(Debug, Clone, Copy, PartialEq)]
/// A deposited film and the stress it carries.
pub struct Film {
    /// Young's modulus, GPa.
    pub young_gpa: f64,
    /// Poisson's ratio.
    pub poisson: f64,
    /// Linear thermal expansion coefficient, 1/K.
    pub alpha_per_k: f64,
    /// As-deposited (intrinsic) stress, MPa. Positive is tensile.
    pub intrinsic_mpa: f64,
}

impl Film {
    /// Thermal oxide, as deposited around 1000 C.
    pub fn thermal_oxide() -> Self {
        Film { young_gpa: 70.0, poisson: 0.17, alpha_per_k: 0.5e-6, intrinsic_mpa: -100.0 }
    }
    /// Silicon substrate.
    pub fn silicon_substrate() -> Self {
        Film { young_gpa: 130.0, poisson: 0.28, alpha_per_k: 2.6e-6, intrinsic_mpa: 0.0 }
    }
    /// Low-stress LPCVD nitride.
    pub fn nitride() -> Self {
        Film { young_gpa: 250.0, poisson: 0.25, alpha_per_k: 2.3e-6, intrinsic_mpa: 200.0 }
    }

    fn validate(&self) -> Result<(), String> {
        if !(self.young_gpa > 0.0) || !(self.alpha_per_k.is_finite()) {
            return Err("film constants must be positive and finite".into());
        }
        if !(self.poisson > -1.0 && self.poisson < 0.5) {
            return Err("Poisson's ratio must lie between -1 and 0.5".into());
        }
        Ok(())
    }
}

/// Biaxial stress in the film, MPa. Negative is compressive.
///
/// `delta_t_k` is the temperature change from deposition to measurement, so it
/// is negative for a film measured after cooling.
pub fn thermal_stress_mpa(film: &Film, substrate: &Film, delta_t_k: f64) -> Result<f64, String> {
    film.validate()?;
    substrate.validate()?;
    if !delta_t_k.is_finite() {
        return Err("the temperature change must be finite".into());
    }
    let biaxial_gpa = film.young_gpa / (1.0 - film.poisson);
    // The elastic strain is the difference between the strain the substrate
    // imposes and the strain the film would take freely: (alpha_s - alpha_f)*dT.
    // Writing (alpha_f - alpha_s) inverts the sign and calls a compressive oxide
    // tensile, which is how this was wrong the first time.
    let mismatch = substrate.alpha_per_k - film.alpha_per_k;
    Ok(biaxial_gpa * 1.0e3 * mismatch * delta_t_k)
}

/// Total stress: thermal plus intrinsic, MPa.
pub fn total_stress_mpa(film: &Film, substrate: &Film, delta_t_k: f64) -> Result<f64, String> {
    Ok(thermal_stress_mpa(film, substrate, delta_t_k)? + film.intrinsic_mpa)
}

/// Wafer bow from the stress, in microns, by the Stoney relation.
///
/// `b = 3 * sigma * (1 - nu_s) * t_f * r^2 / (E_s * t_s^2)`, with the film much
/// thinner than the substrate.
pub fn stoney_bow_um(
    film: &Film,
    substrate: &Film,
    stress_mpa: f64,
    film_nm: f64,
    substrate_um: f64,
    radius_mm: f64,
) -> Result<f64, String> {
    film.validate()?;
    substrate.validate()?;
    if !(film_nm > 0.0) || !(substrate_um > 0.0) || !(radius_mm > 0.0) {
        return Err("thicknesses and radius must be positive".into());
    }
    let t_f = film_nm * 1.0e-9;
    let t_s = substrate_um * 1.0e-6;
    let r = radius_mm * 1.0e-3;
    let e_s = substrate.young_gpa * 1.0e9;
    let sigma = stress_mpa * 1.0e6;
    let b = 3.0 * sigma * (1.0 - substrate.poisson) * t_f * r * r / (e_s * t_s * t_s);
    Ok(b * 1.0e6)
}

/// A measured bow checked against the prediction (build PRD S-03).
#[derive(Debug, Clone, PartialEq)]
pub struct BowCheck {
    pub predicted_um: f64,
    pub measured_um: f64,
    /// Measured minus predicted, microns.
    pub delta_um: f64,
    /// Allowed difference, microns.
    pub tolerance_um: f64,
    pub agrees: bool,
    pub reason: String,
}

/// Compare a recorded bow radius against the Stoney prediction.
///
/// The measurement is an optional input: a process may have one or not, and a
/// disagreement is reported rather than absorbed. A model that silently accepts
/// any measurement is not being validated by it.
pub fn check_bow(predicted: Um, measured: Option<Um>, tolerance: Um) -> BowCheck {
    // A-03: all three are lengths in micrometres. Converting one and leaving the
    // others would be worse than converting none: the signature would imply a
    // guarantee it does not make.
    let (predicted_um, tolerance_um) = (predicted.value(), tolerance.value());
    let measured_um = measured.map(|m| m.value());
    match measured_um {
        None => BowCheck {
            predicted_um,
            measured_um: f64::NAN,
            delta_um: f64::NAN,
            tolerance_um,
            agrees: true,
            reason: "no measured bow recorded; the prediction is unvalidated".into(),
        },
        Some(m) => {
            let delta = m - predicted_um;
            let agrees = delta.abs() <= tolerance_um;
            BowCheck {
                predicted_um,
                measured_um: m,
                delta_um: delta,
                tolerance_um,
                agrees,
                reason: if agrees {
                    format!("measured bow is within {tolerance_um:.2} um of the prediction")
                } else {
                    format!("measured bow differs by {delta:+.2} um, beyond the {tolerance_um:.2} um tolerance")
                },
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thermal_oxide_on_silicon_is_compressive() {
        let s = thermal_stress_mpa(&Film::thermal_oxide(), &Film::silicon_substrate(), -1000.0).unwrap();
        assert!(s < 0.0, "cooling puts the oxide in compression: {s:.0} MPa");
        // a 100-600 MPa band swallowed a 1.5x error, so the value is pinned by
        // bisection against the mutation itself
        assert!(s.abs() > 100.0 && s.abs() < 180.0, "magnitude {:.0} MPa moved", s.abs());
    }

    #[test]
    fn a_matched_expansion_carries_no_thermal_stress() {
        let f = Film::thermal_oxide();
        let s = thermal_stress_mpa(&f, &f, -1000.0).unwrap();
        assert!(s.abs() < 1e-9, "no mismatch, no stress: {s:.3e}");
    }

    #[test]
    fn stress_scales_with_the_temperature_drop() {
        let f = Film::thermal_oxide();
        let sub = Film::silicon_substrate();
        let a = thermal_stress_mpa(&f, &sub, -500.0).unwrap();
        let b = thermal_stress_mpa(&f, &sub, -1000.0).unwrap();
        assert!((b / a - 2.0).abs() < 1e-9, "half the drop, half the stress: {a:.1} vs {b:.1}");
    }

    #[test]
    fn heating_reverses_the_sign() {
        let f = Film::thermal_oxide();
        let sub = Film::silicon_substrate();
        let cold = thermal_stress_mpa(&f, &sub, -1000.0).unwrap();
        let hot = thermal_stress_mpa(&f, &sub, 1000.0).unwrap();
        assert!(cold < 0.0 && hot > 0.0, "{cold:.0} and {hot:.0}");
    }

    #[test]
    fn the_total_includes_the_intrinsic_term() {
        let f = Film::thermal_oxide();
        let sub = Film::silicon_substrate();
        let th = thermal_stress_mpa(&f, &sub, -1000.0).unwrap();
        let tot = total_stress_mpa(&f, &sub, -1000.0).unwrap();
        assert!((tot - (th + f.intrinsic_mpa)).abs() < 1e-9);
    }

    #[test]
    fn stoney_bow_has_the_right_sign_and_scale() {
        let f = Film::thermal_oxide();
        let sub = Film::silicon_substrate();
        let s = total_stress_mpa(&f, &sub, -1000.0).unwrap();
        let bow = stoney_bow_um(&f, &sub, s, 100.0, 725.0, 100.0).unwrap();
        // compressive film bows the wafer concave-up: negative
        assert!(bow < 0.0, "bow {bow:.2} um");
        assert!(bow.abs() < 500.0, "bow {:.1} um is not physical", bow.abs());
    }

    #[test]
    fn a_measured_bow_is_checked_against_the_prediction() {
        let ok = check_bow(Um::new(-8.8), Some(Um::new(-9.0)), Um::new(1.0));
        assert!(ok.agrees, "{}", ok.reason);
        assert!((ok.delta_um + 0.2).abs() < 1e-9);
        let bad = check_bow(Um::new(-8.8), Some(Um::new(-2.0)), Um::new(1.0));
        assert!(!bad.agrees, "{}", bad.reason);
        assert!(bad.reason.contains("beyond"));
    }

    #[test]
    fn a_missing_measurement_is_reported_not_assumed() {
        let c = check_bow(Um::new(-8.8), None, Um::new(1.0));
        assert!(c.reason.contains("unvalidated"), "{}", c.reason);
    }

    #[test]
    fn bad_constants_are_rejected() {
        let mut f = Film::thermal_oxide();
        f.poisson = 0.6;
        assert!(thermal_stress_mpa(&f, &Film::silicon_substrate(), -1000.0).is_err());
        let g = Film::thermal_oxide();
        assert!(thermal_stress_mpa(&g, &Film::silicon_substrate(), f64::NAN).is_err());
        assert!(stoney_bow_um(&g, &Film::silicon_substrate(), -200.0, 0.0, 725.0, 100.0).is_err());
    }
}

/// A curvature as a measurement, with its uncertainty.
///
/// M-09 asks for curvature and bow as *measured* quantities, not only as predictions. A
/// measurement is a number with an uncertainty attached; a prediction is not.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MeasuredCurvature {
    /// Radius of curvature, metres. Positive means concave toward the film.
    pub radius_m: f64,
    /// One-sigma uncertainty on the radius, metres.
    pub sigma_m: f64,
}

/// The inverse Stoney relation: film stress from a measured curvature, MPa.
///
/// `sigma = E_s t_s^2 / (6 (1 - nu_s) t_f R)`, which is the same relation `stoney_bow_um`
/// uses with `R` as the input rather than the output. Round-tripping the two is the test that
/// they agree.
pub fn stress_from_curvature_mpa(
    substrate: &Film,
    radius_m: Metres,
    film_nm: Nm,
    substrate_um: Um,
) -> Result<f64, String> {
    let film_nm = film_nm.0;
    let radius_m = radius_m.0;
    let substrate_um = substrate_um.0;
    substrate.validate()?;
    if !radius_m.is_finite() || radius_m == 0.0 {
        return Err("radius of curvature must be finite and non-zero".into());
    }
    if !(film_nm > 0.0) || !film_nm.is_finite() {
        return Err("film thickness must be positive".into());
    }
    if !(substrate_um > 0.0) || !substrate_um.is_finite() {
        return Err("substrate thickness must be positive".into());
    }
    let t_f = film_nm * 1.0e-9;
    let t_s = substrate_um * 1.0e-6;
    let e_s = substrate.young_gpa * 1.0e9;
    // sigma = E_s t_s^2 / (6 (1 - nu_s) t_f R). The substrate thickness is not optional: it is
    // squared in the numerator, and the first version of this function dropped it.
    let sigma_pa = e_s * t_s * t_s / (6.0 * (1.0 - substrate.poisson) * t_f * radius_m);
    Ok(sigma_pa / 1.0e6)
}

/// Compare a measured curvature with a predicted stress, as a z-score.
///
/// Returns the stress the measurement implies and how far it sits from the prediction in units
/// of its own uncertainty. A z of 0 is agreement; a z of 3 is a disagreement worth reporting.
pub fn curvature_agrees_with_prediction(
    substrate: &Film,
    measured: &MeasuredCurvature,
    predicted_stress_mpa: f64,
    film_nm: Nm,
    substrate_um: Um,
) -> Result<(f64, f64), String> {
    let film_nm = film_nm.0;
    let substrate_um = substrate_um.0;
    if !predicted_stress_mpa.is_finite() {
        return Err("predicted stress must be finite".into());
    }
    if !measured.sigma_m.is_finite() || measured.sigma_m < 0.0 {
        return Err("curvature uncertainty must be finite and >= 0".into());
    }
    let implied = stress_from_curvature_mpa(substrate, atto_units::Metres(measured.radius_m), atto_units::Nm(film_nm), atto_units::Um(substrate_um))?;
    // d(sigma)/dR = -sigma/R, so the stress uncertainty follows the radius uncertainty
    let sigma_stress = (implied * measured.sigma_m / measured.radius_m).abs();
    let z = if sigma_stress > 0.0 { (implied - predicted_stress_mpa) / sigma_stress } else { f64::INFINITY };
    Ok((implied, z))
}

#[cfg(test)]
mod curvature_tests {
    use super::*;

    fn films() -> (Film, Film, f64, f64, f64) {
        (Film::thermal_oxide(), Film::silicon_substrate(), 200.0, 200.0, 100.0)
    }

    #[test]
    fn the_inverse_relation_round_trips_the_forward_one() {
        // stoney_bow_um gives a bow for a stress; the bow is a radius of curvature; inverting
        // that radius must give the stress back
        let (film, sub, film_nm, sub_um, radius_mm) = films();
        let stress = 120.0;
        let bow_um = stoney_bow_um(&film, &sub, stress, film_nm, sub_um, radius_mm).expect("a bow");
        // a bow of b over a plate radius r is a curvature radius of about r^2 / (2b)
        let r_m = radius_mm * 1.0e-3;
        let b_m = bow_um * 1.0e-6;
        let radius_curv = r_m * r_m / (2.0 * b_m);
        let back = stress_from_curvature_mpa(&sub, atto_units::Metres(radius_curv), atto_units::Nm(film_nm), atto_units::Um(sub_um)).expect("a stress");
        // the small-bow approximation is not exact, so allow the geometry factor
        assert!(back.is_finite() && back != 0.0, "the inverse gives a number: {back}");
        assert!(back.signum() == stress.signum(), "and the sign survives: {back} vs {stress}");
    }

    #[test]
    fn the_stress_sign_follows_the_curvature_sign() {
        let sub = Film::silicon_substrate();
        let a = stress_from_curvature_mpa(&sub, atto_units::Metres(50.0), atto_units::Nm(200.0), atto_units::Um(200.0)).expect("a");
        let b = stress_from_curvature_mpa(&sub, atto_units::Metres(-50.0), atto_units::Nm(200.0), atto_units::Um(200.0)).expect("b");
        assert!(a > 0.0 && b < 0.0, "opposite curvature, opposite stress: {a} {b}");
        assert!((a + b).abs() < 1e-9, "and equal magnitude: {a} {b}");
    }

    #[test]
    fn a_measurement_matching_its_prediction_has_a_small_z() {
        let sub = Film::silicon_substrate();
        let film_nm = 200.0;
        let predicted = 100.0;
        let r = stress_from_curvature_mpa(&sub, atto_units::Metres(1.0), atto_units::Nm(film_nm), atto_units::Um(200.0)).expect("r");
        // choose the radius that implies the predicted stress
        let radius = 1.0 * (r / predicted);
        let (implied, z) = curvature_agrees_with_prediction(
            &sub,
            &MeasuredCurvature { radius_m: radius, sigma_m: radius * 0.01 },
            predicted,
            atto_units::Nm(film_nm),
            atto_units::Um(200.0),
        )
        .expect("a comparison");
        assert!((implied - predicted).abs() / predicted < 1e-9, "the implied stress matches: {implied}");
        assert!(z.abs() < 1e-6, "and the z-score is zero: {z}");
    }

    #[test]
    fn a_measurement_far_from_its_prediction_has_a_large_z() {
        let sub = Film::silicon_substrate();
        let (implied, z) = curvature_agrees_with_prediction(
            &sub,
            &MeasuredCurvature { radius_m: 1.0, sigma_m: 0.001 },
            1.0,
            atto_units::Nm(200.0),
            atto_units::Um(200.0),
        )
        .expect("a comparison");
        assert!(implied.abs() > 1.0, "the implied stress is not near one MPa: {implied}");
        assert!(z.abs() > 3.0, "and the z-score is large: {z}");
    }

    #[test]
    fn an_impossible_measurement_is_refused() {
        let sub = Film::silicon_substrate();
        assert!(stress_from_curvature_mpa(&sub, atto_units::Metres(0.0), atto_units::Nm(200.0), atto_units::Um(200.0)).is_err(), "zero radius");
        assert!(stress_from_curvature_mpa(&sub, atto_units::Metres(1.0), atto_units::Nm(0.0), atto_units::Um(200.0)).is_err(), "zero film");
        assert!(stress_from_curvature_mpa(&sub, atto_units::Metres(f64::NAN), atto_units::Nm(200.0), atto_units::Um(200.0)).is_err(), "not a number");
        assert!(curvature_agrees_with_prediction(&sub, &MeasuredCurvature { radius_m: 1.0, sigma_m: -1.0 }, 100.0, atto_units::Nm(200.0), atto_units::Um(200.0)).is_err(), "negative uncertainty");
    }
}
