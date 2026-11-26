// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Ion-energy and angle-dependent sputter yield (E-08).
//!
//! A constant etch rate is wrong in two directions at once. The yield rises with ion energy above
//! a threshold, roughly as `sqrt(E - E_th)`; and it is not highest at normal incidence — it peaks
//! at an oblique angle, typically 60-80 degrees, because a tilted ion deposits its momentum nearer
//! the surface.
//!
//! That peak is the reason a tilted etch is faster than a vertical one, and the reason a mask
//! corner erodes faster than its face. A model with one rate cannot produce either.
use atto_units::{Degrees};

/// Parameters for the yield model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct YieldParams {
    /// Threshold energy below which nothing sputters, eV.
    pub threshold_ev: f64,
    /// Reference energy, eV — the yield is `y_ref` here at normal incidence.
    pub reference_ev: f64,
    /// Yield at the reference energy and normal incidence, atoms per ion.
    pub y_ref: f64,
    /// Angle at which the yield peaks, degrees from normal.
    pub peak_angle_deg: f64,
}

impl YieldParams {
    /// Defaults for a silicon etch.
    ///
    /// **Assumed.** E-08 asks that the yield peak at the published oblique angle and rise with
    /// energy; the shape here does both and the magnitudes are typical rather than fitted.
    /// Registered as `yield.*`.
    pub fn silicon() -> Self {
        YieldParams { threshold_ev: 20.0, reference_ev: 200.0, y_ref: 1.0, peak_angle_deg: 70.0 }
    }
}

/// Why a yield could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub struct YieldError {
    pub what: String,
}

impl std::fmt::Display for YieldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "yield: {}", self.what)
    }
}

fn check(p: &YieldParams) -> Result<(), YieldError> {
    if !p.threshold_ev.is_finite() || p.threshold_ev < 0.0 {
        return Err(YieldError { what: "the threshold energy must be finite and >= 0".into() });
    }
    if !p.reference_ev.is_finite() || p.reference_ev <= p.threshold_ev {
        return Err(YieldError { what: "the reference energy must exceed the threshold".into() });
    }
    if !p.y_ref.is_finite() || p.y_ref <= 0.0 {
        return Err(YieldError { what: "the reference yield must be finite and > 0".into() });
    }
    if !p.peak_angle_deg.is_finite() || p.peak_angle_deg <= 0.0 || p.peak_angle_deg >= 90.0 {
        return Err(YieldError { what: "the peak angle must lie strictly between 0 and 90 degrees".into() });
    }
    Ok(())
}

/// The angle factor, normalised to 1.0 at normal incidence.
///
/// `1/cos(theta)` up to the peak — the path length through the surface layer grows as the ion
/// tilts — and a linear fall to zero at grazing beyond it, where the ion no longer penetrates.
pub fn angle_factor(angle_deg: Degrees, p: &YieldParams) -> Result<f64, YieldError> {
    let angle_deg = angle_deg.0;
    check(p)?;
    if !angle_deg.is_finite() || !(0.0..90.0).contains(&angle_deg) {
        return Err(YieldError { what: "the angle must lie in [0, 90) degrees".into() });
    }
    if angle_deg <= p.peak_angle_deg {
        Ok(1.0 / angle_deg.to_radians().cos())
    } else {
        // fall from the peak value to zero at grazing
        let peak = 1.0 / p.peak_angle_deg.to_radians().cos();
        let t = (angle_deg - p.peak_angle_deg) / (90.0 - p.peak_angle_deg);
        Ok((peak * (1.0 - t)).max(0.0))
    }
}

/// The energy factor, normalised to 1.0 at the reference energy.
///
/// `sqrt((E - E_th) / (E_ref - E_th))`, and zero below the threshold.
pub fn energy_factor(energy_ev: f64, p: &YieldParams) -> Result<f64, YieldError> {
    check(p)?;
    if !energy_ev.is_finite() || energy_ev < 0.0 {
        return Err(YieldError { what: "the energy must be finite and >= 0 eV".into() });
    }
    if energy_ev <= p.threshold_ev {
        return Ok(0.0);
    }
    Ok(((energy_ev - p.threshold_ev) / (p.reference_ev - p.threshold_ev)).sqrt())
}

/// Sputter yield, atoms per ion.
pub fn sputter_yield(energy_ev: f64, angle_deg: Degrees, p: &YieldParams) -> Result<f64, YieldError> {
    let angle_deg = angle_deg.0;
    Ok(p.y_ref * energy_factor(energy_ev, p)? * angle_factor(atto_units::Degrees(angle_deg), p)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn si() -> YieldParams {
        YieldParams::silicon()
    }

    #[test]
    fn the_yield_peaks_at_an_oblique_angle() {
        // E-08's acceptance: a peak away from normal incidence, which a constant rate cannot give
        let p = si();
        let mut best = (0.0, 0.0f64);
        for a in 0..90 {
            let y = angle_factor(atto_units::Degrees(a as f64), &p).expect("a factor");
            if y > best.1 {
                best = (a as f64, y);
            }
        }
        assert!((best.0 - p.peak_angle_deg).abs() < 1.5, "the peak is at the stated angle: {} vs {}", best.0, p.peak_angle_deg);
        assert!(best.1 > 2.0, "and it is well above the normal-incidence value: {}", best.1);
    }

    #[test]
    fn nothing_sputters_below_the_threshold_energy() {
        let p = si();
        assert_eq!(energy_factor(10.0, &p).expect("below"), 0.0, "below threshold");
        assert_eq!(energy_factor(p.threshold_ev, &p).expect("at"), 0.0, "at threshold");
        assert!(energy_factor(30.0, &p).expect("above") > 0.0, "above threshold");
    }

    #[test]
    fn the_yield_rises_with_energy() {
        let p = si();
        let a = sputter_yield(50.0, atto_units::Degrees(0.0), &p).expect("a");
        let b = sputter_yield(200.0, atto_units::Degrees(0.0), &p).expect("b");
        let c = sputter_yield(800.0, atto_units::Degrees(0.0), &p).expect("c");
        assert!(a < b && b < c, "more energy, more yield: {a} {b} {c}");
        assert!((b - p.y_ref).abs() < 1e-12, "the reference energy gives the reference yield: {b}");
    }

    #[test]
    fn a_tilted_etch_is_faster_than_a_vertical_one() {
        let p = si();
        let normal = sputter_yield(200.0, atto_units::Degrees(0.0), &p).expect("normal");
        let tilted = sputter_yield(200.0, atto_units::Degrees(70.0), &p).expect("tilted");
        assert!(tilted > normal, "the oblique yield is higher: {tilted} vs {normal}");
    }

    #[test]
    fn grazing_incidence_sputters_nothing() {
        let p = si();
        assert_eq!(angle_factor(atto_units::Degrees(89.999), &p).expect("grazing").round(), 0.0, "the yield falls to zero at 90 degrees");
    }

    #[test]
    fn impossible_input_is_refused() {
        let p = si();
        assert!(angle_factor(atto_units::Degrees(90.0), &p).is_err(), "90 degrees is not an incidence");
        assert!(angle_factor(atto_units::Degrees(-1.0), &p).is_err(), "negative angle");
        assert!(energy_factor(-1.0, &p).is_err(), "negative energy");
        let mut bad = p;
        bad.reference_ev = bad.threshold_ev;
        assert!(energy_factor(100.0, &bad).is_err(), "reference at the threshold");
    }
}
