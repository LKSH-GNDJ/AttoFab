// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Small-signal and transient response (V-14).
//!
//! A MOS capacitor measured at increasing frequency does not hold its low-frequency value: above
//! the corner set by the series resistance and the capacitance, the measured capacitance falls as
//! `1/f`. The corner is `1/(2 pi R C)` and it is where the reactance equals the resistance — not a
//! fitted parameter, and not a property of the dielectric.
//!
//! The same corner appears in the time domain as a time constant. A C-f sweep and a switching
//! transient measure the same number two ways, and a model that reports them independently is
//! reporting one thing twice.
use atto_units::{Seconds};

use std::f64::consts::PI;

/// A representative series resistance for a C-V measurement, ohms.
///
/// **Assumed.** This is a measurement condition, not a material property: probe contact, cable and
/// substrate add up to something in this range. Registered as `meas.series.r` so a reader can see
/// what the reported corner depends on.
pub const DEFAULT_SERIES_R_OHM: f64 = 100.0;

/// Why a small-signal calculation could not run.
#[derive(Debug, Clone, PartialEq)]
pub struct SmallSignalError {
    pub what: String,
}

impl std::fmt::Display for SmallSignalError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "small signal: {}", self.what)
    }
}

fn check(r_ohm: f64, c_f: f64) -> Result<(), SmallSignalError> {
    if !r_ohm.is_finite() || r_ohm <= 0.0 {
        return Err(SmallSignalError { what: "series resistance must be finite and > 0 ohm".into() });
    }
    if !c_f.is_finite() || c_f <= 0.0 {
        return Err(SmallSignalError { what: "capacitance must be finite and > 0 F".into() });
    }
    Ok(())
}

/// The RC corner frequency, Hz.
pub fn corner_hz(r_ohm: f64, c_f: f64) -> Result<f64, SmallSignalError> {
    check(r_ohm, c_f)?;
    Ok(1.0 / (2.0 * PI * r_ohm * c_f))
}

/// The RC time constant, seconds.
pub fn time_constant_s(r_ohm: f64, c_f: f64) -> Result<f64, SmallSignalError> {
    check(r_ohm, c_f)?;
    Ok(r_ohm * c_f)
}

/// Measured capacitance at a frequency, F.
///
/// The series-RC divider gives `C_meas = C / sqrt(1 + (f/f_c)^2)`, so it is `C` well below the
/// corner and falls as `1/f` well above it.
pub fn capacitance_at(f_hz: f64, r_ohm: f64, c_f: f64) -> Result<f64, SmallSignalError> {
    check(r_ohm, c_f)?;
    if !f_hz.is_finite() || f_hz < 0.0 {
        return Err(SmallSignalError { what: "frequency must be finite and >= 0 Hz".into() });
    }
    let fc = corner_hz(r_ohm, c_f)?;
    let x = f_hz / fc;
    Ok(c_f / (1.0 + x * x).sqrt())
}

/// A charging transient: the voltage across the capacitor at a time after a step, V.
pub fn step_response_v(t_s: Seconds, r_ohm: f64, c_f: f64, v_step: f64) -> Result<f64, SmallSignalError> {
    let t_s = t_s.0;
    check(r_ohm, c_f)?;
    if !t_s.is_finite() || t_s < 0.0 {
        return Err(SmallSignalError { what: "time must be finite and >= 0 s".into() });
    }
    if !v_step.is_finite() {
        return Err(SmallSignalError { what: "the step must be finite".into() });
    }
    let tau = time_constant_s(r_ohm, c_f)?;
    Ok(v_step * (1.0 - (-t_s / tau).exp()))
}

/// Recover the time constant from a sampled transient, seconds.
///
/// Fits `ln(1 - v/v_step)` against `t`, whose slope is `-1/tau`. This is the check that the C-f
/// corner and the transient are the same number: a model that gets one right and the other wrong
/// has two unrelated parameters wearing one name.
pub fn time_constant_from_transient(
    samples: &[(f64, f64)],
    v_step: f64,
) -> Result<f64, SmallSignalError> {
    if samples.len() < 2 {
        return Err(SmallSignalError { what: "at least two samples are needed".into() });
    }
    if !v_step.is_finite() || v_step == 0.0 {
        return Err(SmallSignalError { what: "the step must be finite and non-zero".into() });
    }
    let mut n = 0.0;
    let mut sx = 0.0;
    let mut sy = 0.0;
    let mut sxx = 0.0;
    let mut sxy = 0.0;
    for (t, v) in samples {
        if !t.is_finite() || !v.is_finite() || *t < 0.0 {
            return Err(SmallSignalError { what: "a sample is not a finite pair".into() });
        }
        let ratio = 1.0 - v / v_step;
        if !(ratio > 0.0) {
            return Err(SmallSignalError { what: "a sample has reached or passed the step".into() });
        }
        let y = ratio.ln();
        n += 1.0;
        sx += t;
        sy += y;
        sxx += t * t;
        sxy += t * y;
    }
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-30 {
        return Err(SmallSignalError { what: "the samples do not span any time".into() });
    }
    let slope = (n * sxy - sx * sy) / denom;
    if !(slope < 0.0) {
        return Err(SmallSignalError { what: "the transient does not rise; no time constant".into() });
    }
    Ok(-1.0 / slope)
}

#[cfg(test)]
mod tests {
    use super::*;

    const R: f64 = 1.0e3;
    const C: f64 = 1.0e-9;

    #[test]
    fn the_corner_is_one_over_two_pi_r_c() {
        let fc = corner_hz(R, C).expect("a corner");
        let analytic = 1.0 / (2.0 * PI * 1.0e3 * 1.0e-9);
        assert!((fc - analytic).abs() / analytic < 1e-12, "analytic corner: {fc} vs {analytic}");
        assert!((fc - 159_154.9).abs() < 1.0, "about 159 kHz for 1k and 1nF: {fc}");
    }

    #[test]
    fn the_capacitance_rolls_off_at_the_corner() {
        // V-14's acceptance: at the corner the measured capacitance is 1/sqrt(2) of the low value
        let fc = corner_hz(R, C).expect("a corner");
        let below = capacitance_at(fc / 100.0, R, C).expect("below");
        let at = capacitance_at(fc, R, C).expect("at");
        assert!((below - C).abs() / C < 1e-4, "well below the corner it is C: {below}");
        assert!((at - C / 2.0_f64.sqrt()).abs() / C < 1e-9, "at the corner it is C/sqrt2: {at}");
    }

    #[test]
    fn well_above_the_corner_it_falls_as_one_over_f() {
        let fc = corner_hz(R, C).expect("a corner");
        // well above the corner, where the +1 under the root is negligible. At 10x and 100x the
        // ratio is 9.95, not 10 — the first version of this test asserted 10 there and failed,
        // correctly: 1/sqrt(1+x^2) is only asymptotically 1/x.
        let a = capacitance_at(fc * 100.0, R, C).expect("a");
        let b = capacitance_at(fc * 1000.0, R, C).expect("b");
        assert!((a / b - 10.0).abs() < 0.02, "a decade of frequency, a decade of fall: {}", a / b);
    }

    #[test]
    fn a_transient_recovers_the_analytic_time_constant() {
        // the second half of V-14's acceptance
        let tau = time_constant_s(R, C).expect("tau");
        let mut samples = Vec::new();
        for i in 1..=10 {
            let t = tau * i as f64 * 0.2;
            samples.push((t, step_response_v(atto_units::Seconds(t), R, C, 1.0).expect("v")));
        }
        let recovered = time_constant_from_transient(&samples, 1.0).expect("recovered");
        assert!((recovered - tau).abs() / tau < 1e-9, "the transient gives back tau: {recovered} vs {tau}");
    }

    #[test]
    fn the_corner_and_the_time_constant_agree() {
        // the same number two ways: f_c x tau = 1/(2 pi)
        let fc = corner_hz(R, C).expect("fc");
        let tau = time_constant_s(R, C).expect("tau");
        assert!((fc * tau - 1.0 / (2.0 * PI)).abs() < 1e-12, "they are the same corner: {}", fc * tau);
    }

    #[test]
    fn unusable_input_is_refused() {
        assert!(corner_hz(0.0, C).is_err(), "no resistance");
        assert!(corner_hz(R, 0.0).is_err(), "no capacitance");
        assert!(capacitance_at(-1.0, R, C).is_err(), "negative frequency");
        assert!(time_constant_from_transient(&[(0.0, 0.0)], 1.0).is_err(), "one sample");
        assert!(time_constant_from_transient(&[(0.0, 0.0), (1.0, 2.0)], 1.0).is_err(), "past the step");
    }
}
