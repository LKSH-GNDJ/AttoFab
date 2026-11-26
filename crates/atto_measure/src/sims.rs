// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! SIMS comparison harness (M-08).
//!
//! A secondary-ion mass spectrometry profile is the only measurement that says where the
//! dopant actually is. Comparing one with a simulation is not a subtraction: the two axes
//! are both uncertain. The depth axis is offset by the sputter-rate calibration, and the
//! concentration axis spans four decades, so a linear residual is dominated by the peak and
//! says nothing about the tail.
//!
//! This module therefore does two things before it reports anything: it searches the depth
//! offset that best aligns the two profiles, and it measures the residual in log10
//! concentration, which is the scale the measurement is quoted on.
use atto_units::{Nm};

/// A depth profile: concentration against depth below the surface.
#[derive(Debug, Clone, PartialEq)]
pub struct Profile {
    pub depth_nm: Vec<f64>,
    pub conc_cm3: Vec<f64>,
}

/// Why a comparison could not be made.
#[derive(Debug, Clone, PartialEq)]
pub enum SimsError {
    /// A profile had fewer than two points, or its two axes differ in length.
    BadProfile(&'static str),
    /// Depth did not increase, so the profile cannot be interpolated.
    NotMonotone,
    /// A concentration was zero or negative, which has no logarithm.
    BadConcentration,
    /// The alignment window was not a positive finite length.
    BadWindow,
}

impl std::fmt::Display for SimsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SimsError::BadProfile(w) => write!(f, "unusable profile: {w}"),
            SimsError::NotMonotone => write!(f, "depth must increase from point to point"),
            SimsError::BadConcentration => write!(f, "concentration must be positive everywhere"),
            SimsError::BadWindow => write!(f, "the alignment window must be a positive finite length"),
        }
    }
}

/// How well a simulated profile matches a measured one.
#[derive(Debug, Clone, PartialEq)]
pub struct Residual {
    /// Root-mean-square difference in log10(concentration), decades.
    pub rms_log10: f64,
    /// Largest single difference, decades.
    pub max_log10: f64,
    /// Depth offset that was removed before measuring, nm.
    ///
    /// Positive means the simulation sits this much deeper than the measurement and was pulled
    /// back toward the surface to align.
    pub depth_shift_nm: f64,
    /// Points compared.
    pub points: usize,
    /// Uniform concentration offset that was removed, decades.
    ///
    /// Positive means the simulation is high by that many decades. SIMS absolute concentration
    /// carries a calibration uncertainty, so a factor common to the whole profile is a scale
    /// error, not a shape error; removing it separates "the profile is the wrong shape" from
    /// "the profile is right and the calibration is off".
    pub scale_log10: f64,
}

impl Profile {
    fn check(&self) -> Result<(), SimsError> {
        if self.depth_nm.len() != self.conc_cm3.len() {
            return Err(SimsError::BadProfile("depth and concentration differ in length"));
        }
        if self.depth_nm.len() < 2 {
            return Err(SimsError::BadProfile("fewer than two points"));
        }
        for w in self.depth_nm.windows(2) {
            if !(w[1] > w[0]) {
                return Err(SimsError::NotMonotone);
            }
        }
        if self.depth_nm.iter().chain(self.conc_cm3.iter()).any(|v| !v.is_finite()) {
            return Err(SimsError::BadProfile("a value is not finite"));
        }
        if self.conc_cm3.iter().any(|c| *c <= 0.0) {
            return Err(SimsError::BadConcentration);
        }
        Ok(())
    }

    /// Concentration at a depth, linearly interpolated; `None` outside the profile.
    fn at(&self, depth: f64) -> Option<f64> {
        let n = self.depth_nm.len();
        if depth < self.depth_nm[0] || depth > self.depth_nm[n - 1] {
            return None;
        }
        let mut lo = 0usize;
        while lo + 1 < n && self.depth_nm[lo + 1] < depth {
            lo += 1;
        }
        let (d0, d1) = (self.depth_nm[lo], self.depth_nm[lo + 1]);
        let (c0, c1) = (self.conc_cm3[lo], self.conc_cm3[lo + 1]);
        if d1 <= d0 {
            return Some(c0);
        }
        let t = (depth - d0) / (d1 - d0);
        Some(c0 + t * (c1 - c0))
    }
}

/// Compare a simulated profile against a measured one, removing the best depth offset.
///
/// The offset is searched over `±window_nm` in steps of a quarter of the measured spacing.
/// The residual is in log10 decades, so a factor of two anywhere in the profile costs the
/// same 0.30 whether it is at the peak or in the tail.
pub fn compare(measured: &Profile, simulated: &Profile, window_nm: Nm) -> Result<Residual, SimsError> {
    let window_nm = window_nm.0;
    measured.check()?;
    simulated.check()?;
    if !window_nm.is_finite() || window_nm < 0.0 {
        return Err(SimsError::BadWindow);
    }
    let spacing = {
        let n = measured.depth_nm.len();
        (measured.depth_nm[n - 1] - measured.depth_nm[0]) / (n - 1) as f64
    };
    let step = (spacing / 4.0).max(1e-6);
    let steps = (window_nm / step).round() as isize;

    let mut best: Option<Residual> = None;
    for k in -steps..=steps {
        let shift = k as f64 * step;
        // First pass: the mean log10 difference at this shift, which is the scale error.
        let mut diffs: Vec<f64> = Vec::new();
        for i in 0..measured.depth_nm.len() {
            let d = measured.depth_nm[i];
            let Some(c) = simulated.at(d - shift) else { continue };
            diffs.push(measured.conc_cm3[i].log10() - c.log10());
        }
        if diffs.len() < 2 {
            continue;
        }
        let mean = diffs.iter().sum::<f64>() / diffs.len() as f64;
        // Second pass: the residual after the scale is removed, so a factor of two cannot be
        // mistaken for a depth offset. Without this the search finds shifts that trade one
        // error for the other and reports a better alignment than exists.
        let mut sum = 0.0;
        let mut worst: f64 = 0.0;
        for d in &diffs {
            let r = d - mean;
            sum += r * r;
            if r.abs() > worst {
                worst = r.abs();
            }
        }
        let n = diffs.len();
        let rms = (sum / n as f64).sqrt();
        let candidate = Residual { rms_log10: rms, max_log10: worst, depth_shift_nm: -shift, points: n, scale_log10: -mean };
        let better = match &best {
            None => true,
            Some(b) => rms < b.rms_log10 - 1e-12,
        };
        if better {
            best = Some(candidate);
        }
    }
    best.ok_or(SimsError::BadProfile("no overlap between the two profiles"))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A plausible as-implanted profile: a Gaussian peak below the surface.
    fn synthetic(depth_nm: f64, peak_nm: f64, sigma_nm: f64, peak_cm3: f64) -> Profile {
        let mut depth = Vec::new();
        let mut conc = Vec::new();
        let mut d = 0.0;
        while d <= depth_nm {
            let c = peak_cm3 * (-((d - peak_nm).powi(2)) / (2.0 * sigma_nm * sigma_nm)).exp();
            depth.push(d);
            conc.push(c.max(1.0e14));
            d += 1.0;
        }
        Profile { depth_nm: depth, conc_cm3: conc }
    }

    #[test]
    fn a_profile_against_itself_is_exact() {
        let p = synthetic(200.0, 60.0, 25.0, 1.0e19);
        let r = compare(&p, &p, atto_units::Nm(10.0)).expect("a comparison");
        assert!(r.rms_log10 < 1e-9, "an identical profile has no residual: {}", r.rms_log10);
        assert_eq!(r.depth_shift_nm, 0.0, "and no shift");
        assert!(r.scale_log10.abs() < 1e-9, "and no scale error");
    }

    #[test]
    fn a_known_depth_offset_is_recovered() {
        // the simulation sits 8 nm deeper than the measurement
        let measured = synthetic(200.0, 60.0, 25.0, 1.0e19);
        let simulated = synthetic(200.0, 68.0, 25.0, 1.0e19);
        let r = compare(&measured, &simulated, atto_units::Nm(20.0)).expect("a comparison");
        assert!((r.depth_shift_nm - 8.0).abs() <= 0.5, "the offset is found: {}", r.depth_shift_nm);
        assert!(r.rms_log10 < 0.05, "and removing it leaves almost nothing: {}", r.rms_log10);
    }

    #[test]
    fn a_factor_of_two_costs_the_same_everywhere() {
        // doubling the whole profile is 0.301 decades, and that is what the residual reports
        let measured = synthetic(200.0, 60.0, 25.0, 1.0e19);
        let simulated = Profile {
            depth_nm: measured.depth_nm.clone(),
            conc_cm3: measured.conc_cm3.iter().map(|c| c * 2.0).collect(),
        };
        let r = compare(&measured, &simulated, atto_units::Nm(5.0)).expect("a comparison");
        assert!((r.scale_log10 - 0.30103).abs() < 1e-4, "log10(2) is 0.301 decades: {}", r.scale_log10);
        assert!(r.rms_log10 < 1e-6, "and once the scale is removed the shape is exact: {}", r.rms_log10);
        assert_eq!(r.depth_shift_nm, 0.0, "a scale error must not be reported as a shift");
    }

    #[test]
    fn a_tail_error_is_visible_even_when_the_peak_is_right() {
        // same peak, one decade low in the tail: a linear residual would call this a good fit
        let measured = synthetic(300.0, 60.0, 40.0, 1.0e19);
        let simulated = Profile {
            depth_nm: measured.depth_nm.clone(),
            conc_cm3: measured
                .conc_cm3
                .iter()
                .enumerate()
                .map(|(i, c)| if i > 200 { c * 0.1 } else { *c })
                .collect(),
        };
        let r = compare(&measured, &simulated, atto_units::Nm(5.0)).expect("a comparison");
        assert!(r.max_log10 > 0.5, "a decade of tail error is visible after the scale is removed: {}", r.max_log10);
    }

    #[test]
    fn unusable_input_is_refused() {
        let good = synthetic(100.0, 30.0, 20.0, 1.0e19);
        let mut short = good.clone();
        short.depth_nm.truncate(1);
        short.conc_cm3.truncate(1);
        assert!(compare(&short, &good, atto_units::Nm(5.0)).is_err(), "one point is not a profile");

        let mut backwards = good.clone();
        backwards.depth_nm[3] = backwards.depth_nm[1];
        assert!(compare(&backwards, &good, atto_units::Nm(5.0)).is_err(), "depth must increase");

        let mut negative = good.clone();
        negative.conc_cm3[5] = -1.0;
        assert!(compare(&negative, &good, atto_units::Nm(5.0)).is_err(), "no logarithm of a negative");

        assert!(compare(&good, &good, atto_units::Nm(-1.0)).is_err(), "a negative window");
    }
}
