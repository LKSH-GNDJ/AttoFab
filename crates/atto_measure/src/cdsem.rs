// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! CD-SEM secondary-electron edge model (PRD M-06).
//!
//! A CD-SEM does not see an edge. It sees secondary electrons escaping from an
//! interaction volume, so the signal is the edge convolved with a resolution
//! kernel — and a critical dimension measured at a fixed fraction of that signal
//! is biased by the kernel width. Geometric extraction reports the drawn width and
//! hides the bias; this reports the width the instrument would give.
//!
//! The model is deliberately simple: a box edge blurred by a Gaussian, with the
//! edge brightness enhancement that makes a real SEM image look the way it does.

use atto_units::{Nm};

/// The instrument.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SemSettings {
    /// Edge resolution, nm. The interaction volume's half-width.
    pub sigma_nm: f64,
    /// Signal fraction at which the edge is located, 0..1.
    pub threshold: f64,
    /// Edge brightness enhancement: the peak above the plateau, as a fraction.
    pub edge_gain: f64,
}

impl Default for SemSettings {
    fn default() -> Self {
        SemSettings { sigma_nm: 5.0, threshold: 0.5, edge_gain: 0.35 }
    }
}

impl SemSettings {
    fn valid(&self) -> bool {
        self.sigma_nm.is_finite() && self.sigma_nm >= 0.0
            && self.threshold > 0.0 && self.threshold < 1.0
            && self.edge_gain.is_finite() && self.edge_gain >= 0.0
    }
}

fn erf(x: f64) -> f64 {
    // Abramowitz & Stegun 7.1.26, |error| < 1.5e-7
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t
            * (-x * x).exp();
    if x >= 0.0 { y } else { -y }
}

/// Secondary-electron signal across a line of width `line_nm` on a pitch, on a
/// grid of `nx` points spanning `field_nm`.
///
/// A bright feature on a dark background: the signal is high over the line.
pub fn signal(line: Nm, field: Nm, nx: usize, s: &SemSettings) -> Option<Vec<f64>> {
    let (line_nm, field_nm) = (line.value(), field.value());
    if !s.valid() || !(line_nm > 0.0) || !(field_nm > 0.0) || nx < 2 {
        return None;
    }
    let dx = field_nm / nx as f64;
    let x0 = field_nm * 0.5 - line_nm * 0.5;
    let x1 = x0 + line_nm;
    let sigma = s.sigma_nm.max(1e-9);
    let k = 1.0 / (2.0f64.sqrt() * sigma);
    let mut out = Vec::with_capacity(nx);
    for i in 0..nx {
        let x = (i as f64 + 0.5) * dx;
        // fraction of the interaction volume that sits over the feature
        let cover = 0.5 * (erf((x - x0) * k) - erf((x - x1) * k));
        // the edge enhancement peaks where the volume straddles an edge
        let near_left = (-((x - x0).powi(2)) / (2.0 * sigma * sigma)).exp();
        let near_right = (-((x - x1).powi(2)) / (2.0 * sigma * sigma)).exp();
        out.push(cover + s.edge_gain * cover.max(1e-12).min(1.0) * (near_left + near_right) * 0.5);
    }
    Some(out)
}

/// Width measured from the signal at the threshold fraction, nm.
///
/// The crossings are linearly interpolated between samples, so the result does not
/// quantise to the pixel grid.
pub fn measure_cd(signal: &[f64], field: Nm, threshold: f64) -> Option<f64> {
    let field_nm = field.value();
    if signal.len() < 3 || !(field_nm > 0.0) || !(threshold > 0.0 && threshold < 1.0) {
        return None;
    }
    let hi = signal.iter().cloned().fold(f64::MIN, f64::max);
    let lo = signal.iter().cloned().fold(f64::MAX, f64::min);
    if !(hi > lo) {
        return None;
    }
    let level = lo + threshold * (hi - lo);
    let dx = field_nm / signal.len() as f64;
    let mut crossings = Vec::new();
    for i in 0..signal.len() - 1 {
        let (a, b) = (signal[i], signal[i + 1]);
        if (a - level) * (b - level) < 0.0 {
            let f = (level - a) / (b - a);
            crossings.push((i as f64 + f) * dx);
        }
    }
    if crossings.len() < 2 {
        return None;
    }
    // the outermost crossings bound the feature
    Some(crossings[crossings.len() - 1] - crossings[0])
}

/// The CD a geometric extraction would report, for comparison.
pub fn geometric_cd(line: Nm, s: &SemSettings) -> Option<f64> {
    let line_nm = line.value();
    if !s.valid() || !(line_nm > 0.0) {
        return None;
    }
    Some(line_nm)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn cd_of(line: f64, s: &SemSettings) -> f64 {
        let field = 400.0;
        let sig = signal(Nm::new(line), Nm::new(field), 2000, s).expect("a signal");
        measure_cd(&sig, Nm::new(field), s.threshold).expect("a CD")
    }

    #[test]
    fn a_zero_width_kernel_recovers_the_geometric_cd() {
        let s = SemSettings { sigma_nm: 0.01, threshold: 0.5, edge_gain: 0.0 };
        let cd = cd_of(100.0, &s);
        assert!((cd - 100.0).abs() < 1.0, "sharp instrument gives the drawn width: {cd:.2}");
    }

    #[test]
    fn a_symmetric_blur_does_not_move_the_half_maximum_crossing() {
        // A property worth asserting because it is easy to assume the opposite:
        // blurring a box symmetrically leaves the 50 % point exactly on the edge.
        // So a CD read at half maximum is *insensitive* to resolution — the bias
        // has to come from somewhere else.
        let sharp = SemSettings { sigma_nm: 0.01, threshold: 0.5, edge_gain: 0.0 };
        let blunt = SemSettings { sigma_nm: 20.0, threshold: 0.5, edge_gain: 0.0 };
        let a = cd_of(100.0, &sharp);
        let b = cd_of(100.0, &blunt);
        assert!((b - a).abs() < 2.0, "blur alone barely moves it: {a:.2} vs {b:.2}");
    }

    #[test]
    fn the_bias_comes_from_the_threshold_and_the_edge_peak() {
        // with a threshold away from the half maximum, resolution does bias the
        // reading — the wider the kernel, the more the signal spreads
        let sharp = SemSettings { sigma_nm: 0.01, threshold: 0.25, edge_gain: 0.0 };
        let blunt = SemSettings { sigma_nm: 20.0, threshold: 0.25, edge_gain: 0.0 };
        let a = cd_of(100.0, &sharp);
        let b = cd_of(100.0, &blunt);
        assert!((b - a).abs() > 1.0, "off half-maximum the blur matters: {a:.2} vs {b:.2}");
    }

    #[test]
    fn the_measured_cd_grows_with_the_drawn_width() {
        let s = SemSettings::default();
        let mut prev = 0.0;
        for w in [50.0, 100.0, 200.0] {
            let cd = cd_of(w, &s);
            assert!(cd > prev, "width {w} measured {cd:.2}");
            prev = cd;
        }
    }

    #[test]
    fn the_edge_enhancement_moves_the_crossing() {
        let plain = SemSettings { sigma_nm: 10.0, threshold: 0.5, edge_gain: 0.0 };
        let bright = SemSettings { sigma_nm: 10.0, threshold: 0.5, edge_gain: 0.8 };
        let a = cd_of(100.0, &plain);
        let b = cd_of(100.0, &bright);
        assert!((a - b).abs() > 1e-6, "the edge peak changes the reading: {a:.2} vs {b:.2}");
    }

    #[test]
    fn bad_input_is_refused() {
        let s = SemSettings::default();
        assert!(signal(Nm::new(0.0), Nm::new(400.0), 100, &s).is_none());
        assert!(signal(Nm::new(100.0), Nm::new(400.0), 1, &s).is_none());
        let bad = SemSettings { sigma_nm: -1.0, ..s };
        assert!(signal(Nm::new(100.0), Nm::new(400.0), 100, &bad).is_none());
        assert!(measure_cd(&[1.0, 1.0, 1.0], Nm::new(400.0), 0.5).is_none(), "a flat signal has no edge");
        assert!(measure_cd(&[], Nm::new(400.0), 0.5).is_none());
        assert!(geometric_cd(Nm::new(0.0), &s).is_none());
    }
}
