// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Rule-based optical proximity correction (PRD L-07).
//!
//! A dense array and an isolated line of the same drawn width do not print the
//! same. The dense one receives more of its neighbours' light, so it prints wider
//! — and a process window sized on one will fail on the other. Rule-based OPC
//! corrects the *mask* by a width that depends on the local pitch, so the printed
//! result is closer to constant.
//!
//! This is the simple form the requirement asks for: a bias table indexed by
//! pitch, not a full model-based inverse. It corrects the systematic part and
//! leaves the rest visible.

use atto_units::{Nm};

/// A bias rule: how much to shrink the mask, nm, as a function of pitch.
#[derive(Debug, Clone, PartialEq)]
pub struct BiasRule {
    /// Pitch below which a feature is treated as dense, nm.
    pub dense_pitch_nm: f64,
    /// Pitch above which a feature is treated as isolated, nm.
    pub isolated_pitch_nm: f64,
    /// Bias applied to a dense feature, nm. Positive shrinks the mask.
    pub dense_bias_nm: f64,
    /// Bias applied to an isolated feature, nm.
    pub isolated_bias_nm: f64,
}

impl Default for BiasRule {
    fn default() -> Self {
        // a dense feature prints wider, so its mask is pulled in more
        BiasRule {
            dense_pitch_nm: 300.0,
            isolated_pitch_nm: 800.0,
            dense_bias_nm: 12.0,
            isolated_bias_nm: 0.0,
        }
    }
}

impl BiasRule {
    fn valid(&self) -> bool {
        self.dense_pitch_nm.is_finite() && self.isolated_pitch_nm.is_finite()
            && self.dense_pitch_nm > 0.0
            && self.isolated_pitch_nm > self.dense_pitch_nm
            && self.dense_bias_nm.is_finite() && self.isolated_bias_nm.is_finite()
    }

    /// The mask bias for a feature on this pitch, nm. Positive shrinks the mask.
    ///
    /// Interpolated linearly between the two anchors rather than stepped, because
    /// a step in the bias is a step in the printed CD at the same pitch.
    pub fn bias_nm(&self, pitch_nm: f64) -> Option<f64> {
        if !self.valid() || !pitch_nm.is_finite() || pitch_nm <= 0.0 {
            return None;
        }
        if pitch_nm <= self.dense_pitch_nm {
            return Some(self.dense_bias_nm);
        }
        if pitch_nm >= self.isolated_pitch_nm {
            return Some(self.isolated_bias_nm);
        }
        let t = (pitch_nm - self.dense_pitch_nm) / (self.isolated_pitch_nm - self.dense_pitch_nm);
        Some(self.dense_bias_nm + t * (self.isolated_bias_nm - self.dense_bias_nm))
    }

    /// The drawn width after the bias, nm. Never below a minimum, because a
    /// negative width is a mask defect rather than a correction.
    pub fn corrected_width_nm(&self, drawn: Nm, pitch_nm: f64) -> Option<f64> {
        let drawn_nm = drawn.value();
        let b = self.bias_nm(pitch_nm)?;
        if !(drawn_nm > 0.0) || !drawn_nm.is_finite() {
            return None;
        }
        Some((drawn_nm - b).max(1.0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_dense_feature_is_pulled_in_more_than_an_isolated_one() {
        let r = BiasRule::default();
        let dense = r.bias_nm(250.0).unwrap();
        let isolated = r.bias_nm(1000.0).unwrap();
        assert!(dense > isolated, "dense {dense:.1} vs isolated {isolated:.1}");
    }

    #[test]
    fn the_bias_interpolates_rather_than_stepping() {
        let r = BiasRule::default();
        let a = r.bias_nm(300.0).unwrap();
        let m = r.bias_nm(550.0).unwrap();
        let b = r.bias_nm(800.0).unwrap();
        assert!(m < a && m > b, "the midpoint sits between: {a:.1} > {m:.1} > {b:.1}");
        // and it is monotone across the range
        let mut prev = a;
        for p in [350.0, 450.0, 550.0, 650.0, 750.0] {
            let v = r.bias_nm(p).unwrap();
            assert!(v <= prev + 1e-12, "bias rose at pitch {p}");
            prev = v;
        }
    }

    #[test]
    fn the_corrected_width_never_goes_negative() {
        let r = BiasRule { dense_bias_nm: 500.0, ..BiasRule::default() };
        let w = r.corrected_width_nm(Nm::new(100.0), 250.0).unwrap();
        assert!(w >= 1.0, "a correction cannot make a negative feature: {w:.1}");
    }

    #[test]
    fn the_same_drawn_width_comes_out_different_on_two_pitches() {
        // the whole point: before correction, pitch changes the mask width
        let r = BiasRule::default();
        let dense = r.corrected_width_nm(Nm::new(125.0), 250.0).unwrap();
        let isolated = r.corrected_width_nm(Nm::new(125.0), 1000.0).unwrap();
        assert!(dense < isolated, "dense mask {dense:.1} vs isolated {isolated:.1}");
    }

    #[test]
    fn bad_input_is_refused() {
        let r = BiasRule::default();
        assert!(r.bias_nm(0.0).is_none());
        assert!(r.bias_nm(-10.0).is_none());
        assert!(r.bias_nm(f64::NAN).is_none());
        assert!(r.corrected_width_nm(Nm::new(0.0), 250.0).is_none());
        let bad = BiasRule { isolated_pitch_nm: 100.0, ..BiasRule::default() };
        assert!(bad.bias_nm(250.0).is_none(), "inverted anchors are a configuration error");
    }
}
