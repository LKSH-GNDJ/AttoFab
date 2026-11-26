// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Step coverage from aspect ratio and a sticking coefficient (PRD S-01).
//!
//! A precursor that sticks on first contact coats the top of a feature and
//! starves the bottom; one with a low sticking probability bounces, redistributes
//! and coats more conformally. The model is the standard two-parameter form: the
//! bottom coverage falls with aspect ratio and with how readily the species
//! sticks.

#[derive(Debug, Clone, Copy, PartialEq)]
/// A precursor: how readily it sticks and how far it travels.
pub struct Precursor {
    /// Sticking probability, 0 (fully re-emitting) to 1 (sticking on contact).
    pub sticking: f64,
    /// Aspect ratio at which coverage falls to a half with `sticking = 1`.
    pub ar_half: f64,
}

impl Default for Precursor {
    fn default() -> Self {
        Precursor { sticking: 0.5, ar_half: 4.0 }
    }
}

impl Precursor {
    fn validate(&self) -> Result<(), String> {
        if !(0.0..=1.0).contains(&self.sticking) || !self.sticking.is_finite() {
            return Err("sticking must lie between 0 and 1".into());
        }
        if !(self.ar_half > 0.0) || !self.ar_half.is_finite() {
            return Err("the half-coverage aspect ratio must be positive".into());
        }
        Ok(())
    }

    /// Bottom coverage as a fraction of the top, 0..1.
    pub fn bottom_coverage(&self, aspect_ratio: f64) -> Result<f64, String> {
        self.validate()?;
        if !(aspect_ratio >= 0.0) || !aspect_ratio.is_finite() {
            return Err("the aspect ratio must be non-negative".into());
        }
        // at zero aspect ratio the bottom is the top; coverage then falls with
        // aspect ratio, more steeply the more readily the species sticks
        let x = self.sticking * aspect_ratio / self.ar_half;
        Ok(1.0 / (1.0 + x))
    }

    /// Sidewall coverage. A low sticking coefficient coats the sidewall better
    /// than the bottom, because the flux that reaches it has bounced.
    pub fn sidewall_coverage(&self, aspect_ratio: f64) -> Result<f64, String> {
        let bottom = self.bottom_coverage(aspect_ratio)?;
        // geometric mean of the top and the bottom: between the two, always
        Ok(bottom.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_flat_surface_is_fully_covered() {
        let p = Precursor::default();
        assert!((p.bottom_coverage(0.0).unwrap() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn coverage_falls_with_aspect_ratio() {
        let p = Precursor::default();
        let mut prev = 1.0;
        for ar in [0.5, 1.0, 2.0, 4.0, 8.0] {
            let c = p.bottom_coverage(ar).unwrap();
            assert!(c < prev, "coverage rose at AR {ar}: {prev:.3} -> {c:.3}");
            prev = c;
        }
    }

    #[test]
    fn a_less_sticky_precursor_coats_more_conformally() {
        let conformal = Precursor { sticking: 0.1, ar_half: 4.0 };
        let sticky = Precursor { sticking: 1.0, ar_half: 4.0 };
        let a = conformal.bottom_coverage(4.0).unwrap();
        let b = sticky.bottom_coverage(4.0).unwrap();
        assert!(a > b, "a low sticking coefficient fills better: {a:.3} vs {b:.3}");
        // at the half-coverage aspect ratio a fully sticky species gives 0.5
        assert!((b - 0.5).abs() < 1e-12, "AR 4 with sticking 1 gives {b:.3}");
    }

    #[test]
    fn the_sidewall_sits_between_the_top_and_the_bottom() {
        let p = Precursor::default();
        for ar in [1.0, 4.0, 10.0] {
            let b = p.bottom_coverage(ar).unwrap();
            let s = p.sidewall_coverage(ar).unwrap();
            assert!(s >= b && s <= 1.0, "AR {ar}: side {s:.3}, bottom {b:.3}");
        }
    }

    #[test]
    fn bad_input_is_rejected() {
        let p = Precursor::default();
        assert!(p.bottom_coverage(-1.0).is_err());
        assert!(p.bottom_coverage(f64::NAN).is_err());
        let mut q = p;
        q.sticking = 2.0;
        assert!(q.bottom_coverage(1.0).is_err());
        q.sticking = 0.5;
        q.ar_half = 0.0;
        assert!(q.bottom_coverage(1.0).is_err());
    }
}
