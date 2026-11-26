// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Thick-mask (EMF) effects: the focus shift a thin mask cannot produce (L-08).
//!
//! A mask is not a plane. At 193 nm the absorber is tens of nanometres thick, and light entering
//! the space between two absorber lines at an angle sees a tunnel, not an aperture. The result is
//! a best focus that is **not** where the thin-mask model puts it, and the shift grows with the
//! mask thickness and with how obliquely the light arrives.
//!
//! This is why a scanner has to refocus when a mask is swapped for another of the same feature
//! size but a different absorber stack, and why a 6% attenuated mask focuses differently from a
//! binary one.

/// A mask, as a physical stack rather than a plane.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mask3d {
    /// Absorber thickness, nm.
    pub thickness_nm: f64,
    /// Transmission of the clear areas, 0..1. 1.0 is a binary mask, 0.06 is a 6% attenuated one.
    pub transmission: f64,
}

impl Mask3d {
    /// A binary mask with a typical absorber.
    pub fn binary() -> Self {
        Mask3d { thickness_nm: 60.0, transmission: 1.0 }
    }
    /// A 6% attenuated mask.
    pub fn attenuated_6() -> Self {
        Mask3d { thickness_nm: 60.0, transmission: 0.06 }
    }
}

/// Why a thick-mask shift could not be computed.
#[derive(Debug, Clone, PartialEq)]
pub struct Mask3dError {
    pub what: String,
}

impl std::fmt::Display for Mask3dError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "mask 3-D: {}", self.what)
    }
}

/// Coefficient of the focus shift, nm of focus per nm of absorber at NA 1.
///
/// **Assumed.** L-08 asks that a 6% attenuated mask show a focus shift against a thin-mask
/// baseline; the magnitude here is chosen to put that shift in the tens of nanometres, which is
/// the observed order. Registered as `mask3d.k`.
pub const MASK3D_K: f64 = 0.35;

/// Best-focus shift from a thick mask, nm.
///
/// `shift = k x thickness x (1 - transmission) / NA`. The dependence on `1 - transmission` is what
/// separates an attenuated mask from a binary one: a binary mask's clear areas carry full
/// amplitude and its shift is zero in this model, an attenuated mask's do not.
pub fn best_focus_shift_nm(mask: &Mask3d, na: f64) -> Result<f64, Mask3dError> {
    if !mask.thickness_nm.is_finite() || mask.thickness_nm < 0.0 {
        return Err(Mask3dError { what: "the absorber thickness must be finite and >= 0 nm".into() });
    }
    if !mask.transmission.is_finite() || !(0.0..=1.0).contains(&mask.transmission) {
        return Err(Mask3dError { what: "transmission must lie in 0..=1".into() });
    }
    if !na.is_finite() || na <= 0.0 {
        return Err(Mask3dError { what: "the numerical aperture must be finite and > 0".into() });
    }
    Ok(MASK3D_K * mask.thickness_nm * (1.0 - mask.transmission) / na)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_binary_mask_has_no_shift_in_this_model() {
        // the thin-mask baseline: full transmission, so the 1-t term vanishes
        let s = best_focus_shift_nm(&Mask3d::binary(), 0.93).expect("a shift");
        assert_eq!(s, 0.0, "a binary mask is the baseline: {s}");
    }

    #[test]
    fn a_six_percent_attenuated_mask_shifts_the_focus() {
        // L-08's acceptance: a shift against the thin-mask baseline
        let s = best_focus_shift_nm(&Mask3d::attenuated_6(), 0.93).expect("a shift");
        assert!(s > 10.0, "the shift is tens of nanometres, not a rounding error: {s}");
        assert!(s.is_finite());
    }

    #[test]
    fn the_shift_grows_with_the_absorber_thickness() {
        let thin = best_focus_shift_nm(&Mask3d { thickness_nm: 30.0, transmission: 0.06 }, 0.93).expect("thin");
        let thick = best_focus_shift_nm(&Mask3d { thickness_nm: 90.0, transmission: 0.06 }, 0.93).expect("thick");
        assert!(thick > thin, "a thicker absorber shifts more: {thick} vs {thin}");
        assert!((thick / thin - 3.0).abs() < 1e-12, "and linearly: {}", thick / thin);
    }

    #[test]
    fn a_higher_na_shifts_less() {
        // the same mask at a larger aperture: the light is steeper, the tunnel is shorter
        let low = best_focus_shift_nm(&Mask3d::attenuated_6(), 0.5).expect("low");
        let high = best_focus_shift_nm(&Mask3d::attenuated_6(), 1.0).expect("high");
        assert!(high < low, "more NA, less shift: {high} vs {low}");
    }

    #[test]
    fn a_less_transmissive_mask_shifts_more() {
        let a = best_focus_shift_nm(&Mask3d { thickness_nm: 60.0, transmission: 0.2 }, 0.93).expect("a");
        let b = best_focus_shift_nm(&Mask3d { thickness_nm: 60.0, transmission: 0.06 }, 0.93).expect("b");
        assert!(b > a, "darker clear areas shift more: {b} vs {a}");
    }

    #[test]
    fn impossible_input_is_refused() {
        assert!(best_focus_shift_nm(&Mask3d { thickness_nm: -1.0, transmission: 0.06 }, 0.93).is_err(), "negative thickness");
        assert!(best_focus_shift_nm(&Mask3d { thickness_nm: 60.0, transmission: 1.5 }, 0.93).is_err(), "transmission above one");
        assert!(best_focus_shift_nm(&Mask3d::binary(), 0.0).is_err(), "no aperture");
        assert!(best_focus_shift_nm(&Mask3d { thickness_nm: f64::NAN, transmission: 0.06 }, 0.93).is_err(), "not a number");
    }
}
