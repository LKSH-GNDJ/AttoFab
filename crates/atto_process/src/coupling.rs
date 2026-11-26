// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Cross-stage coupling: one stage's output as the next stage's error (P-04).
//!
//! Every stage in this tree has been modelled as if the wafer arrived ideal. A real flow
//! does not: a polish leaves the field a few tens of nanometres out of plane, and the next
//! lithography step prints that error as a focus error, which shows up as a CD shift.
//!
//! This module is deliberately small. It carries one coupling — planarity error to printed
//! CD — because that is the one P-04 asks for and the first of its kind here. It does not
//! model resist thickness variation over topography, nor the swing curve, nor the fact that a
//! focus error and an exposure error are not separable in a single measurement.

use crate::litho::{cd_from_threshold, Grating, LithoError, Optics};
use atto_units::Nm;

/// Focus error caused by a planarity error, nm.
///
/// To first order the resist top sits exactly `erosion` below where the scanner expects it,
/// so the focus error equals the height error. The sign is negative: a surface that has been
/// polished back is *below* the focal plane.
pub fn defocus_from_erosion(erosion: Nm) -> Nm {
    Nm(-erosion.0)
}

/// The printed CD of a grating when the field carries a planarity error, nm.
///
/// The optics are copied with the focus shifted and the aerial image is recomputed; the CD
/// comes from the same threshold crossing the flat-field path uses, so the two are directly
/// comparable.
pub fn cd_with_erosion(
    optics: &Optics,
    grating: &Grating,
    threshold: f64,
    positive: bool,
    erosion: Nm,
    samples: usize,
) -> Result<f64, LithoError> {
    if !erosion.0.is_finite() {
        return Err(LithoError::BadInput("planarity error must be a finite length".into()));
    }
    if !threshold.is_finite() || threshold <= 0.0 {
        return Err(LithoError::BadInput("threshold must be finite and > 0".into()));
    }
    let mut shifted = *optics;
    shifted.defocus_nm = optics.defocus_nm + defocus_from_erosion(erosion).0;
    let img = shifted.image(grating, samples)?;
    Ok(cd_from_threshold(&img, Nm::new(grating.pitch_nm), threshold, positive))
}

/// How much a planarity error moves the printed CD, nm, signed.
///
/// Positive means the printed feature is wider than it would be on a flat field.
pub fn cd_shift_from_erosion(
    optics: &Optics,
    grating: &Grating,
    threshold: f64,
    positive: bool,
    erosion: Nm,
    samples: usize,
) -> Result<f64, LithoError> {
    let flat = cd_with_erosion(optics, grating, threshold, positive, Nm::new(0.0), samples)?;
    let tilted = cd_with_erosion(optics, grating, threshold, positive, erosion, samples)?;
    Ok(tilted - flat)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn optics() -> Optics {
        Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 }
    }
    fn grating() -> Grating {
        Grating { pitch_nm: 160.0, line_nm: 80.0 }
    }

    #[test]
    fn a_flat_field_prints_the_same_cd() {
        let shift = cd_shift_from_erosion(&optics(), &grating(), 0.3, true, Nm::new(0.0), 256).expect("a shift");
        assert!(shift.abs() < 1e-9, "no erosion, no shift: {shift}");
    }

    #[test]
    fn a_hundred_nanometre_erosion_degrades_the_printed_cd() {
        // P-04's acceptance: a 100 nm erosion field degrades the printed CD by the defocus it
        // causes. The number is whatever the optics give; what the requirement asks for is that
        // it is not zero and that it is the defocus doing it.
        let shift = cd_shift_from_erosion(&optics(), &grating(), 0.3, true, Nm::new(100.0), 256).expect("a shift");
        assert!(shift.abs() > 1.0, "100 nm of erosion must move the CD by more than a nanometre: {shift}");
    }

    #[test]
    fn the_shift_grows_with_the_erosion() {
        let o = optics();
        let g = grating();
        let mut last = 0.0f64;
        for e in [0.0, 25.0, 50.0, 100.0] {
            let s = cd_shift_from_erosion(&o, &g, 0.3, true, Nm::new(e), 256).expect("a shift");
            assert!(s.abs() >= last.abs() - 1e-9, "more erosion cannot print better: {e} nm gave {s}, after {last}");
            last = s;
        }
        assert!(last.abs() > 0.0, "the largest erosion moved the CD");
    }

    #[test]
    fn the_focus_error_has_the_sign_of_the_polished_surface() {
        // a polished-back surface sits below the focal plane
        assert_eq!(defocus_from_erosion(Nm::new(100.0)).0, -100.0);
        assert_eq!(defocus_from_erosion(Nm::new(0.0)).0, 0.0);
    }

    #[test]
    fn a_field_that_cannot_exist_is_refused() {
        let o = optics();
        let g = grating();
        assert!(cd_with_erosion(&o, &g, 0.3, true, Nm::new(f64::NAN), 256).is_err(), "not a length");
        assert!(cd_with_erosion(&o, &g, 0.0, true, Nm::new(0.0), 256).is_err(), "no threshold");
        assert!(cd_with_erosion(&o, &g, f64::NAN, true, Nm::new(0.0), 256).is_err(), "threshold not a number");
    }
}
