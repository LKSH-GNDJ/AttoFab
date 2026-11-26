// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! What interface charge does to a C-V curve (V-13).
//!
//! Two different things happen, and they are easy to confuse:
//!
//! * **Fixed oxide charge** shifts the whole curve sideways. It is a voltage offset: the gate
//!   must supply `-Qf/Cox` more or less to reach the same surface potential, so the flat-band
//!   voltage moves and the curve does not change shape.
//! * **Interface states** stretch the curve along the voltage axis. Charging them takes gate
//!   voltage that does not move the surface potential, so the same change in surface potential
//!   needs more volts — the curve is stretched, not shifted.
//!
//! A model that reports only a flat-band shift will fit the offset and still be wrong about the
//! shape, which is what a real C-V comparison shows first.
use atto_units::{FPerCm2, PerCm2, Volts};

use crate::mos::QE;

/// Why an interface calculation could not run.
#[derive(Debug, Clone, PartialEq)]
pub enum InterfaceError {
    /// A charge, capacitance or potential was not usable.
    BadInput(&'static str),
}

impl std::fmt::Display for InterfaceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            InterfaceError::BadInput(w) => write!(f, "interface effect: {w}"),
        }
    }
}

/// The two ways interface charge moves a C-V.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceEffect {
    /// Flat-band shift from fixed charge, volts. Positive charge gives a negative shift.
    pub vfb_shift_v: f64,
    /// Stretch-out from interface states at this surface potential, volts.
    pub stretch_out_v: f64,
    /// The two together, volts.
    pub total_v: f64,
}

/// The C-V effect of a fixed charge and an interface state density.
///
/// `cox_f_per_cm2` is the oxide capacitance per unit area; both effects are that capacitance
/// absorbing charge, so both scale as `1/Cox`.
pub fn interface_effect(
    qf_cm2: PerCm2,
    dit_cm2ev: f64,
    cox_f_per_cm2: FPerCm2,
    psi_s_v: Volts,
) -> Result<InterfaceEffect, InterfaceError> {
    let cox_f_per_cm2 = cox_f_per_cm2.0;
    let psi_s_v = psi_s_v.0;
    let qf_cm2 = qf_cm2.0;
    if !qf_cm2.is_finite() {
        return Err(InterfaceError::BadInput("fixed charge must be finite"));
    }
    if !dit_cm2ev.is_finite() || dit_cm2ev < 0.0 {
        return Err(InterfaceError::BadInput("interface state density must be finite and >= 0"));
    }
    if !cox_f_per_cm2.is_finite() || cox_f_per_cm2 <= 0.0 {
        return Err(InterfaceError::BadInput("oxide capacitance must be finite and > 0"));
    }
    if !psi_s_v.is_finite() {
        return Err(InterfaceError::BadInput("surface potential must be finite"));
    }
    let vfb_shift = -QE * qf_cm2 / cox_f_per_cm2;
    let stretch = QE * dit_cm2ev * psi_s_v / cox_f_per_cm2;
    Ok(InterfaceEffect {
        vfb_shift_v: vfb_shift,
        stretch_out_v: stretch,
        total_v: vfb_shift + stretch,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A 96 nm oxide: 3.9 * 8.854e-14 / 96e-7 F/cm^2
    fn cox() -> f64 {
        3.9 * 8.854_187_817e-14 / 96.0e-7
    }

    #[test]
    fn a_clean_interface_moves_nothing() {
        let e = interface_effect(atto_units::PerCm2(0.0), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).expect("no charge");
        assert_eq!(e.vfb_shift_v, 0.0);
        assert_eq!(e.stretch_out_v, 0.0);
        assert_eq!(e.total_v, 0.0);
    }

    #[test]
    fn positive_fixed_charge_shifts_flat_band_negative() {
        // the textbook sign: positive oxide charge moves V_fb to the left
        let e = interface_effect(atto_units::PerCm2(1.0e12), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).expect("charge");
        assert!(e.vfb_shift_v < 0.0, "positive charge gives a negative shift: {}", e.vfb_shift_v);
        assert_eq!(e.stretch_out_v, 0.0, "fixed charge does not stretch the curve");
    }

    #[test]
    fn the_shift_is_linear_in_the_charge() {
        let a = interface_effect(atto_units::PerCm2(1.0e12), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).expect("a");
        let b = interface_effect(atto_units::PerCm2(2.0e12), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).expect("b");
        assert!((b.vfb_shift_v / a.vfb_shift_v - 2.0).abs() < 1e-12, "twice the charge, twice the shift");
    }

    #[test]
    fn interface_states_stretch_rather_than_shift() {
        // Dit alone must not move flat-band, and must grow with surface potential
        let small = interface_effect(atto_units::PerCm2(0.0), 5.0e11, atto_units::FPerCm2(cox()), atto_units::Volts(0.1)).expect("small");
        let large = interface_effect(atto_units::PerCm2(0.0), 5.0e11, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).expect("large");
        assert_eq!(small.vfb_shift_v, 0.0, "interface states do not shift flat-band");
        assert!(large.stretch_out_v > small.stretch_out_v, "they stretch more at a larger surface potential: {} vs {}", large.stretch_out_v, small.stretch_out_v);
        assert!((large.stretch_out_v / small.stretch_out_v - 3.0).abs() < 1e-12, "and linearly so");
    }

    #[test]
    fn a_thicker_oxide_shifts_more() {
        let thin = interface_effect(atto_units::PerCm2(1.0e12), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).expect("thin");
        let thick = interface_effect(atto_units::PerCm2(1.0e12), 0.0, atto_units::FPerCm2(cox() / 2.0), atto_units::Volts(0.3)).expect("thick");
        assert!(thick.vfb_shift_v.abs() > thin.vfb_shift_v.abs(), "half the capacitance, twice the shift");
    }

    #[test]
    fn unusable_input_is_refused() {
        assert!(interface_effect(atto_units::PerCm2(f64::NAN), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).is_err(), "not a number");
        assert!(interface_effect(atto_units::PerCm2(0.0), -1.0, atto_units::FPerCm2(cox()), atto_units::Volts(0.3)).is_err(), "negative Dit");
        assert!(interface_effect(atto_units::PerCm2(0.0), 0.0, atto_units::FPerCm2(0.0), atto_units::Volts(0.3)).is_err(), "no capacitance");
        assert!(interface_effect(atto_units::PerCm2(0.0), 0.0, atto_units::FPerCm2(cox()), atto_units::Volts(f64::NAN)).is_err(), "not a potential");
    }
}
