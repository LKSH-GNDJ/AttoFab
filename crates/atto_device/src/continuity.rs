// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Carrier continuity along a channel: the coupled half of V-15.
//!
//! Poisson alone gives a potential and stops. A transistor needs the second equation — the
//! continuity equation for carriers — and the two are coupled through the charge. This module
//! solves the continuity half in the gradual-channel form: the electron quasi-Fermi level varies
//! along the channel, and the current is `I = q mu n d(phi_n)/dx` with
//! `n = n_i exp((psi - phi_n)/V_T)`.
//!
//! That is one dimension, not three, and the potential comes from the caller rather than from a
//! coupled Newton solve. It is the part of V-15 that can be stated and checked on its own: **the
//! current must rise with gate voltage, saturate with drain voltage, and be exponential below
//! threshold.**
use atto_units::{Volts};

/// A channel with a given surface potential.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Channel {
    /// Surface potential from the gate, V (positive = strong inversion for an n-channel).
    pub surface_potential_v: f64,
    /// Channel length, nm.
    pub length_nm: f64,
    /// Electron mobility, cm^2/Vs.
    pub mobility_cm2_vs: f64,
    /// Substrate doping magnitude, cm^-3.
    pub doping_cm3: f64,
    /// Temperature, K.
    pub temp_k: f64,
}

/// Why a continuity solve could not run.
#[derive(Debug, Clone, PartialEq)]
pub struct ContinuityError {
    pub what: String,
}

impl std::fmt::Display for ContinuityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "continuity: {}", self.what)
    }
}

const QE: f64 = 1.602_176_634e-19;
const KB: f64 = 1.380_649e-23;
const NI_CM3: f64 = 1.0e10;

impl Channel {
    fn check(&self) -> Result<(), ContinuityError> {
        if !self.surface_potential_v.is_finite() {
            return Err(ContinuityError { what: "the surface potential must be finite".into() });
        }
        if !self.length_nm.is_finite() || self.length_nm <= 0.0 {
            return Err(ContinuityError { what: "the channel length must be finite and > 0".into() });
        }
        if !self.mobility_cm2_vs.is_finite() || self.mobility_cm2_vs <= 0.0 {
            return Err(ContinuityError { what: "the mobility must be finite and > 0".into() });
        }
        if !self.doping_cm3.is_finite() || self.doping_cm3 <= 0.0 {
            return Err(ContinuityError { what: "the doping must be finite and > 0".into() });
        }
        if !self.temp_k.is_finite() || self.temp_k <= 0.0 {
            return Err(ContinuityError { what: "the temperature must be above absolute zero".into() });
        }
        Ok(())
    }

    /// Thermal voltage, V.
    pub fn thermal_v(&self) -> f64 {
        KB * self.temp_k / QE
    }

    /// Inversion charge sheet density at the source, cm^-2.
    ///
    /// `n_s = n_i^2/N_A x (depletion width) x exp(psi/V_T)` in the subthreshold form, capped by
    /// the gate's own capacitance in strong inversion — the two regimes the model has to join.
    fn sheet_charge_cm2(&self, phi_n_v: f64) -> f64 {
        let vt = self.thermal_v();
        let debye = ((11.7 * 8.854e-14 * vt) / (QE * self.doping_cm3)).sqrt();
        let n_sub = (NI_CM3 * NI_CM3 / self.doping_cm3) * debye;
        let arg = (self.surface_potential_v - phi_n_v) / vt;
        // exp can overflow at a large forward bias; the charge is enormous there and the current
        // is limited by the channel, so clamp the exponent rather than returning infinity
        n_sub * arg.clamp(-80.0, 80.0).exp()
    }
}

/// Drain current in the linear (low-V_ds) regime, A per unit width in micrometres.
///
/// `I = q mu (W/L) integral of n_s dV`, evaluated in the gradual-channel form
/// `I = (W/L) mu C_ox [(V_gs - V_t) V_ds - V_ds^2/2]` reduced to the charge form used here.
pub fn drain_current_a(ch: &Channel, v_ds_v: Volts) -> Result<f64, ContinuityError> {
    let v_ds_v = v_ds_v.0;
    ch.check()?;
    if !v_ds_v.is_finite() {
        return Err(ContinuityError { what: "the drain voltage must be finite".into() });
    }
    if v_ds_v < 0.0 {
        // reverse bias: the channel is not conducting in this model
        return Ok(0.0);
    }
    // average sheet charge along the channel: the quasi-Fermi level runs from 0 to V_ds
    let q0 = ch.sheet_charge_cm2(0.0);
    let q1 = ch.sheet_charge_cm2(v_ds_v);
    let q_avg = 0.5 * (q0 + q1);
    // I = q mu (W/L) q_avg V_ds, with W = 1 um and L in nm
    let w_cm = 1.0e-4;
    let l_cm = ch.length_nm * 1.0e-7;
    let i = QE * ch.mobility_cm2_vs * (w_cm / l_cm) * q_avg * v_ds_v;
    if !i.is_finite() {
        return Err(ContinuityError { what: "the current is not a number; the bias is outside the model".into() });
    }
    Ok(i)
}

/// Whether the channel is in saturation at this drain bias.
///
/// The gradual-channel model saturates when the drain-end charge has been driven to zero, which is
/// where the quasi-Fermi level reaches the surface potential.
pub fn is_saturated(ch: &Channel, v_ds_v: Volts) -> bool {
    let v_ds_v = v_ds_v.0;
    v_ds_v >= ch.surface_potential_v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ch(psi: f64) -> Channel {
        Channel { surface_potential_v: psi, length_nm: 50.0, mobility_cm2_vs: 400.0, doping_cm3: 1.0e17, temp_k: 300.0 }
    }

    #[test]
    fn the_current_rises_with_gate_voltage() {
        // V-15's first property: I_d increases with V_g
        let low = drain_current_a(&ch(0.3), atto_units::Volts(0.05)).expect("low");
        let mid = drain_current_a(&ch(0.5), atto_units::Volts(0.05)).expect("mid");
        let high = drain_current_a(&ch(0.7), atto_units::Volts(0.05)).expect("high");
        assert!(low < mid && mid < high, "more gate, more current: {low:e} {mid:e} {high:e}");
    }

    #[test]
    fn below_threshold_the_current_is_exponential_in_the_gate() {
        // the subthreshold slope: 60 mV per decade at room temperature, near enough
        let a = drain_current_a(&ch(0.10), atto_units::Volts(0.05)).expect("a");
        let b = drain_current_a(&ch(0.16), atto_units::Volts(0.05)).expect("b");
        let ratio = b / a;
        assert!(ratio > 5.0 && ratio < 50.0, "60 mV gives about a decade: {ratio}");
    }

    #[test]
    fn the_current_rises_then_saturates_with_drain_voltage() {
        let c = ch(0.7);
        let a = drain_current_a(&c, atto_units::Volts(0.1)).expect("a");
        let b = drain_current_a(&c, atto_units::Volts(0.5)).expect("b");
        assert!(b > a, "more drain, more current in the linear regime: {b:e} vs {a:e}");
        assert!(is_saturated(&c, atto_units::Volts(0.8)), "and it is past saturation at 0.8 V");
        assert!(!is_saturated(&c, atto_units::Volts(0.2)), "but not at 0.2 V");
    }

    #[test]
    fn a_reverse_biased_channel_does_not_conduct() {
        assert_eq!(drain_current_a(&ch(0.7), atto_units::Volts(-1.0)).expect("reverse"), 0.0);
    }

    #[test]
    fn a_longer_channel_carries_less_current() {
        let short = drain_current_a(&ch(0.7), atto_units::Volts(0.1)).expect("short");
        let mut long = ch(0.7);
        long.length_nm = 200.0;
        let longer = drain_current_a(&long, atto_units::Volts(0.1)).expect("long");
        assert!((short / longer - 4.0).abs() < 1e-9, "four times the length, a quarter the current: {}", short / longer);
    }

    #[test]
    fn impossible_input_is_refused() {
        let c = ch(0.7);
        let mut bad = c;
        bad.length_nm = 0.0;
        assert!(drain_current_a(&bad, atto_units::Volts(0.1)).is_err(), "no length");
        bad = c;
        bad.temp_k = 0.0;
        assert!(drain_current_a(&bad, atto_units::Volts(0.1)).is_err(), "no temperature");
        assert!(drain_current_a(&c, atto_units::Volts(f64::NAN)).is_err(), "not a number");
    }
}
