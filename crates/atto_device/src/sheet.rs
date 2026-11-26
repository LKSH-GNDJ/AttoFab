// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Sheet resistance and the four-point probe (M-07).
//!
//! A four-point probe is the cheapest electrical measurement in a fab and the one most often
//! compared against a simulation. It reports ohms per square, which for a doped layer is
//! `1 / (q mu N t)` — the conductivity integrated through the film, with no thickness in it.
//!
//! The probe geometry enters as a constant: for a thin layer with probe spacing far larger
//! than the film thickness, `R_s = (pi / ln 2) (V / I)`, whatever the spacing is. That is the
//! only geometric fact this module uses, and it is exact for the thin-layer limit.

use atto_units::{Amps, PerCm2, Volts};

/// The elementary charge, coulombs — the same constant the MOS electrostatics uses.
pub use crate::mos::QE;

/// A typical drift mobility for a moderately doped silicon layer, cm^2/Vs.
///
/// **Assumed.** A real mobility depends on the species, the concentration and the temperature,
/// and this module does not model that. Registered as `sheet.mu.default` so a reader can find
/// it and see what it is.
pub const DEFAULT_MOBILITY_CM2_VS: f64 = 500.0;

/// `pi / ln 2` — the four-point probe's geometric factor for a thin layer.
///
/// Independent of probe spacing and of film thickness in the thin-layer limit, which is why a
/// four-point measurement needs no length at all to be converted to ohms per square.
pub const FOUR_POINT_FACTOR: f64 = std::f64::consts::PI / std::f64::consts::LN_2;

/// Sheet resistance of a uniformly doped layer, ohms per square.
///
/// `R_s = 1 / (q mu N t)`. The areal dose `N t` is what enters; the thickness is already inside
/// it, so a 100 nm layer at 1e17 cm^-3 and a 200 nm layer at 5e16 cm^-3 give the same answer.
pub fn sheet_resistance(active_dose: PerCm2, mobility_cm2_vs: f64) -> Result<f64, String> {
    if !active_dose.0.is_finite() || active_dose.0 <= 0.0 {
        return Err("active dose must be finite and > 0 cm^-2".into());
    }
    if !mobility_cm2_vs.is_finite() || mobility_cm2_vs <= 0.0 {
        return Err("mobility must be finite and > 0 cm^2/Vs".into());
    }
    Ok(1.0 / (QE * mobility_cm2_vs * active_dose.0))
}

/// The voltage a four-point probe reads for a sheet resistance at a given current, volts.
///
/// `V = I R_s / (pi / ln 2)`.
pub fn four_point_voltage(sheet_ohm_sq: f64, current_a: Amps) -> Result<f64, String> {
    let current_a = current_a.0;
    if !sheet_ohm_sq.is_finite() || sheet_ohm_sq <= 0.0 {
        return Err("sheet resistance must be finite and > 0 ohm/sq".into());
    }
    if !current_a.is_finite() || current_a <= 0.0 {
        return Err("probe current must be finite and > 0 A".into());
    }
    Ok(current_a * sheet_ohm_sq / FOUR_POINT_FACTOR)
}

/// Sheet resistance from a measured four-point voltage, ohms per square.
///
/// The inverse of [`four_point_voltage`], so a measured V and I can be compared with a
/// simulated profile without either side owning the conversion.
pub fn from_four_point(voltage_v: Volts, current_a: Amps) -> Result<f64, String> {
    let current_a = current_a.0;
    let voltage_v = voltage_v.0;
    if !voltage_v.is_finite() || voltage_v <= 0.0 {
        return Err("probe voltage must be finite and > 0 V".into());
    }
    if !current_a.is_finite() || current_a <= 0.0 {
        return Err("probe current must be finite and > 0 A".into());
    }
    Ok(FOUR_POINT_FACTOR * voltage_v / current_a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_known_uniform_layer_matches_the_analytic_value() {
        // 1e17 cm^-3 over 100 nm is 1e12 cm^-2; at 500 cm^2/Vs the analytic sheet resistance
        // is 1 / (1.602176634e-19 x 500 x 1e12) = 12 484.0 ohm/sq
        let dose = PerCm2::new(1.0e17 * 1.0e-5);
        let rs = sheet_resistance(dose, 500.0).expect("a valid layer");
        let analytic = 1.0 / (1.602_176_634e-19 * 500.0 * 1.0e12);
        assert!((rs - analytic).abs() / analytic < 0.05, "within 5% of analytic: {rs} vs {analytic}");
        assert!((rs - 12_484.0).abs() / 12_484.0 < 0.05, "and near the handbook number: {rs}");
    }

    #[test]
    fn thickness_and_concentration_trade_off_at_equal_dose() {
        // the same areal dose twice, reached two ways
        let a = sheet_resistance(PerCm2::new(1.0e12), 500.0).expect("a");
        let b = sheet_resistance(PerCm2::new(1.0e12), 500.0).expect("b");
        assert!((a - b).abs() < 1e-12, "the dose is the only thing that enters");
    }

    #[test]
    fn doubling_the_mobility_halves_the_resistance() {
        let a = sheet_resistance(PerCm2::new(1.0e12), 250.0).expect("a");
        let b = sheet_resistance(PerCm2::new(1.0e12), 500.0).expect("b");
        assert!((a / b - 2.0).abs() < 1e-12, "R_s falls as mobility rises");
    }

    #[test]
    fn the_probe_round_trips_through_voltage() {
        let rs = sheet_resistance(PerCm2::new(1.0e12), 500.0).expect("a layer");
        let v = four_point_voltage(rs, atto_units::Amps(1.0e-3)).expect("a voltage");
        let back = from_four_point(atto_units::Volts(v), atto_units::Amps(1.0e-3)).expect("back");
        assert!((back - rs).abs() / rs < 1e-12, "V and R_s are inverses: {back} vs {rs}");
    }

    #[test]
    fn the_geometric_factor_is_the_textbook_number() {
        assert!((FOUR_POINT_FACTOR - 4.532_360).abs() < 1e-5, "pi/ln2 is about 4.5324");
    }

    #[test]
    fn a_layer_that_cannot_exist_is_refused() {
        assert!(sheet_resistance(PerCm2::new(0.0), 500.0).is_err(), "no dose");
        assert!(sheet_resistance(PerCm2::new(1.0e12), 0.0).is_err(), "no mobility");
        assert!(sheet_resistance(PerCm2::new(f64::NAN), 500.0).is_err(), "not a number");
        assert!(four_point_voltage(0.0, atto_units::Amps(1e-3)).is_err(), "no resistance");
        assert!(four_point_voltage(100.0, atto_units::Amps(0.0)).is_err(), "no current");
        assert!(from_four_point(atto_units::Volts(-1.0), atto_units::Amps(1e-3)).is_err(), "negative voltage");
    }
}

#[cfg(test)]
mod activation_coupling_tests {
    use super::*;
    use atto_dopant::activation::activate;
    use atto_params::Species;
    use atto_units::{Kelvin, PerCm3};

    #[test]
    fn sheet_resistance_sees_the_active_dose_not_the_total() {
        // D-14's second half: a 1e20 implant at 1000 C, and R_s from what is actually active.
        // The reference is the analytic 1/(q mu N_active t) with the same constants, so this
        // tests the chain — activation feeding sheet resistance — not a published table.
        let total = PerCm3(1.0e20);
        let t_nm = 100.0;
        let a = activate(total, Species::Boron, Kelvin::new(1273.15)).expect("an activation");
        assert!(a.active.0 < total.0, "active below total: {} vs {}", a.active.0, total.0);

        let dose_active = PerCm2::new(a.active.0 * t_nm * 1.0e-7);
        let dose_total = PerCm2::new(total.0 * t_nm * 1.0e-7);
        let rs_active = sheet_resistance(dose_active, 500.0).expect("active");
        let rs_total = sheet_resistance(dose_total, 500.0).expect("total");

        assert!(rs_active > rs_total, "fewer carriers means more resistance: {rs_active} vs {rs_total}");

        let analytic = 1.0 / (1.602_176_634e-19 * 500.0 * a.active.0 * t_nm * 1.0e-7);
        assert!((rs_active - analytic).abs() / analytic < 0.10, "within 10% of the analytic reference: {rs_active} vs {analytic}");
    }
}
