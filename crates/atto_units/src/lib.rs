// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Dimensioned scalars. Every physical quantity crossing a module boundary
//! carries its unit in its type, so the class of defect the the earlier build audit found
//! (a cm^2/nm^2 factor, a metres/nm radius, a hardcoded database unit) cannot
//! be expressed.


use std::fmt;
use std::ops::{Add, Div, Mul, Sub};

macro_rules! unit {
    ($name:ident, $unit:literal, $doc:literal) => {
        #[doc = $doc]
        #[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
        /// The value in this quantity's own unit.
        pub struct $name(pub f64);

        impl $name {
            /// The zero of this quantity.
            pub const ZERO: $name = $name(0.0);
            /// Build from the value in this quantity's own unit.
            pub fn new(v: f64) -> Self { $name(v) }
            /// The value in this quantity's own unit.
            pub fn value(self) -> f64 { self.0 }
            /// Whether the value is finite, for validation at a boundary.
            pub fn is_finite(self) -> bool { self.0.is_finite() }
            /// The smaller of two values.
            pub fn min(self, o: Self) -> Self { $name(self.0.min(o.0)) }
            /// The larger of two values.
            pub fn max(self, o: Self) -> Self { $name(self.0.max(o.0)) }
            /// The magnitude, discarding the sign.
            pub fn abs(self) -> Self { $name(self.0.abs()) }
        }
        impl Add for $name { type Output = Self; fn add(self, o: Self) -> Self { $name(self.0 + o.0) } }
        impl Sub for $name { type Output = Self; fn sub(self, o: Self) -> Self { $name(self.0 - o.0) } }
        impl Mul<f64> for $name { type Output = Self; fn mul(self, k: f64) -> Self { $name(self.0 * k) } }
        impl Div<f64> for $name { type Output = Self; fn div(self, k: f64) -> Self { $name(self.0 / k) } }
        impl Div for $name { type Output = f64; fn div(self, o: Self) -> f64 { self.0 / o.0 } }
        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                write!(f, "{} {}", self.0, $unit)
            }
        }
    };
}

unit!(Nm, "nm", "Length in nanometres.");
unit!(Um, "um", "Length in micrometres, for a wafer-scale bow.");
unit!(Mm, "mm", "Length in millimetres, for a wafer diameter.");
unit!(KeV, "keV", "Ion energy in kiloelectronvolts.");
unit!(PerCm3, "cm^-3", "Number density per cubic centimetre.");
unit!(PerCm2, "cm^-2", "Areal density (implant dose) per square centimetre.");
unit!(Watts, "W", "RF power.");
unit!(MTorr, "mTorr", "Chamber pressure.");
unit!(MPa, "MPa", "Stress.");
unit!(Minutes, "min", "Process time in minutes.");
unit!(Seconds, "s", "Time in seconds, for a lifetime or a rate constant.");
unit!(Kelvin, "K", "Absolute temperature.");
unit!(Volts, "V", "Electric potential.");
unit!(Celsius, "C", "Temperature in degrees Celsius.");
unit!(Degrees, "deg", "An angle in degrees.");
unit!(Amps, "A", "Electric current.");
unit!(VPerCm, "V/cm", "Electric field.");
unit!(FPerCm2, "F/cm^2", "Capacitance per unit area.");
unit!(NmPerMin, "nm/min", "A rate in nanometres per minute.");
unit!(MjPerCm2, "mJ/cm^2", "Exposure dose.");
unit!(Metres, "m", "Length in metres, for a curvature radius.");
unit!(Fraction, "1", "A dimensionless fraction between 0 and 1.");
unit!(KPa, "kPa", "Pressure in kilopascals, for a polishing down-force.");
unit!(MetresPerSecond, "m/s", "Relative pad-wafer velocity in metres per second.");

/// 1 cm = 1e7 nm, so 1 cm^2 = 1e14 nm^2 and 1 cm^3 = 1e21 nm^3.
pub const NM2_PER_CM2: f64 = 1.0e14;
pub const NM3_PER_CM3: f64 = 1.0e21;

/// A Gaussian has this factor in its normalisation.
pub const SQRT_2PI: f64 = 2.506_628_274_631_000_7;

impl Nm {
    /// Convert an area in nm^2 to cm^2.
    pub fn nm2_to_cm2(a_nm2: f64) -> f64 { a_nm2 / NM2_PER_CM2 }
    /// Convert a volume in nm^3 to cm^3.
    pub fn nm3_to_cm3(v_nm3: f64) -> f64 { v_nm3 / NM3_PER_CM3 }
}

impl PerCm3 {
    /// Integrate a concentration (cm^-3) over a length (nm) to an areal dose (cm^-2).
    /// This is the conversion the earlier build performed inline as `* 1e-7`; here it is named once.
    pub fn times_length(self, d: Nm) -> PerCm2 { PerCm2(self.0 * d.0 * 1.0e-7) }
}

impl PerCm2 {
    /// Areal dose to the peak of a Gaussian of the given straggle.
    pub fn gaussian_peak(self, sigma: Nm) -> PerCm3 {
        PerCm3(self.0 / (SQRT_2PI * sigma.0 * 1.0e-7))
    }
}

impl Kelvin {
    /// Room temperature, the reference for anneal parameters.
    pub const ROOM: Kelvin = Kelvin(293.15);
    pub fn from_celsius(c: f64) -> Self { Kelvin(c + 273.15) }
    /// The temperature in degrees Celsius.
    pub fn celsius(self) -> f64 { self.0 - 273.15 }
}

/// Boltzmann constant in eV/K.
pub const KB_EV: f64 = 8.617_333_262e-5;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn area_conversion_round_trips() {
        assert!((Nm::nm2_to_cm2(1.0e14) - 1.0).abs() < 1e-12);
        assert!((Nm::nm3_to_cm3(1.0e21) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn dose_from_concentration_and_length() {
        // 1e18 cm^-3 over 100 nm = 1e13 cm^-2.
        let d = PerCm3(1.0e18).times_length(Nm(100.0));
        assert!((d.value() - 1.0e13).abs() < 1e6, "{d}");
    }

    #[test]
    fn gaussian_peak_matches_hand_calculation() {
        // 1e14 cm^-2 with sigma 55.7 nm -> 7.2e18 cm^-3 (the earlier build default doping case).
        let p = PerCm2(1.0e14).gaussian_peak(Nm(55.7));
        assert!((p.value() - 7.17e18).abs() / 7.17e18 < 0.01, "{p}");
    }

    #[test]
    fn temperature_conversion() {
        assert!((Kelvin::from_celsius(1000.0).celsius() - 1000.0).abs() < 1e-9);
        assert!((Kelvin::ROOM.celsius() - 20.0).abs() < 1e-9);
    }

    #[test]
    fn units_do_not_mix() {
        // The point of the newtype: this does not compile.
        // let _ = Nm(1.0) + PerCm3(1.0);
        let a = Nm(2.0) + Nm(3.0);
        assert_eq!(a.value(), 5.0);
    }
}
