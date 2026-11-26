// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Dopant partitioning at an interface (D-15).
//!
//! At a Si/SiO2 boundary a dopant prefers one side. The ratio of the concentrations either side
//! at equilibrium is the segregation coefficient `m = C_Si / C_ox`:
//!
//! * `m > 1` — the dopant prefers silicon; the oxide is a poor sink and the surface layer keeps
//!   its dose.
//! * `m < 1` — the dopant prefers the oxide; a growing oxide **pulls dopant out of the silicon**,
//!   which is why a boron surface layer depletes under oxidation and a phosphorus one does not.
//!
//! Whatever the value of `m`, the total is conserved. That is the property worth asserting: a
//! segregation model that invents or destroys dopant is wrong however plausible its profile looks.
use atto_units::{Nm, PerCm2};

/// Segregation coefficient for boron at a Si/SiO2 interface.
///
/// **Assumed.** Boron prefers the oxide (`m < 1`), which is why a boron surface layer depletes
/// under oxidation. The magnitude is typical of the quoted range and is not fitted to a measured
/// profile; it is registered as `seg.m.boron` so a reader can find it.
pub const BORON_SEGREGATION_M: f64 = 0.3;

/// Where the dopant ended up after an interface was established.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segregation {
    /// Concentration in the silicon at the interface, cm^-3.
    pub c_si_cm3: f64,
    /// Concentration in the oxide at the interface, cm^-3.
    pub c_ox_cm3: f64,
    /// Dose left in the silicon, cm^-2.
    pub dose_si_cm2: f64,
    /// Dose taken up by the oxide, cm^-2.
    pub dose_ox_cm2: f64,
    /// Fraction of the original dose that stayed in the silicon, 0..1.
    pub retained_fraction: f64,
}

impl Segregation {
    /// One line for a report.
    pub fn label(&self) -> String {
        format!(
            "{:.0}% of the dose stays in Si (m = {:.2}, C_Si {:.1e}, C_ox {:.1e})",
            self.retained_fraction * 100.0,
            if self.c_ox_cm3 > 0.0 { self.c_si_cm3 / self.c_ox_cm3 } else { f64::INFINITY },
            self.c_si_cm3,
            self.c_ox_cm3
        )
    }
}

/// Partition a dose between a silicon layer and an oxide layer at equilibrium.
///
/// With `m = C_Si / C_ox`, conservation gives
/// `C_ox = dose / (m t_si + t_ox)` and `C_Si = m C_ox`, where the thicknesses are in cm.
pub fn segregate(
    dose_cm2: PerCm2,
    si_thickness_nm: Nm,
    ox_thickness_nm: Nm,
    segregation_coefficient: f64,
) -> Result<Segregation, String> {
    let dose_cm2 = dose_cm2.0;
    let ox_thickness_nm = ox_thickness_nm.0;
    let si_thickness_nm = si_thickness_nm.0;
    if !dose_cm2.is_finite() || dose_cm2 <= 0.0 {
        return Err("dose must be finite and > 0 cm^-2".into());
    }
    if !si_thickness_nm.is_finite() || si_thickness_nm <= 0.0 {
        return Err("the silicon thickness must be finite and > 0 nm".into());
    }
    if !ox_thickness_nm.is_finite() || ox_thickness_nm < 0.0 {
        return Err("the oxide thickness must be finite and >= 0 nm".into());
    }
    if !segregation_coefficient.is_finite() || segregation_coefficient <= 0.0 {
        return Err("the segregation coefficient must be finite and > 0".into());
    }
    let t_si = si_thickness_nm * 1.0e-7;
    let t_ox = ox_thickness_nm * 1.0e-7;
    let denom = segregation_coefficient * t_si + t_ox;
    if !(denom > 0.0) {
        return Err("the layers have no thickness to hold the dose".into());
    }
    let c_ox = dose_cm2 / denom;
    let c_si = segregation_coefficient * c_ox;
    let dose_si = c_si * t_si;
    let dose_ox = c_ox * t_ox;
    Ok(Segregation {
        c_si_cm3: c_si,
        c_ox_cm3: c_ox,
        dose_si_cm2: dose_si,
        dose_ox_cm2: dose_ox,
        retained_fraction: dose_si / dose_cm2,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_dose_is_conserved_across_the_interface() {
        // D-15's acceptance, asserted as an identity
        let s = segregate(atto_units::PerCm2(1.0e14), atto_units::Nm(100.0), atto_units::Nm(20.0), 0.3).expect("segregation");
        let total = s.dose_si_cm2 + s.dose_ox_cm2;
        assert!((total - 1.0e14).abs() / 1.0e14 < 1e-12, "nothing is created or lost: {total} vs 1e14");
    }

    #[test]
    fn a_coefficient_below_one_depletes_the_silicon() {
        // boron in oxide: m < 1, so the oxide pulls dopant out of the silicon
        let s = segregate(atto_units::PerCm2(1.0e14), atto_units::Nm(100.0), atto_units::Nm(100.0), 0.3).expect("segregation");
        assert!(s.retained_fraction < 0.5, "most of it leaves: {}", s.retained_fraction);
        assert!(s.c_ox_cm3 > s.c_si_cm3, "the oxide is the sink: {} vs {}", s.c_ox_cm3, s.c_si_cm3);
    }

    #[test]
    fn a_coefficient_above_one_keeps_it_in_the_silicon() {
        let s = segregate(atto_units::PerCm2(1.0e14), atto_units::Nm(100.0), atto_units::Nm(20.0), 10.0).expect("segregation");
        assert!(s.retained_fraction > 0.9, "almost all of it stays: {}", s.retained_fraction);
        assert!(s.c_si_cm3 > s.c_ox_cm3, "silicon is the sink");
    }

    #[test]
    fn a_thicker_oxide_takes_more_when_the_coefficient_favours_it() {
        let thin = segregate(atto_units::PerCm2(1.0e14), atto_units::Nm(100.0), atto_units::Nm(10.0), 0.3).expect("thin");
        let thick = segregate(atto_units::PerCm2(1.0e14), atto_units::Nm(100.0), atto_units::Nm(100.0), 0.3).expect("thick");
        assert!(thick.retained_fraction < thin.retained_fraction, "a bigger sink takes more: {} vs {}", thick.retained_fraction, thin.retained_fraction);
    }

    #[test]
    fn no_oxide_keeps_everything() {
        let s = segregate(atto_units::PerCm2(1.0e14), atto_units::Nm(100.0), atto_units::Nm(0.0), 0.3).expect("segregation");
        assert!((s.retained_fraction - 1.0).abs() < 1e-12, "nowhere else to go: {}", s.retained_fraction);
    }

    #[test]
    fn a_recipe_that_cannot_exist_is_refused() {
        assert!(segregate(atto_units::PerCm2(0.0), atto_units::Nm(100.0), atto_units::Nm(20.0), 0.3).is_err(), "no dose");
        assert!(segregate(atto_units::PerCm2(1e14), atto_units::Nm(0.0), atto_units::Nm(20.0), 0.3).is_err(), "no silicon");
        assert!(segregate(atto_units::PerCm2(1e14), atto_units::Nm(100.0), atto_units::Nm(-1.0), 0.3).is_err(), "negative oxide");
        assert!(segregate(atto_units::PerCm2(1e14), atto_units::Nm(100.0), atto_units::Nm(20.0), 0.0).is_err(), "zero coefficient");
        assert!(segregate(atto_units::PerCm2(f64::NAN), atto_units::Nm(100.0), atto_units::Nm(20.0), 0.3).is_err(), "not a number");
    }
}
