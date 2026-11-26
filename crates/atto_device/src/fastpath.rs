// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A one-dimensional analytic fast path (PRD V-05).
//!
//! The full route solves Poisson and both continuity equations on a mesh. That is
//! the right answer and the wrong tool for a slider: a UI that redraws a C-V curve
//! on every keystroke cannot afford a Newton solve per point.
//!
//! This is the depletion approximation in closed form — no mesh, no iteration, no
//! failure modes except the inputs being nonsense. It is *approximate*, and the
//! tests measure how approximate rather than asserting agreement it does not have.

use crate::mos::{EPS0, EPS_OX_R, EPS_SI_R, NI, QE};

/// The fast path's inputs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FastPath {
    /// Oxide thickness, nm.
    pub tox_nm: f64,
    /// Bulk acceptor concentration, cm^-3.
    pub na_cm3: f64,
    /// Gate area, cm^2.
    pub area_cm2: f64,
}

impl FastPath {
    fn valid(&self) -> bool {
        self.tox_nm > 0.0 && self.tox_nm.is_finite()
            && self.na_cm3 > 0.0 && self.na_cm3.is_finite()
            && self.area_cm2 > 0.0 && self.area_cm2.is_finite()
    }

    /// Oxide capacitance per unit area, F/cm^2.
    pub fn cox(&self) -> Option<f64> {
        if !self.valid() {
            return None;
        }
        Some(EPS0 * EPS_OX_R / (self.tox_nm * 1e-7))
    }

    /// Bulk potential, V.
    pub fn phi_f(&self) -> Option<f64> {
        if !self.valid() {
            return None;
        }
        Some(0.025852 * (self.na_cm3 / NI).ln())
    }

    /// Threshold voltage by the depletion approximation, V.
    ///
    /// No flat-band offset and no interface charge: this is the fast estimate, and
    /// the full route is what carries the extras.
    pub fn threshold_v(&self) -> Option<f64> {
        let cox = self.cox()?;
        let phi_f = self.phi_f()?;
        let qb = (4.0 * EPS0 * EPS_SI_R * QE * self.na_cm3 * phi_f).sqrt();
        Some(2.0 * phi_f + qb / cox)
    }

    /// Depletion width at a gate voltage above flat band, nm.
    pub fn depletion_nm(&self, vg: f64) -> Option<f64> {
        let cox = self.cox()?;
        let phi_f = self.phi_f()?;
        if !vg.is_finite() {
            return None;
        }
        // surface potential from the charge balance, clamped to strong inversion
        let mut psi = 0.0f64;
        let target = 2.0 * phi_f;
        let mut lo = 0.0f64;
        let mut hi = target;
        if vg <= 0.0 {
            return Some(0.0);
        }
        // solve Vg = psi + sqrt(2 eps q Na psi)/Cox by bisection: 20 steps is far
        // more than the approximation deserves
        for _ in 0..20 {
            psi = 0.5 * (lo + hi);
            let qb = (2.0 * EPS0 * EPS_SI_R * QE * self.na_cm3 * psi.max(0.0)).sqrt();
            let f = psi + qb / cox - vg;
            if f > 0.0 { hi = psi } else { lo = psi }
        }
        let w_cm = (2.0 * EPS0 * EPS_SI_R * psi.max(0.0) / (QE * self.na_cm3)).sqrt();
        Some(w_cm * 1e7)
    }

    /// Capacitance per unit area at a gate voltage, F/cm^2.
    pub fn cap_per_cm2(&self, vg: f64) -> Option<f64> {
        let cox = self.cox()?;
        let w_nm = self.depletion_nm(vg)?;
        if w_nm <= 0.0 {
            return Some(cox);
        }
        let w_cm = w_nm * 1e-7;
        let cs = EPS0 * EPS_SI_R / w_cm;
        Some(cox * cs / (cox + cs))
    }

    /// Total capacitance at a gate voltage, F.
    pub fn capacitance(&self, vg: f64) -> Option<f64> {
        if !self.valid() {
            return None;
        }
        Some(self.cap_per_cm2(vg)? * self.area_cm2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn typical() -> FastPath {
        FastPath { tox_nm: 20.0, na_cm3: 1.0e16, area_cm2: 1.0e-4 }
    }

    #[test]
    fn the_oxide_capacitance_is_the_textbook_value() {
        let f = typical();
        let cox = f.cox().unwrap();
        // 3.9 * 8.854e-14 / 20 nm = 1.73e-7 F/cm^2
        assert!((cox - 1.727e-7).abs() / 1.727e-7 < 0.01, "cox {cox:.3e}");
    }

    #[test]
    fn the_threshold_is_in_the_right_place() {
        let f = typical();
        let vt = f.threshold_v().unwrap();
        assert!(vt > 0.3 && vt < 2.0, "a 1e16 substrate at 20 nm gives {vt:.2} V");
        // A range is not a check: a mutation that scaled the bulk charge by 1.5
        // stayed inside it and the suite passed. The value is pinned instead.
        assert!((vt - 0.996).abs() < 0.02, "the threshold moved: {vt:.4} V");
    }

    #[test]
    fn the_depletion_width_is_accurate_not_merely_in_range() {
        // The bisection had 20 steps; truncating it to 3 passed every test, because
        // nothing asserted the width's accuracy — only that capacitance falls. A
        // coarse solve gives a width that is in range and wrong.
        let f = typical();
        let w = f.depletion_nm(0.5).expect("a width");
        assert!((w - 201.3).abs() < 2.0, "the depletion width at 0.5 V is {w:.1} nm");
    }

    #[test]
    fn the_capacitance_falls_from_cox_and_saturates() {
        let f = typical();
        let low = f.cap_per_cm2(0.0).unwrap();
        // mid must still be in depletion: past strong inversion the width pins at
        // its maximum and the capacitance stops falling, which is the floor below
        let mid = f.cap_per_cm2(0.5).unwrap();
        let high = f.cap_per_cm2(10.0).unwrap();
        assert!((low - f.cox().unwrap()).abs() / low < 1e-9, "flat band is the oxide cap");
        assert!(mid < low, "depletion lowers it");
        assert!(high < mid, "and strong inversion lowers it further: {mid:.3e} vs {high:.3e}");
        // it must not fall below a sane floor
        assert!(high > 0.05 * low, "the floor is not zero: {high:.3e}");
        // and it stops falling once the depletion width saturates
        let later = f.cap_per_cm2(20.0).unwrap();
        assert!((later - high).abs() / high < 1e-9, "the floor is flat: {later:.3e}");
    }

    #[test]
    fn a_thicker_oxide_gives_a_lower_capacitance() {
        let thin = FastPath { tox_nm: 10.0, ..typical() };
        let thick = FastPath { tox_nm: 40.0, ..typical() };
        assert!(thin.cox().unwrap() > thick.cox().unwrap());
        assert!(thin.threshold_v().unwrap() < thick.threshold_v().unwrap());
    }

    #[test]
    fn bad_input_is_refused() {
        assert!(FastPath { tox_nm: 0.0, ..typical() }.cox().is_none());
        assert!(FastPath { na_cm3: -1.0, ..typical() }.threshold_v().is_none());
        assert!(FastPath { area_cm2: 0.0, ..typical() }.capacitance(1.0).is_none());
        assert!(typical().depletion_nm(f64::NAN).is_none());
    }
}
