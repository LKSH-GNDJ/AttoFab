// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Long-channel MOS electrostatics (build PRD V-03).
//!
//! The depletion approximation: enough to give a threshold voltage from the
//! process state, and enough to prove that the cross-section carries the doping,
//! the oxide and the interface charge a device solve needs.

use atto_units::{Nm, PerCm2};

use crate::cross_section::{CrossSection, CrossSectionError};

/// Permittivity of free space, F/cm.
pub const EPS0: f64 = 8.854_187_812_8e-14;
/// Relative permittivity of silicon and of thermal oxide.
pub const EPS_SI_R: f64 = 11.7;
pub const EPS_OX_R: f64 = 3.9;
/// Elementary charge, C.
pub const QE: f64 = 1.602_176_634e-19;
/// Boltzmann constant, J/K.
pub const KB: f64 = 1.380_649e-23;
/// Intrinsic carrier concentration at 300 K, cm^-3.
pub const NI: f64 = 1.0e10;

#[derive(Debug, Clone, Copy)]
/// The gate stack and its interface charge.
pub struct MosParams {
    pub temp_k: f64,
    /// Metal-semiconductor work function difference, V. n+ polysilicon on
    /// p-silicon is about -0.55 V.
    pub phi_ms: f64,
    /// Interface state density, cm^-2 eV^-1. Zero is a clean interface and leaves the
    /// threshold where it was before this field existed.
    pub dit_cm2ev: f64,
}

impl Default for MosParams {
    fn default() -> Self {
        MosParams { temp_k: 300.15, phi_ms: -0.55, dit_cm2ev: 0.0 }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// The threshold voltage and what went into it.
pub struct ThresholdResult {
    /// Threshold voltage, V.
    pub vt: f64,
    /// Flat-band voltage, V.
    pub vfb: f64,
    /// Bulk Fermi potential, V (magnitude).
    pub phi_f: f64,
    /// Areal oxide capacitance, F/cm^2.
    pub cox: f64,
    /// Depletion charge at the onset of strong inversion, C/cm^2.
    pub qdep: f64,
    /// Oxide thickness used, nm.
    pub tox_nm: f64,
    /// Bulk doping used, cm^-3 (signed, donors positive).
    pub n_bulk: f64,
}

impl ThresholdResult {
    /// A one-line description of the result.
    pub fn label(&self) -> String {
        format!("Vt {:.3} V (Vfb {:.3}, phi_f {:.3}, tox {:.0} nm, N {:.1e})", self.vt, self.vfb, self.phi_f, self.tox_nm, self.n_bulk)
    }
}

/// Threshold voltage from the three numbers a long-channel model needs, with no
/// cross-section involved. This is what a design-of-experiments sweep calls.
pub fn threshold_from(tox: Nm, n_bulk: f64, qf: PerCm2, p: &MosParams) -> Option<ThresholdResult> {
    // A-03: the two arguments that carry a unit are dimensioned. A caller that
    // passed micrometres for nanometres now fails to compile rather than
    // producing a plausible wrong number.
    let (tox_nm, qf_cm2) = (tox.value(), qf.value());
    if !(tox_nm > 0.0) || !tox_nm.is_finite() || !n_bulk.is_finite() || n_bulk == 0.0 {
        return None;
    }
    let n_abs = n_bulk.abs();
    let vt_thermal = KB * p.temp_k / QE;
    // Heavy doping narrows the gap, which raises the effective intrinsic concentration and so
    // lowers phi_f and the threshold. At light doping this is exactly NI and nothing moves.
    let ni_eff = crate::heavy::ni_effective_cm3(atto_units::PerCm3(NI), atto_units::PerCm3(n_abs), atto_units::Kelvin(p.temp_k)).unwrap_or(NI);
    let phi_f = vt_thermal * (n_abs / ni_eff).ln();
    if !(phi_f > 0.0) {
        return None;
    }
    let cox = EPS_OX_R * EPS0 / (tox_nm * 1.0e-7);
    let qdep = (2.0 * EPS_SI_R * EPS0 * QE * n_abs * 2.0 * phi_f).max(0.0).sqrt();
    // Fixed charge and interface states both act through the oxide capacitance, so both are
    // computed in one place. The inline `qf * QE / cox` that used to be here is now the
    // module's `vfb_shift_v`, and the stretch-out V-13 asks about is added to the threshold.
    let effect = crate::interface_effect::interface_effect(atto_units::PerCm2(qf_cm2), p.dit_cm2ev, atto_units::FPerCm2(cox), atto_units::Volts(2.0 * phi_f)).ok()?;
    let vfb = p.phi_ms + effect.vfb_shift_v;
    let vt = vfb + 2.0 * phi_f + qdep / cox + effect.stretch_out_v;
    Some(ThresholdResult { vt, vfb, phi_f, cox, qdep, tox_nm, n_bulk })
}

/// Threshold voltage of the MOS structure in one column.
pub fn threshold_voltage(cs: &CrossSection, ix: usize, p: &MosParams) -> Result<ThresholdResult, CrossSectionError> {
    if ix >= cs.nx {
        return Err(CrossSectionError::NoOxide(ix));
    }
    let tox_nm = cs.oxide_nm[ix];
    if !(tox_nm > 0.0) {
        return Err(CrossSectionError::NoOxide(ix));
    }
    let n_bulk = cs.bulk_doping(ix).ok_or(CrossSectionError::NoOxide(ix))?;
    threshold_from(Nm::new(tox_nm), n_bulk, PerCm2::new(cs.qf_cm2[ix]), p).ok_or(CrossSectionError::NoOxide(ix))
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;
    use atto_dopant::{DopantField, DopantRecord, Type};
    use atto_geom::column::{set_interval, solid_top, ColumnGrid, ProcessWindow};
    use atto_geom::svo::{Material, Wafer};
    use atto_params::Species;
    use atto_process::oxidation::InterfaceCharge;

    fn build(n_bulk: f64, tox_nm: f64, qf: f64) -> CrossSection {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] };
        let g = ColumnGrid::sample(&w, &win).unwrap();
        let mut f = DopantField::new(atto_units::Nm::new(16.0), w.bounds.min);
        f.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, n_bulk));
        let mut w2 = w.clone();
        let mut gg = g.clone();
        for iy in 0..gg.ny {
            for ix in 0..gg.nx {
                let k = gg.idx(ix, iy);
                let top = solid_top(&gg.cols[k]).unwrap().0;
                set_interval(&mut gg.cols[k], top, top + tox_nm, Material::SiliconDioxide);
            }
        }
        gg.write_back(&g, &mut w2);
        CrossSection::build(&w2, &f, &win, InterfaceCharge { qf_cm2: qf, dit_cm2ev: 1e10 }).unwrap()
    }

    #[test]
    fn matches_a_hand_calculation() {
        // N_A = 1e15, t_ox = 32 nm (two 16 nm cells, so no quantisation error),
        // no fixed charge, n+ polysilicon gate.
        let cs = build(1.0e15, 32.0, 0.0);
        let r = threshold_voltage(&cs, 0, &MosParams::default()).unwrap();
        assert!((r.tox_nm - 32.0).abs() < 1e-9, "tox {:.1}", r.tox_nm);
        // phi_f = kT/q ln(N_A/n_i) = 0.2978 V
        // Cox   = 3.9 eps0 / 32 nm = 1.0791e-7 F/cm^2
        // Qdep  = sqrt(2 eps_s q N_A 2 phi_f) = 1.406e-8 C/cm^2
        // Vt    = -0.55 + 2 phi_f + Qdep/Cox = 0.1759 V
        assert!((r.phi_f - 0.2978).abs() < 0.002, "phi_f {:.4}", r.phi_f);
        assert!((r.cox - 1.0791e-7).abs() / 1.0791e-7 < 0.01, "Cox {:.4e}", r.cox);
        assert!((r.qdep - 1.4060e-8).abs() / 1.4060e-8 < 0.01, "Qdep {:.4e}", r.qdep);
        assert!((r.vt - 0.1759).abs() < 0.005, "Vt {:.4}", r.vt);
    }

    #[test]
    fn positive_fixed_charge_lowers_the_threshold() {
        let a = threshold_voltage(&build(1.0e15, 20.0, 0.0), 0, &MosParams::default()).unwrap();
        let b = threshold_voltage(&build(1.0e15, 20.0, 1.0e11), 0, &MosParams::default()).unwrap();
        assert!(b.vt < a.vt, "Qf 1e11 lowers Vt: {:.3} vs {:.3}", b.vt, a.vt);
        // the shift is Qf/Cox
        let expect = 1.0e11 * QE / a.cox;
        assert!((a.vt - b.vt - expect).abs() < 0.005, "shift {:.3} vs {:.3}", a.vt - b.vt, expect);
    }

    #[test]
    fn a_thicker_oxide_raises_the_threshold() {
        let a = threshold_voltage(&build(1.0e15, 10.0, 0.0), 0, &MosParams::default()).unwrap();
        let b = threshold_voltage(&build(1.0e15, 40.0, 0.0), 0, &MosParams::default()).unwrap();
        assert!(b.vt > a.vt, "thicker oxide raises Vt: {:.3} vs {:.3}", b.vt, a.vt);
        assert!(b.cox < a.cox);
    }

    #[test]
    fn heavier_doping_raises_the_threshold() {
        let a = threshold_voltage(&build(1.0e15, 20.0, 0.0), 0, &MosParams::default()).unwrap();
        let b = threshold_voltage(&build(1.0e17, 20.0, 0.0), 0, &MosParams::default()).unwrap();
        assert!(b.vt > a.vt, "heavier doping raises Vt: {:.3} vs {:.3}", b.vt, a.vt);
        assert!(b.phi_f > a.phi_f);
    }

    #[test]
    fn the_standalone_threshold_agrees_with_the_column_one() {
        let cs = build(1.0e15, 32.0, 0.0);
        let a = threshold_voltage(&cs, 0, &MosParams::default()).unwrap();
        let b = threshold_from(Nm::new(a.tox_nm), a.n_bulk, PerCm2::new(0.0), &MosParams::default()).unwrap();
        assert!((a.vt - b.vt).abs() < 1e-12);
    }

    #[test]
    fn the_standalone_threshold_rejects_impossible_inputs() {
        assert!(threshold_from(Nm::new(0.0), -1.0e15, PerCm2::new(0.0), &MosParams::default()).is_none());
        assert!(threshold_from(Nm::new(f64::NAN), -1.0e15, PerCm2::new(0.0), &MosParams::default()).is_none());
        assert!(threshold_from(Nm::new(20.0), 0.0, PerCm2::new(0.0), &MosParams::default()).is_none());
    }

    #[test]
    fn no_oxide_is_an_error_not_a_number() {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] };
        let mut f = DopantField::new(atto_units::Nm::new(16.0), w.bounds.min);
        f.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e15));
        let cs = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        assert!(threshold_voltage(&cs, 0, &MosParams::default()).is_err());
    }
}

#[cfg(test)]
mod interface_wiring_tests {
    use super::*;

    #[test]
    fn a_clean_interface_leaves_the_threshold_where_it_was() {
        // the inline `qf * QE / cox` this replaced must give the same answer at Dit = 0
        let p = MosParams { temp_k: 300.15, phi_ms: -0.55, dit_cm2ev: 0.0 };
        let cox = EPS_OX_R * EPS0 / (96.0 * 1.0e-7);
        let qf = 1.0e11;
        let effect = crate::interface_effect::interface_effect(atto_units::PerCm2(qf), p.dit_cm2ev, atto_units::FPerCm2(cox), atto_units::Volts(0.6)).expect("effect");
        let old_style = -qf * QE / cox;
        assert!((effect.vfb_shift_v - old_style).abs() < 1e-15, "the module reproduces the inline term: {} vs {}", effect.vfb_shift_v, old_style);
        assert_eq!(effect.stretch_out_v, 0.0, "and a clean interface does not stretch");
    }

    #[test]
    fn interface_states_raise_the_threshold() {
        // a stretch-out means the gate needs more volts for the same surface potential
        let cox = EPS_OX_R * EPS0 / (96.0 * 1.0e-7);
        let clean = crate::interface_effect::interface_effect(atto_units::PerCm2(0.0), 0.0, atto_units::FPerCm2(cox), atto_units::Volts(0.6)).expect("clean");
        let dirty = crate::interface_effect::interface_effect(atto_units::PerCm2(0.0), 5.0e11, atto_units::FPerCm2(cox), atto_units::Volts(0.6)).expect("dirty");
        assert!(dirty.stretch_out_v > 0.0, "interface states stretch the curve: {}", dirty.stretch_out_v);
        assert!(dirty.stretch_out_v > clean.stretch_out_v);
    }
}

#[cfg(test)]
mod heavy_wiring_tests {
    use super::*;

    #[test]
    fn light_doping_keeps_the_intrinsic_concentration() {
        // below the narrowing threshold ni_eff is exactly NI, so phi_f is unchanged
        let ni = crate::heavy::ni_effective_cm3(atto_units::PerCm3(NI), atto_units::PerCm3(1.0e15), atto_units::Kelvin(300.15)).expect("ni");
        assert!((ni - NI).abs() / NI < 1e-12, "light doping leaves n_i alone: {ni}");
    }

    #[test]
    fn narrowing_lowers_the_fermi_potential_at_fixed_doping() {
        // Compare the SAME doping with and without narrowing. Comparing 1e17 against 1e20 does
        // not show the effect — a heavier well has a larger phi_f anyway, and narrowing only
        // reduces it relative to the same well unnarrowed. This is the third test in this
        // session to make exactly that mistake; the pattern is worth naming.
        let vt = KB * 300.15 / QE;
        let n = 1.0e20f64;
        let without = vt * (n / NI).ln();
        let ni_eff = crate::heavy::ni_effective_cm3(atto_units::PerCm3(NI), atto_units::PerCm3(n), atto_units::Kelvin(300.15)).expect("ni");
        let with = vt * (n / ni_eff).ln();
        assert!(with < without, "narrowing lowers phi_f: {with} vs {without}");
        assert!(with > 0.0, "and it stays positive: {with}");
        assert!(without - with > 0.05, "the shift is large enough to matter: {}", without - with);
    }
}
