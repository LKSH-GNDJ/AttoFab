// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Heavy doping: bandgap narrowing and degeneracy statistics (V-10).
//!
//! Above roughly 1e17 cm^-3 two things stop being negligible. The band gap shrinks because the
//! dopant states perturb the lattice (bandgap narrowing), which raises the effective intrinsic
//! carrier concentration; and the carriers become degenerate, so Boltzmann statistics
//! over-predicts how many of them there are.
//!
//! Both push the same way on a C-V: the built-in potential falls, so the flat-band voltage
//! shifts. A simulator that ignores them reports a C-V that is too far to the left, and the
//! error grows exactly where modern devices live.

use atto_units::{Kelvin, PerCm3};

/// Boltzmann constant, eV/K.
const KB_EV: f64 = 8.617_333_262e-5;

/// Why a heavy-doping calculation could not run.
#[derive(Debug, Clone, PartialEq)]
pub enum HeavyError {
    /// A concentration or temperature was not usable.
    BadInput(&'static str),
}

impl std::fmt::Display for HeavyError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            HeavyError::BadInput(w) => write!(f, "heavy-doping model: {w}"),
        }
    }
}

/// Reference doping for the narrowing expression, cm^-3.
pub const N_REF_CM3: f64 = 1.0e17;
/// Slotboom's reference energy, eV.
pub const E_REF_EV: f64 = 9.0e-3;

/// Bandgap narrowing, eV (Slotboom–de Graaff form).
///
/// `dEg = E_ref [ ln(u) + sqrt(ln(u)^2 + 0.5) ]` with `u = N / N_ref`. It is zero at `N_ref`
/// and rises roughly logarithmically above it, which is the shape every measurement shows.
pub fn band_gap_narrowing_ev(doping: PerCm3) -> Result<f64, HeavyError> {
    let doping_cm3 = doping.0;
    if !doping_cm3.is_finite() || doping_cm3 <= 0.0 {
        return Err(HeavyError::BadInput("doping must be finite and > 0 cm^-3"));
    }
    if doping_cm3 <= N_REF_CM3 {
        return Ok(0.0);
    }
    let ln_u = (doping_cm3 / N_REF_CM3).ln();
    Ok(E_REF_EV * (ln_u + (ln_u * ln_u + 0.5).sqrt()))
}

/// Effective intrinsic carrier concentration with narrowing, cm^-3.
///
/// `n_i,eff^2 = n_i^2 exp(dEg / kT)`. The narrowing raises it, which is why a heavily doped
/// junction has a smaller built-in potential than the intrinsic value would give.
pub fn ni_effective_cm3(ni: PerCm3, doping: PerCm3, temp: Kelvin) -> Result<f64, HeavyError> {
    let (ni_cm3, doping_cm3, temp_k) = (ni.0, doping.0, temp.0);
    if !ni_cm3.is_finite() || ni_cm3 <= 0.0 {
        return Err(HeavyError::BadInput("intrinsic concentration must be finite and > 0"));
    }
    if !temp_k.is_finite() || temp_k <= 0.0 {
        return Err(HeavyError::BadInput("temperature must be positive absolute"));
    }
    let d_eg = band_gap_narrowing_ev(PerCm3(doping_cm3))?;
    Ok(ni_cm3 * (d_eg / (KB_EV * temp_k)).exp())
}

/// The degeneracy correction to Boltzmann statistics.
///
/// Returns `n / Nc` given `eta = (Ef - Ec) / kT`. Boltzmann gives `exp(eta)`; the
/// Joyce–Dixon approximation adds the first correction, and for a degenerate electron gas it
/// returns *less* than Boltzmann for the same Fermi level — which is the point.
pub fn reduced_concentration(eta: f64) -> Result<f64, HeavyError> {
    if !eta.is_finite() {
        return Err(HeavyError::BadInput("reduced Fermi level must be finite"));
    }
    let u = eta.exp();
    if !u.is_finite() {
        return Err(HeavyError::BadInput("reduced Fermi level is too large to exponentiate"));
    }
    // Joyce-Dixon: n/Nc = exp(eta) / (1 + 0.25 exp(eta) / sqrt(2) ...) — first order:
    // n/Nc ~ exp(eta) * (1 - exp(eta)/(2 sqrt(2))) for small exp(eta); use the standard form
    let b = u / (2.0 * 2.0_f64.sqrt());
    let v = u / (1.0 + b);
    if !v.is_finite() || v < 0.0 {
        return Err(HeavyError::BadInput("the degeneracy correction left its domain"));
    }
    Ok(v)
}

/// Built-in potential of an abrupt junction, volts, with narrowing included.
///
/// `V_bi = (kT/q) ln(Na Nd / n_i,eff^2)`, where `n_i,eff` is evaluated at the heavier side
/// because that is where the narrowing acts.
pub fn built_in_potential_v(
    acceptor_cm3: PerCm3,
    donor_cm3: PerCm3,
    ni_cm3: PerCm3,
    temp_k: Kelvin,
) -> Result<f64, HeavyError> {
    let acceptor_cm3 = acceptor_cm3.0;
    let donor_cm3 = donor_cm3.0;
    let ni_cm3 = ni_cm3.0;
    let temp_k = temp_k.0;
    for (name, v) in [("acceptor", acceptor_cm3), ("donor", donor_cm3)] {
        if !v.is_finite() || v <= 0.0 {
            return Err(HeavyError::BadInput(if name == "acceptor" { "acceptor must be finite and > 0" } else { "donor must be finite and > 0" }));
        }
    }
    let heavier = acceptor_cm3.max(donor_cm3);
    let ni_eff = ni_effective_cm3(PerCm3(ni_cm3), PerCm3(heavier), Kelvin(temp_k))?;
    let ratio = acceptor_cm3 * donor_cm3 / (ni_eff * ni_eff);
    if !(ratio > 0.0) {
        return Err(HeavyError::BadInput("the junction has no built-in potential"));
    }
    Ok(KB_EV * temp_k * ratio.ln())
}

#[cfg(test)]
mod tests {
    use super::*;

    const NI: f64 = 1.0e10;
    const T: f64 = 300.0;

    #[test]
    fn narrowing_is_zero_below_the_reference_and_grows_above_it() {
        assert_eq!(band_gap_narrowing_ev(PerCm3(1.0e16)).unwrap(), 0.0, "light doping does not narrow");
        let a = band_gap_narrowing_ev(PerCm3(1.0e18)).unwrap();
        let b = band_gap_narrowing_ev(PerCm3(1.0e20)).unwrap();
        assert!(a > 0.0 && b > a, "narrowing grows with doping: {a} then {b}");
    }

    #[test]
    fn narrowing_raises_the_effective_intrinsic_concentration() {
        let light = ni_effective_cm3(PerCm3(NI), PerCm3(1.0e15), Kelvin(T)).unwrap();
        let heavy = ni_effective_cm3(PerCm3(NI), PerCm3(1.0e20), Kelvin(T)).unwrap();
        assert!((light - NI).abs() / NI < 1e-12, "light doping leaves it alone: {light}");
        assert!(heavy > NI, "heavy doping raises it: {heavy}");
    }

    #[test]
    fn bandgap_narrowing_lowers_the_built_in_potential_at_fixed_doping() {
        // V-10 asks for the flat-band shift narrowing causes. Comparing two different doping
        // levels does NOT show it — a heavier junction has a larger V_bi anyway, and the first
        // version of this test asserted the wrong thing. The effect is a reduction at the same
        // doping, measured against the same junction computed without narrowing.
        let n = 1.0e20;
        let with = built_in_potential_v(atto_units::PerCm3(n), atto_units::PerCm3(n), atto_units::PerCm3(NI), atto_units::Kelvin(T)).unwrap();
        // the same junction with the intrinsic concentration left unnarrowed
        let without = KB_EV * T * ((n * n) / (NI * NI)).ln();
        assert!(with < without, "narrowing pulls V_bi down: {with} vs {without}");
        assert!(with > 0.0, "and it stays positive: {with}");
        let shift = without - with;
        assert!(shift > 0.1, "the shift is large enough to matter in a C-V: {shift} V");
    }

    #[test]
    fn degeneracy_statistics_undershoot_boltzmann() {
        // at a high Fermi level, Boltzmann over-counts
        let eta = 5.0f64;
        let boltzmann = eta.exp();
        let fd = reduced_concentration(eta).unwrap();
        assert!(fd < boltzmann, "degenerate statistics give fewer carriers: {fd} vs {boltzmann}");
        assert!(fd > 0.0, "and a positive number: {fd}");
    }

    #[test]
    fn at_a_low_fermi_level_the_correction_vanishes() {
        let eta = -10.0f64;
        let boltzmann = eta.exp();
        let fd = reduced_concentration(eta).unwrap();
        assert!((fd - boltzmann).abs() / boltzmann < 1e-4, "non-degenerate agrees with Boltzmann to the first-order correction (~1.6e-5 at eta=-10): {fd} vs {boltzmann}");
    }

    #[test]
    fn unusable_input_is_refused() {
        assert!(band_gap_narrowing_ev(PerCm3(0.0)).is_err(), "no doping");
        assert!(band_gap_narrowing_ev(PerCm3(f64::NAN)).is_err(), "not a number");
        assert!(ni_effective_cm3(PerCm3(0.0), PerCm3(1e18), Kelvin(T)).is_err(), "no intrinsic concentration");
        assert!(ni_effective_cm3(PerCm3(NI), PerCm3(1e18), Kelvin(0.0)).is_err(), "zero kelvin");
        assert!(reduced_concentration(f64::INFINITY).is_err(), "infinite Fermi level");
        assert!(built_in_potential_v(atto_units::PerCm3(0.0), atto_units::PerCm3(1e17), atto_units::PerCm3(NI), atto_units::Kelvin(T)).is_err(), "no acceptor");
    }
}
