// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The Massoud thin regime (O-07).
//!
//! Deal–Grove is wrong for the first tens of nanometres of a dry oxide. The measured growth
//! rate there is higher than the linear-parabolic law predicts and falls back onto it as the
//! film thickens, so a 10 nm oxide grown from bare silicon comes out noticeably too thin if
//! only Deal–Grove is used.
//!
//! Massoud's correction adds a decaying term to the growth rate:
//!
//! ```text
//! dx/dt = B / (A + 2x)  +  C exp(-x / L)
//! ```
//!
//! `C` is an acceleration rate and `L` a characteristic length of order the film thickness at
//! which the anomaly disappears. Both are empirical, per orientation, temperature and ambient.

use crate::oxidation::oxidation_nm;
use atto_params::Ambient;
use atto_units::{Kelvin, Nm};

/// The empirical pair behind the thin-regime acceleration.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Massoud {
    /// Acceleration rate `C`, micrometres per hour.
    pub c_um_per_h: f64,
    /// Characteristic length `L`, micrometres.
    pub l_um: f64,
}

impl Massoud {
    /// Defaults for dry oxygen on (100) silicon.
    ///
    /// **These are assumed, not fitted.** O-07 asks for 5/10/20 nm at 800 C within 10% of the
    /// published Massoud fit; no such fit is in this tree and none could be checked in the
    /// sandbox this was written in, so the pair below reproduces the *shape* and the right
    /// order of magnitude of the acceleration, not a published curve. Calibrating it is a
    /// separate block with the paper's table in hand.
    pub fn bundled() -> Self {
        Massoud { c_um_per_h: 0.012, l_um: 0.007 }
    }

    /// Whether the correction is usable at all.
    pub fn validate(&self) -> Result<(), String> {
        if !self.c_um_per_h.is_finite() || self.c_um_per_h < 0.0 {
            return Err("the Massoud rate must be finite and >= 0".into());
        }
        if !self.l_um.is_finite() || self.l_um <= 0.0 {
            return Err("the Massoud length must be finite and > 0".into());
        }
        Ok(())
    }
}

/// Dry-oxide growth with the thin regime included, nm.
///
/// The rate is integrated explicitly in small steps rather than solved in closed form, because
/// the Massoud term makes the equation non-separable. The step is 1% of the characteristic
/// length, which keeps the integration error well below the uncertainty in `C` and `L`.
pub fn oxidation_nm_thin(
    ambient: Ambient,
    temp: Kelvin,
    minutes: f64,
    initial_nm: Nm,
    massoud: &Massoud,
) -> Result<Nm, String> {
    let initial_nm = initial_nm.0;
    massoud.validate()?;
    if !minutes.is_finite() || minutes < 0.0 {
        return Err("growth time must be finite and >= 0".into());
    }
    if !initial_nm.is_finite() || initial_nm < 0.0 {
        return Err("initial oxide must be finite and >= 0".into());
    }
    if minutes == 0.0 {
        return Ok(Nm::new(initial_nm));
    }

    let l_nm = massoud.l_um * 1000.0;
    // Operator split: each step advances the Deal-Grove solution by one step from the current
    // thickness, then adds the Massoud increment at that thickness. The first version divided
    // the full-time growth by a probe interval, which is dimensionally a slope but physically
    // meaningless, and the test caught it.
    let steps = ((minutes * 20.0).clamp(8.0, 20_000.0)) as usize;
    let dt_min = minutes / steps as f64;
    let mut x_nm = initial_nm;
    for _ in 0..steps {
        let dg = oxidation_nm(ambient, temp, dt_min, x_nm).map_err(|e| e.to_string())?;
        let accel_um_per_h = massoud.c_um_per_h * (-(x_nm / l_nm)).exp();
        let extra_nm = (accel_um_per_h * (dt_min / 60.0) * 1000.0).max(0.0);
        x_nm = dg + extra_nm;
    }
    Ok(Nm::new(x_nm))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dry_800() -> (Ambient, Kelvin) {
        (Ambient::Dry, Kelvin::new(1073.15))
    }

    #[test]
    fn the_thin_regime_grows_faster_than_deal_grove() {
        // this is the whole point of O-07: at 10 nm Deal-Grove is short
        let (a, t) = dry_800();
        let m = Massoud::bundled();
        let dg = oxidation_nm(a, t, 20.0, 0.0).expect("deal-grove");
        let thin = oxidation_nm_thin(a, t, 20.0, atto_units::Nm(0.0), &m).expect("thin regime");
        assert!(thin.0 > dg, "the thin regime must exceed Deal-Grove: {} vs {}", thin.0, dg);
    }

    #[test]
    fn the_acceleration_fades_as_the_film_thickens() {
        let (a, t) = dry_800();
        let m = Massoud::bundled();
        let dg_thick = oxidation_nm(a, t, 600.0, 0.0).expect("dg");
        let thin_thick = oxidation_nm_thin(a, t, 600.0, atto_units::Nm(0.0), &m).expect("thin");
        let rel_thin = ((thin_thick.0 - dg_thick) / dg_thick).abs();
        let dg_early = oxidation_nm(a, t, 10.0, 0.0).expect("dg");
        let thin_early = oxidation_nm_thin(a, t, 10.0, atto_units::Nm(0.0), &m).expect("thin");
        let rel_early = ((thin_early.0 - dg_early) / dg_early).abs();
        assert!(rel_thin < rel_early, "the correction must shrink with thickness: {} vs {}", rel_thin, rel_early);
    }

    #[test]
    fn growth_is_monotone_in_time() {
        let (a, t) = dry_800();
        let m = Massoud::bundled();
        let mut last = 0.0;
        for minutes in [1.0, 5.0, 20.0, 60.0, 240.0] {
            let x = oxidation_nm_thin(a, t, minutes, atto_units::Nm(0.0), &m).expect("growth").0;
            assert!(x > last, "oxide cannot shrink: {x} after {last} at {minutes} min");
            last = x;
        }
    }

    #[test]
    fn zero_time_changes_nothing() {
        let (a, t) = dry_800();
        let x = oxidation_nm_thin(a, t, 0.0, atto_units::Nm(12.0), &Massoud::bundled()).expect("no growth");
        assert_eq!(x.0, 12.0, "no time, no oxide");
    }

    #[test]
    fn a_correction_that_cannot_exist_is_refused() {
        let (a, t) = dry_800();
        assert!(Massoud { c_um_per_h: -1.0, l_um: 0.007 }.validate().is_err(), "negative rate");
        assert!(Massoud { c_um_per_h: 0.012, l_um: 0.0 }.validate().is_err(), "zero length");
        assert!(oxidation_nm_thin(a, t, -1.0, atto_units::Nm(0.0), &Massoud::bundled()).is_err(), "negative time");
        assert!(oxidation_nm_thin(a, t, 10.0, atto_units::Nm(f64::NAN), &Massoud::bundled()).is_err(), "not a number");
    }
}
