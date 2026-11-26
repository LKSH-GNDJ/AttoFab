// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Resist development: how fast the developer eats the resist, as a function of dose.
//!
//! The simulation used to decide "removed or not" with a single threshold on intensity.
//! That gives vertical sidewalls by construction: a column is either gone or untouched, so
//! the resist CD is exactly the mask CD and the profile has no slope at all.
//!
//! A real developer removes resist at a rate that depends on how much photoactive compound
//! is left. Where the dose was high the rate is fast; where it was low the rate is slow but
//! not zero. Given a development time, that produces a *depth* per column instead of a flag,
//! and depth varying across the edge of a feature is what a sidewall angle is.
//!
//! The rate law is Mack's: `R(M) = R_max (a+1)(1-M)^n / (a + (1-M)^n)`, with the PAC
//! fraction `M` falling exponentially with dose.

/// Mack development parameters for one resist.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DevelopModel {
    /// Rate where the resist is fully exposed, nm/min.
    pub r_max_nm_per_min: f64,
    /// Mack `a`: the rate at zero exposure relative to `R_max`, as a fraction.
    pub a: f64,
    /// Mack `n`: selectivity. Higher is a sharper response to dose.
    pub n: f64,
    /// Dill `C`, per unit dose ratio: how fast the PAC is consumed.
    pub dill_c: f64,
}

impl Default for DevelopModel {
    fn default() -> Self {
        // typical positive i-line / KrF numbers; replace with measured values per resist
        DevelopModel { r_max_nm_per_min: 300.0, a: 0.02, n: 8.0, dill_c: 2.5 }
    }
}

impl DevelopModel {
    /// A model is usable only if every parameter is finite and the rate can be positive.
    pub fn validate(&self) -> Result<(), String> {
        if !self.r_max_nm_per_min.is_finite() || self.r_max_nm_per_min <= 0.0 {
            return Err("R_max must be finite and > 0".into());
        }
        if !self.a.is_finite() || self.a < 0.0 {
            return Err("a must be finite and >= 0".into());
        }
        if !self.n.is_finite() || self.n <= 0.0 {
            return Err("n must be finite and > 0".into());
        }
        if !self.dill_c.is_finite() || self.dill_c < 0.0 {
            return Err("Dill C must be finite and >= 0".into());
        }
        Ok(())
    }

    /// Fraction of photoactive compound still present after this dose ratio.
    ///
    /// 1.0 is unexposed, 0.0 is fully bleached. `dose_ratio` is exposure over the dose to
    /// clear, so 1.0 is the nominal clearing dose.
    pub fn pac_fraction(&self, dose_ratio: f64) -> f64 {
        if !dose_ratio.is_finite() || dose_ratio <= 0.0 {
            return 1.0;
        }
        (-self.dill_c * dose_ratio).exp().clamp(0.0, 1.0)
    }

    /// Development rate at a given PAC fraction, nm/min.
    pub fn rate_nm_per_min(&self, pac: f64) -> f64 {
        let m = pac.clamp(0.0, 1.0);
        let exposed = 1.0 - m;
        let num = (self.a + 1.0) * exposed.powf(self.n);
        let den = self.a + exposed.powf(self.n);
        if den <= 0.0 {
            return 0.0;
        }
        self.r_max_nm_per_min * num / den
    }

    /// Rate at a given dose ratio, nm/min.
    pub fn rate_at_dose(&self, dose_ratio: f64) -> f64 {
        self.rate_nm_per_min(self.pac_fraction(dose_ratio))
    }

    /// Depth developed into the resist after `time_min`, nm.
    ///
    /// Clipped to the film: the developer cannot remove resist that is not there, and a
    /// column that clears early does not keep eating the film below it.
    pub fn depth_nm(&self, dose_ratio: f64, time_min: atto_units::Minutes, thickness: atto_units::Nm) -> atto_units::Nm {
        let (t, th) = (time_min.0, thickness.0);
        if !t.is_finite() || t <= 0.0 || !th.is_finite() || th <= 0.0 {
            return atto_units::Nm::new(0.0);
        }
        atto_units::Nm::new((self.rate_at_dose(dose_ratio) * t).clamp(0.0, th))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m() -> DevelopModel { DevelopModel::default() }

    #[test]
    fn a_bigger_dose_develops_faster() {
        let d = m();
        let low = d.rate_at_dose(0.5);
        let high = d.rate_at_dose(2.0);
        assert!(high > low, "more dose must develop faster: low={low} high={high}");
    }

    #[test]
    fn an_unexposed_column_barely_develops() {
        // the dark rate is small but not zero: that is why an unexposed area loses a little
        // resist and why resist thickness is never exactly the coat thickness
        let d = m();
        let r = d.rate_at_dose(1e-6);
        assert!(r < d.r_max_nm_per_min * 0.05, "dark rate must be small: {r}");
        assert!(r > 0.0, "but not zero");
    }

    #[test]
    fn a_fully_exposed_column_develops_at_the_maximum_rate() {
        let d = m();
        let r = d.rate_at_dose(50.0);
        assert!((r - d.r_max_nm_per_min).abs() / d.r_max_nm_per_min < 0.01, "near R_max: {r}");
    }

    #[test]
    fn the_developed_depth_is_clipped_to_the_film() {
        let d = m();
        assert!((d.depth_nm(10.0, atto_units::Minutes::new(100.0), atto_units::Nm::new(150.0)).0 - 150.0).abs() < 1e-9, "cannot remove more than the film");
        let part = d.depth_nm(1e-6, atto_units::Minutes::new(0.01), atto_units::Nm::new(150.0)).0;
        assert!(part >= 0.0 && part < 150.0, "a short develop is partial: {part}");
    }

    #[test]
    fn a_sharper_resist_has_a_bigger_ratio_between_exposed_and_dark() {
        let mut soft = m(); soft.n = 2.0;
        let mut sharp = m(); sharp.n = 12.0;
        let ratio = |d: &DevelopModel| d.rate_at_dose(1.0) / d.rate_at_dose(1e-3).max(1e-12);
        assert!(ratio(&sharp) > ratio(&soft), "higher n must be more selective");
    }

    #[test]
    fn a_model_that_cannot_be_used_is_refused() {
        let mut d = m(); d.r_max_nm_per_min = 0.0;
        assert!(d.validate().is_err(), "zero rate");
        let mut d = m(); d.n = f64::NAN;
        assert!(d.validate().is_err(), "not a number");
        let mut d = m(); d.a = -1.0;
        assert!(d.validate().is_err(), "negative a");
        assert!(m().validate().is_ok(), "the default is usable");
    }
}
