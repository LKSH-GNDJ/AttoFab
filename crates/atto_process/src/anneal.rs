// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Anneal: numerical 1-D diffusion and thermal budget (build PRD D-10, D-11).
//!
//! the earlier build broadened the implant profile analytically with `sigma^2 += 2 D t`, which
//! cannot express a concentration-dependent diffusivity, a finite wafer, or a
//! multi-step thermal history. This solves the diffusion equation instead, on
//! the column profile the field already holds, with an implicit tridiagonal
//! scheme (unconditionally stable, conservative under the no-flux boundaries).

use atto_units::{Kelvin, Nm};
use atto_dopant::{DopantField, DopantRecord, Type};
use atto_geom::column::ColumnGrid;
use atto_geom::svo::Material;
use atto_params::{Species, KB_EV};

#[derive(Debug, Clone, PartialEq)]
/// Why an anneal step was refused.
pub enum AnnealError {
    NonFiniteInput,
    InvalidTemperature,
    Empty,
}

impl std::fmt::Display for AnnealError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            AnnealError::NonFiniteInput => write!(f, "anneal inputs must be finite and non-negative"),
            AnnealError::InvalidTemperature => write!(f, "temperature must be above absolute zero"),
            AnnealError::Empty => write!(f, "the column has no dopant to diffuse"),
        }
    }
}
impl std::error::Error for AnnealError {}

/// One thermal step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ThermalStep {
    pub temp_c: f64,
    pub minutes: f64,
}

/// The thermal history of a wafer. Every anneal appends; nothing is overwritten.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ThermalHistory {
    pub steps: Vec<ThermalStep>,
}

impl ThermalHistory {
    /// Check the history's invariants, returning every violation (PRD A-11).
    ///
    /// A thermal budget that goes backwards is the failure this catches: a step
    /// appended out of order makes every accumulated total wrong from there on,
    /// and nothing downstream would notice.
    pub fn check(&self) -> Vec<String> {
        let mut out = Vec::new();
        let mut prev: Option<f64> = None;
        for (i, step) in self.steps.iter().enumerate() {
            if !step.temp_c.is_finite() || step.temp_c <= -273.15 {
                out.push(format!("step {i} has a temperature of {} C", step.temp_c));
            }
            if !step.minutes.is_finite() || step.minutes < 0.0 {
                out.push(format!("step {i} has a duration of {} min", step.minutes));
            }
            if let Some(t) = prev {
                if step.temp_c < t - 1e-9 {
                    out.push(format!("step {i} is colder than the step before it"));
                }
            }
            prev = Some(step.temp_c);
        }
        out
    }

    /// An empty history.
    pub fn new() -> Self { ThermalHistory { steps: Vec::new() } }

    /// Append a step, refusing one that goes backwards in time.
    pub fn push(&mut self, step: ThermalStep) -> Result<(), AnnealError> {
        if !(step.temp_c.is_finite() && step.minutes.is_finite()) || step.minutes < 0.0 {
            return Err(AnnealError::NonFiniteInput);
        }
        if step.temp_c + 273.15 <= 0.0 {
            return Err(AnnealError::InvalidTemperature);
        }
        self.steps.push(step);
        Ok(())
    }

    /// The total time at temperature.
    pub fn total_minutes(&self) -> f64 { self.steps.iter().map(|s| s.minutes).sum() }

    /// The hottest step.
    pub fn peak_c(&self) -> f64 {
        self.steps.iter().map(|s| s.temp_c).fold(f64::NEG_INFINITY, f64::max)
    }

    /// Intrinsic diffusivity at a temperature, cm^2/s.
    pub fn diffusivity_at(species: Species, temp_c: f64) -> f64 {
        let (d0, ea) = species.diffusivity();
        let t_k = temp_c + 273.15;
        if t_k <= 0.0 {
            return 0.0;
        }
        d0 * (-ea / (KB_EV * t_k)).exp()
    }

    /// Accumulated thermal budget in cm^2: sum of D(T_i) * t_i.
    ///
    /// This is the quantity that makes two 5-minute anneals equivalent to one
    /// 10-minute anneal for a first-order process, and it is what the earlier build had no
    /// representation for at all.
    pub fn accumulated_dt_cm2(&self, species: Species) -> f64 {
        self.steps
            .iter()
            .map(|s| Self::diffusivity_at(species, s.temp_c) * s.minutes * 60.0)
            .sum()
    }
}

/// Anneal configuration for one step.
#[derive(Debug, Clone, Copy)]
pub struct AnnealConfig {
    pub species: Species,
    pub kind: Type,
    pub temp_c: f64,
    pub minutes: f64,
    /// Time steps in the implicit solve.
    pub steps: usize,
    /// Concentration dependence: `D = D_i * (1 + alpha * C / ni)`. `alpha` is a
    /// declared assumption, not a measured parameter (see the manifest).
    pub alpha: f64,
    /// Intrinsic carrier concentration at the anneal temperature, cm^-3.
    pub ni_cm3: f64,
}

impl AnnealConfig {
    /// One thermal step.
    pub fn new(species: Species, kind: Type, temp: Kelvin, minutes: f64) -> Self {
        let temp_c = temp.celsius();
        AnnealConfig { species, kind, temp_c, minutes, steps: 24, alpha: 1.0, ni_cm3: 1.0e19 }
    }

    fn validate(&self) -> Result<(), AnnealError> {
        if ![self.temp_c, self.minutes, self.alpha, self.ni_cm3].iter().all(|v| v.is_finite()) {
            return Err(AnnealError::NonFiniteInput);
        }
        if self.minutes < 0.0 || self.alpha < 0.0 || self.ni_cm3 <= 0.0 {
            return Err(AnnealError::NonFiniteInput);
        }
        if self.temp_c + 273.15 <= 0.0 {
            return Err(AnnealError::InvalidTemperature);
        }
        if self.steps == 0 {
            return Err(AnnealError::NonFiniteInput);
        }
        Ok(())
    }
}

/// Diffuse one profile in place.
///
/// Backward Euler on a cell-centred grid with no-flux ends:
/// `-a[i-1] c[i-1] + (1 + a[i-1] + a[i]) c[i] - a[i] c[i+1] = c[i]^n`
/// where `a[i] = D[i+1/2] dt / h^2`. Solved by the Thomas algorithm.
pub fn diffuse(c: &mut [f64], spacing: Nm, cfg: &AnnealConfig) -> Result<(), AnnealError> {
    // A-03, and the argument this project has the most reason to get right: a
    // depth in nanometres read as centimetres is exactly the bug that made the
    // dose-loss integration return zero everywhere.
    let h_cm = spacing.value() * 1e-7;
    cfg.validate()?;
    let n = c.len();
    if n == 0 {
        return Err(AnnealError::Empty);
    }
    if n == 1 {
        return Ok(());
    }
    let total_s = cfg.minutes * 60.0;
    if total_s <= 0.0 {
        return Ok(());
    }
    let dt = total_s / cfg.steps as f64;
    let d_int = ThermalHistory::diffusivity_at(cfg.species, cfg.temp_c);
    if !(d_int.is_finite() && d_int >= 0.0) {
        return Err(AnnealError::InvalidTemperature);
    }
    let mut next = c.to_vec();
    let mut a = vec![0.0f64; n - 1];
    let mut cp = vec![0.0f64; n];
    let mut dp = vec![0.0f64; n];
    for _ in 0..cfg.steps {
        // face diffusivities from the current profile (Picard lagging)
        for i in 0..n - 1 {
            let d_i = d_int * (1.0 + cfg.alpha * c[i] / cfg.ni_cm3);
            let d_j = d_int * (1.0 + cfg.alpha * c[i + 1] / cfg.ni_cm3);
            a[i] = 0.5 * (d_i + d_j) * dt / (h_cm * h_cm);
        }
        // Thomas algorithm for lower[i] = -a[i-1], diag[i] = 1 + a[i-1] + a[i],
        // upper[i] = -a[i]. cp carries upper/m, which is *negative* here; the
        // back substitution subtracts it, so the two signs must not be flipped
        // together (that cancels and silently returns the wrong field).
        let mut m = 1.0 + a[0];
        cp[0] = -a[0] / m;
        dp[0] = c[0] / m;
        for i in 1..n - 1 {
            m = 1.0 + a[i - 1] + a[i] + a[i - 1] * cp[i - 1];
            cp[i] = -a[i] / m;
            dp[i] = (c[i] + a[i - 1] * dp[i - 1]) / m;
        }
        m = 1.0 + a[n - 2] + a[n - 2] * cp[n - 2];
        dp[n - 1] = (c[n - 1] + a[n - 2] * dp[n - 2]) / m;
        // back substitution
        next[n - 1] = dp[n - 1];
        for i in (0..n - 1).rev() {
            next[i] = dp[i] - cp[i] * next[i + 1];
        }
        c.copy_from_slice(&next);
    }
    Ok(())
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// What an anneal moved and how far.
pub struct AnnealStats {
    pub columns: usize,
    pub leaves: usize,
    /// Accumulated Dt for this step, cm^2.
    pub dt_cm2: f64,
    /// Total dose before and after, for a conservation check.
    pub dose_before_cm2: f64,
    pub dose_after_cm2: f64,
}

/// Anneal the dopant field: diffuse every column's profile for one species.
pub fn anneal_field(field: &mut DopantField, grid: &ColumnGrid, cfg: &AnnealConfig) -> Result<AnnealStats, AnnealError> {
    cfg.validate()?;
    let leaf = field.leaf_nm();
    let h_cm = leaf * 1.0e-7;
    let mut st = AnnealStats { dt_cm2: ThermalHistory::diffusivity_at(cfg.species, cfg.temp_c) * cfg.minutes * 60.0, ..Default::default() };

    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            if !grid.active[k] {
                continue;
            }
            let col = &grid.cols[k];
            let Some(si) = col.iter().rev().find(|r| r.material == Material::Silicon) else {
                continue;
            };
            let x = grid.x0 + (ix as f64 + 0.5) * grid.pitch;
            let y = grid.y0 + (iy as f64 + 0.5) * grid.pitch;
            // gather the column's leaves, top down
            let mut keys = Vec::new();
            let mut d = leaf * 0.5;
            while si.z1 - d > si.z0 {
                let key = field.key_for(x, y, si.z1 - d);
                let has = field.records(key).iter().any(|r| r.species == cfg.species && r.kind == cfg.kind);
                if has {
                    keys.push(key);
                } else if !keys.is_empty() {
                    break; // the doped region is contiguous from the top
                }
                d += leaf;
            }
            if keys.len() < 2 {
                continue;
            }
            let mut c: Vec<f64> = keys
                .iter()
                .map(|k| match field.records(*k).iter().find(|r| r.species == cfg.species && r.kind == cfg.kind) {
                    // no record means the leaf carries none of this species,
                    // which is a real zero rather than a missing value
                    Some(r) => r.concentration_cm3,
                    None => 0.0,
                })
                .collect();
            let before: f64 = c.iter().sum::<f64>() * h_cm;
            diffuse(&mut c, Nm::new(h_cm * 1e7), cfg)?;
            let after: f64 = c.iter().sum::<f64>() * h_cm;
            for (key, v) in keys.iter().zip(c.iter()) {
                field.set(*key, DopantRecord::new(cfg.species, cfg.kind, *v));
            }
            st.columns += 1;
            st.leaves += keys.len();
            st.dose_before_cm2 += before;
            st.dose_after_cm2 += after;
        }
    }
    Ok(st)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A Gaussian centred in the middle of the domain. Centring matters: a
    /// profile sitting on the no-flux boundary is a half-Gaussian, and its
    /// second moment is not sigma^2.
    fn gaussian(n: usize, h_nm: f64, sigma_nm: f64, peak: f64) -> Vec<f64> {
        let mid = n as f64 * 0.5;
        (0..n)
            .map(|i| {
                let z = (i as f64 + 0.5 - mid) * h_nm;
                peak * (-(z * z) / (2.0 * sigma_nm * sigma_nm)).exp()
            })
            .collect()
    }

    fn second_moment_sigma(c: &[f64], h_nm: f64) -> f64 {
        let m0: f64 = c.iter().sum();
        let m1: f64 = c.iter().enumerate().map(|(i, v)| v * (i as f64 + 0.5) * h_nm).sum::<f64>() / m0;
        let m2: f64 = c.iter().enumerate().map(|(i, v)| v * ((i as f64 + 0.5) * h_nm - m1).powi(2)).sum::<f64>() / m0;
        m2.sqrt()
    }

    #[test]
    fn dilute_limit_matches_the_analytic_gaussian() {
        // sigma^2 = sigma0^2 + 2 D t. Keep the peak low so D stays intrinsic.
        let h = 4.0;
        let n = 200;
        let sigma0 = 40.0;
        let mut c = gaussian(n, h, sigma0, 1.0e15);
        let cfg = AnnealConfig { alpha: 0.0, ..AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 10.0) };
        diffuse(&mut c, Nm::new(h), &cfg).unwrap();
        let d = ThermalHistory::diffusivity_at(Species::Boron, 1000.0);
        let expect = (sigma0 * sigma0 + 2.0 * d * 600.0 * 1e14).sqrt();
        let got = second_moment_sigma(&c, h);
        assert!((got - expect).abs() / expect < 0.05, "sigma {got:.1} vs analytic {expect:.1} nm");
    }

    #[test]
    fn mass_is_conserved_under_the_no_flux_boundaries() {
        let mut c = gaussian(120, 4.0, 30.0, 1.0e18);
        let before: f64 = c.iter().sum();
        let cfg = AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 30.0);
        diffuse(&mut c, Nm::new(4e-7 * 1e7), &cfg).unwrap();
        let after: f64 = c.iter().sum();
        assert!((after - before).abs() / before < 1e-9, "mass {before} -> {after}");
    }

    #[test]
    fn a_longer_or_hotter_anneal_broadens_more() {
        let base = {
            let mut c = gaussian(200, 4.0, 30.0, 1.0e16);
            diffuse(&mut c, Nm::new(4e-7 * 1e7), &AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 10.0)).unwrap();
            second_moment_sigma(&c, 4.0)
        };
        let longer = {
            let mut c = gaussian(200, 4.0, 30.0, 1.0e16);
            diffuse(&mut c, Nm::new(4e-7 * 1e7), &AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 60.0)).unwrap();
            second_moment_sigma(&c, 4.0)
        };
        let hotter = {
            let mut c = gaussian(200, 4.0, 30.0, 1.0e16);
            diffuse(&mut c, Nm::new(4e-7 * 1e7), &AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1100.0), 10.0)).unwrap();
            second_moment_sigma(&c, 4.0)
        };
        assert!(longer > base, "{longer} > {base}");
        assert!(hotter > base, "{hotter} > {base}");
    }

    #[test]
    fn two_short_anneals_match_one_long_one() {
        let cfg5 = AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 5.0);
        let mut split = gaussian(200, 4.0, 30.0, 1.0e16);
        diffuse(&mut split, Nm::new(4e-7 * 1e7), &cfg5).unwrap();
        diffuse(&mut split, Nm::new(4e-7 * 1e7), &cfg5).unwrap();
        let mut single = gaussian(200, 4.0, 30.0, 1.0e16);
        diffuse(&mut single, Nm::new(4e-7 * 1e7), &AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 10.0)).unwrap();
        let a = second_moment_sigma(&split, 4.0);
        let b = second_moment_sigma(&single, 4.0);
        assert!((a - b).abs() / b < 0.02, "5+5 min ({a:.2}) vs 10 min ({b:.2})");
    }

    #[test]
    fn concentration_dependence_broadens_more_than_the_dilute_case() {
        let cfg_conc = AnnealConfig { alpha: 2.0, ni_cm3: 1.0e18, ..AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 20.0) };
        let cfg_dilute = AnnealConfig { alpha: 0.0, ..cfg_conc };
        let mut a = gaussian(200, 4.0, 30.0, 1.0e19);
        let mut b = a.clone();
        diffuse(&mut a, Nm::new(4e-7 * 1e7), &cfg_conc).unwrap();
        diffuse(&mut b, Nm::new(4e-7 * 1e7), &cfg_dilute).unwrap();
        assert!(second_moment_sigma(&a, 4.0) > second_moment_sigma(&b, 4.0));
    }

    #[test]
    fn thermal_budget_is_additive() {
        let mut one = ThermalHistory::new();
        one.push(ThermalStep { temp_c: 1000.0, minutes: 10.0 }).unwrap();
        let mut two = ThermalHistory::new();
        two.push(ThermalStep { temp_c: 1000.0, minutes: 5.0 }).unwrap();
        two.push(ThermalStep { temp_c: 1000.0, minutes: 5.0 }).unwrap();
        let a = one.accumulated_dt_cm2(Species::Boron);
        let b = two.accumulated_dt_cm2(Species::Boron);
        assert!((a - b).abs() / a < 1e-12, "{a} vs {b}");
        assert!((one.total_minutes() - 10.0).abs() < 1e-12);
        assert_eq!(one.peak_c(), 1000.0);
    }

    #[test]
    fn diffusivity_matches_the_textbook_value() {
        // B in Si at 1000 C is about 1.5e-14 cm^2/s.
        let d = ThermalHistory::diffusivity_at(Species::Boron, 1000.0);
        assert!((d - 1.5e-14).abs() / 1.5e-14 < 0.1, "D = {d:.2e}");
    }

    #[test]
    fn a_healthy_history_passes_its_own_check() {
        let mut h = ThermalHistory::new();
        h.push(ThermalStep { temp_c: 900.0, minutes: 30.0 }).unwrap();
        h.push(ThermalStep { temp_c: 1000.0, minutes: 10.0 }).unwrap();
        assert!(h.check().is_empty(), "{:?}", h.check());
    }

    #[test]
    fn a_history_that_goes_backwards_is_reported() {
        let mut h = ThermalHistory::new();
        h.push(ThermalStep { temp_c: 1000.0, minutes: 10.0 }).unwrap();
        h.steps.push(ThermalStep { temp_c: 800.0, minutes: 10.0 });
        let v = h.check();
        assert!(!v.is_empty(), "a colder step after a hotter one must be reported");
        assert!(v[0].contains("colder"), "{v:?}");
    }

    #[test]
    fn bad_inputs_are_rejected() {
        let mut c = vec![1.0, 2.0, 3.0];
        assert!(diffuse(&mut c, Nm::new(4e-7 * 1e7), &AnnealConfig { temp_c: f64::NAN, ..AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(1000.0), 10.0) }).is_err());
        assert!(diffuse(&mut c, Nm::new(4e-7 * 1e7), &AnnealConfig::new(Species::Boron, Type::Acceptor, Kelvin::from_celsius(-400.0), 10.0)).is_err());
        assert!(ThermalHistory::new().push(ThermalStep { temp_c: 1000.0, minutes: -1.0 }).is_err());
    }
}
