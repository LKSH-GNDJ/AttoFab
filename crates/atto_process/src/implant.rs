// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Implant profiles from tabulated moments (build PRD D-05, D-08).
//!
//! The profile is a Pearson IV distribution rebuilt from (Rp, dRp, skewness,
//! kurtosis) and normalised numerically. Where the moments do not admit a
//! Pearson IV the code falls back to a Gaussian and says so, rather than
//! producing a distribution with poles.

use atto_units::{Nm, PerCm2};
use atto_params::{moments_at, Species};

#[derive(Debug, Clone, PartialEq)]
/// Why an implant profile could not be built.
pub enum ProfileError {
    NonFiniteInput,
    Empty,
}

impl std::fmt::Display for ProfileError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProfileError::NonFiniteInput => write!(f, "implant inputs must be finite"),
            ProfileError::Empty => write!(f, "profile has no support"),
        }
    }
}
impl std::error::Error for ProfileError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The analytic shape a profile takes.
pub enum Shape {
    /// Moments admitted a Pearson IV.
    PearsonIv,
    /// Moments did not admit one; a Gaussian of the same Rp and dRp was used.
    GaussianFallback,
}

/// A normalised implant profile in the standardised variable z = (x - Rp)/dRp.
#[derive(Debug, Clone)]
pub struct Profile {
    pub rp_nm: f64,
    pub drp_nm: f64,
    pub skew: f64,
    pub kurt: f64,
    pub shape: Shape,
    /// (z, f(z)) on a uniform grid, normalised so that the trapezoid integral is 1.
    grid: Vec<(f64, f64)>,
    /// Fraction of the dose that channels along <100> (build PRD D-09).
    pub channel_fraction: f64,
    /// Decay length of the channelling tail, nm.
    pub channel_decay_nm: f64,
}

impl Profile {
    /// Build from species and energy, using the tabulated moments.
    pub fn from_species(s: Species, energy_kev: f64) -> Result<Self, ProfileError> {
        let m = moments_at(s, energy_kev);
        Self::from_moments(Nm::new(m.rp_nm), Nm::new(m.drp_nm), m.skew, m.kurt)
    }

    /// A profile from four moments.
    pub fn from_moments(rp: Nm, drp: Nm, skew: f64, kurt: f64) -> Result<Self, ProfileError> {
    let (rp_nm, drp_nm) = (rp.value(), drp.value());
        if ![rp_nm, drp_nm, skew, kurt].iter().all(|v| v.is_finite()) {
            return Err(ProfileError::NonFiniteInput);
        }
        if !(rp_nm > 0.0 && drp_nm > 0.0) {
            return Err(ProfileError::NonFiniteInput);
        }
        let (grid, shape) = match pearson_iv_grid(skew, kurt) {
            Some(g) => (standardise(g), Shape::PearsonIv),
            None => (gaussian_grid(), Shape::GaussianFallback),
        };
        Ok(Profile { rp_nm, drp_nm, skew, kurt, shape, grid, channel_fraction: 0.0, channel_decay_nm: 0.0 })
    }

    /// The density f(z), linearly interpolated and zero outside the grid.
    pub fn density(&self, z: f64) -> f64 {
        let n = self.grid.len();
        if n < 2 {
            return 0.0;
        }
        let (z0, _) = self.grid[0];
        let (zn, _) = self.grid[n - 1];
        if z <= z0 || z >= zn {
            return 0.0;
        }
        let dz = (zn - z0) / (n - 1) as f64;
        let t = (z - z0) / dz;
        let i = t.floor() as usize;
        let f = t - i as f64;
        let a = self.grid[i].1;
        let b = self.grid[(i + 1).min(n - 1)].1;
        a + (b - a) * f
    }

    /// Concentration (cm^-3) at a depth below the surface.
    pub fn concentration_cm3(&self, depth: Nm, dose_cm2: f64) -> f64 {
        let depth_nm = depth.value();
        if depth_nm < 0.0 {
            return 0.0;
        }
        let z = (depth_nm - self.rp_nm) / self.drp_nm;
        // N(x) = dose * f(z) / dRp, with dRp converted to cm.
        let amorphous = dose_cm2 * self.density(z) / (self.drp_nm * 1.0e-7);
        if self.channel_fraction <= 0.0 || self.channel_decay_nm <= 0.0 {
            return amorphous;
        }
        // The channelled share of the dose runs ahead of the amorphous peak and
        // decays exponentially along <100>. It is taken from the amorphous part
        // so the total dose is unchanged.
        let x0 = self.rp_nm;
        if depth_nm <= x0 {
            return amorphous * (1.0 - self.channel_fraction);
        }
        let tail = dose_cm2 * self.channel_fraction
            * (-(depth_nm - x0) / self.channel_decay_nm).exp()
            / (self.channel_decay_nm * 1.0e-7);
        amorphous * (1.0 - self.channel_fraction) + tail
    }

    /// Fraction of the dose that lands beyond `depth_nm`, by trapezoid over the
    /// stored profile. The profile is normalised, so this is a retained dose.
    pub fn fraction_beyond(&self, depth_nm: f64) -> f64 {
        if self.grid.len() < 2 || !depth_nm.is_finite() {
            return 0.0;
        }
        // The stored grid is in normalised z = (x - rp)/drp, not nanometres.
        // Comparing a nanometre cutoff against it directly returns zero for
        // every depth, which is how this was wrong the first time.
        let z_cut = (depth_nm - self.rp_nm) / self.drp_nm;
        let dz = self.grid[1].0 - self.grid[0].0;
        let mut acc = 0.0;
        for (z, f) in &self.grid {
            if *z + dz * 0.5 > z_cut {
                acc += f * dz;
            }
        }
        acc.clamp(0.0, 1.0)
    }

    /// Fraction of the dose that reaches the silicon through a screen oxide
    /// (build PRD D-12). Dopant stopping in the oxide is lost to the device.
    pub fn retained_in_silicon(&self, oxide_nm: f64) -> f64 {
        if !(oxide_nm > 0.0) || !oxide_nm.is_finite() {
            return 1.0;
        }
        self.fraction_beyond(oxide_nm)
    }

    /// Apply channelling for a tilt angle, in degrees.
    ///
    /// An implant aligned with the <100> axis channels strongly; tilting a few
    /// degrees scatters the ions off the open channels and the tail collapses.
    /// Both constants are **assumed** and registered in the manifest.
    pub fn with_tilt(mut self, tilt_deg: f64) -> Self {
        let t = tilt_deg.max(0.0);
        // Channelling needs alignment to within about a degree of the <100> axis,
        // so the roll-off is sharp: a 7 degree tilt, the standard off-axis cut,
        // suppresses it almost entirely. A gentle cosine roll-off leaves 94% of
        // the tail at 7 degrees, which is not what a tilted implant does.
        let roll = (-(t / 1.2) * (t / 1.2)).exp();
        self.channel_fraction = 0.25 * roll;
        self.channel_decay_nm = 60.0;
        self
    }

    /// Depth at which the concentration falls to `background_cm3`, searching the
    /// deep tail. Returns None when the peak never exceeds the background.
    pub fn junction_depth_nm(&self, dose: PerCm2, background_cm3: f64) -> Option<f64> {
        let dose_cm2 = dose.value();
        let peak = self.concentration_cm3(Nm::new(self.rp_nm), dose_cm2);
        if !(peak > background_cm3) {
            return None;
        }
        let mut lo = self.rp_nm;
        let mut hi = self.rp_nm + 40.0 * self.drp_nm;
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if self.concentration_cm3(Nm::new(mid), dose_cm2) > background_cm3 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some(lo)
    }

    /// Fraction of the profile that lies above the surface and is therefore
    /// reflected rather than implanted.
    pub fn reflected_fraction(&self) -> f64 {
        let cut = -self.rp_nm / self.drp_nm;
        let n = self.grid.len();
        if n < 2 {
            return 0.0;
        }
        let dz = self.grid[1].0 - self.grid[0].0;
        let mut s = 0.0;
        for (z, f) in self.grid.iter() {
            if *z < cut {
                s += f * dz;
            }
        }
        s.clamp(0.0, 1.0)
    }

    /// Numerically integrate the concentration to check dose conservation.
    pub fn integrated_dose_cm2(&self, dose_cm2: f64) -> f64 {
        let zmax = self.rp_nm + 14.0 * self.drp_nm;
        let n = 4000;
        let dx = zmax / n as f64;
        let mut s = 0.0;
        for i in 0..n {
            let x = (i as f64 + 0.5) * dx;
            s += self.concentration_cm3(Nm::new(x), dose_cm2) * dx;
        }
        s * 1.0e-7
    }

    /// The first two moments of the reconstructed density, for self-checking.
    pub fn computed_moments(&self) -> (f64, f64) {
        let n = self.grid.len();
        if n < 2 {
            return (0.0, 3.0);
        }
        let (m0, m1, m2, m3, m4) = self.raw_moments();
        let _ = n;
        let var = m2 - m1 * m1;
        if var <= 0.0 {
            return (0.0, 3.0);
        }
        let sd = var.sqrt();
        let skew = (m3 - 3.0 * m1 * m2 + 2.0 * m1 * m1 * m1) / (sd * sd * sd);
        let kurt = (m4 - 4.0 * m1 * m3 + 6.0 * m1 * m1 * m2 - 3.0 * m1.powi(4)) / (var * var);
        let _ = m0;
        (skew, kurt)
    }

    fn raw_moments(&self) -> (f64, f64, f64, f64, f64) {
        let n = self.grid.len();
        let (z0, _) = self.grid[0];
        let (zn, _) = self.grid[n - 1];
        let dz = (zn - z0) / (n - 1) as f64;
        let (mut m0, mut m1, mut m2, mut m3, mut m4) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for (z, f) in self.grid.iter() {
            let w = f * dz;
            m0 += w;
            m1 += w * z;
            m2 += w * z * z;
            m3 += w * z * z * z;
            m4 += w * z * z * z * z;
        }
        (m0, m1, m2, m3, m4)
    }
}

/// Shift and scale so the density has zero mean and unit variance. The Pearson
/// parameters give the right *shape*; this puts it on the standardised axis.
fn standardise(mut v: Vec<(f64, f64)>) -> Vec<(f64, f64)> {
    let n = v.len();
    if n < 2 {
        return v;
    }
    let dz = v[1].0 - v[0].0;
    let (mut m0, mut m1, mut m2) = (0.0, 0.0, 0.0);
    for (z, f) in v.iter() {
        let w = f * dz;
        m0 += w;
        m1 += w * z;
        m2 += w * z * z;
    }
    if !(m0 > 0.0) {
        return v;
    }
    let mean = m1 / m0;
    let var = m2 / m0 - mean * mean;
    if !(var > 1e-12) {
        return v;
    }
    let sd = var.sqrt();
    let lo = (v[0].0 - mean) / sd;
    let hi = (v[n - 1].0 - mean) / sd;
    let dz_new = (hi - lo) / (n - 1) as f64;
    for (i, p) in v.iter_mut().enumerate() {
        p.0 = lo + i as f64 * dz_new;
    }
    let mut sum = 0.0;
    for (_, f) in v.iter() {
        sum += f * dz_new;
    }
    for p in v.iter_mut() {
        p.1 /= sum;
    }
    v
}

const Z_LO: f64 = -14.0;
const Z_HI: f64 = 14.0;
const NZ: usize = 5601;

fn gaussian_grid() -> Vec<(f64, f64)> {
    let mut v = Vec::with_capacity(NZ);
    let mut sum = 0.0;
    let dz = (Z_HI - Z_LO) / (NZ - 1) as f64;
    for i in 0..NZ {
        let z = Z_LO + i as f64 * dz;
        let f = (-0.5 * z * z).exp();
        v.push((z, f));
        sum += f * dz;
    }
    for p in v.iter_mut() {
        p.1 /= sum;
    }
    v
}

/// Pearson IV rebuilt by integrating ln f' = -(z+a)/D, with D = b0 + b1 z + b2 z^2.
/// Returns None when the moments do not admit a valid Pearson IV.
fn pearson_iv_grid(skew: f64, kurt: f64) -> Option<Vec<(f64, f64)>> {
    let g = skew;
    let b = kurt;
    let delta = 10.0 * b - 12.0 * g * g - 18.0;
    if !(delta.is_finite() && delta > 1e-9) {
        return None;
    }
    // Pearson parameterisation with p'/p = (a - z) / (b0 + b1 z + b2 z^2).
    // Check against the Gaussian limit (gamma = 0, beta = 3): a = 0, b0 = 1,
    // b2 = 0, so p'/p = -z and p = exp(-z^2/2).
    let a = g * (b + 3.0) / delta;
    let b0 = (4.0 * b - 3.0 * g * g) / delta;
    let b1 = a;
    let b2 = (2.0 * b - 3.0 * g * g - 6.0) / delta;
    if b2.abs() < 1e-9 {
        return None;
    }
    // D = b0 + b1 z + b2 z^2 must never vanish, which for b2 > 0 means a
    // negative discriminant. A real root would put a pole in the density.
    let disc = b1 * b1 - 4.0 * b0 * b2;
    if !(b2 > 0.0) || disc >= -1e-9 {
        return None;
    }
    let dz = (Z_HI - Z_LO) / (NZ - 1) as f64;
    let mut lnf = Vec::with_capacity(NZ);
    let mut acc = 0.0;
    let mut prev = {
        let z = Z_LO;
        (a - z) / (b0 + b1 * z + b2 * z * z)
    };
    lnf.push(0.0);
    for i in 1..NZ {
        let z = Z_LO + i as f64 * dz;
        let d = b0 + b1 * z + b2 * z * z;
        if !(d.is_finite() && d.abs() > 1e-12) {
            return None;
        }
        let cur = (a - z) / d;
        acc += 0.5 * (prev + cur) * dz;
        prev = cur;
        lnf.push(acc);
    }
    let mx = lnf.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !mx.is_finite() {
        return None;
    }
    let mut v: Vec<(f64, f64)> = Vec::with_capacity(NZ);
    let mut sum = 0.0;
    for i in 0..NZ {
        let z = Z_LO + i as f64 * dz;
        let f = (lnf[i] - mx).exp();
        if !f.is_finite() {
            return None;
        }
        v.push((z, f));
        sum += f * dz;
    }
    if !(sum.is_finite() && sum > 0.0) {
        return None;
    }
    for p in v.iter_mut() {
        p.1 /= sum;
    }
    Some(v)
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;

    #[test]
    fn gaussian_limit_reproduces_unit_variance() {
        let p = Profile::from_moments(Nm::new(100.0), Nm::new(30.0), 0.0, 3.0).unwrap();
        // (0, 3) is the degenerate case: the Pearson IV collapses to a Gaussian,
        // so the Gaussian path is taken and labelled as such.
        assert_eq!(p.shape, Shape::GaussianFallback);
        let (skew, kurt) = p.computed_moments();
        assert!(skew.abs() < 1e-6, "skew {skew}");
        assert!((kurt - 3.0).abs() < 1e-3, "kurtosis {kurt}");
        // and it matches exp(-z^2/2)/sqrt(2pi)
        let expect = (-0.5f64).exp() / (2.0 * std::f64::consts::PI).sqrt();
        assert!((p.density(1.0) - expect).abs() < 1e-4, "{} vs {expect}", p.density(1.0));
    }

    #[test]
    fn pearson_iv_reproduces_the_input_moments() {
        // Boron-like asymmetry: negative skew, kurtosis above 3.
        for (skew, kurt) in [(-0.95, 5.0), (-1.09, 5.6), (0.4, 3.6)] {
            let p = Profile::from_moments(Nm::new(100.0), Nm::new(36.0), skew, kurt).unwrap();
            assert_eq!(p.shape, Shape::PearsonIv, "moments ({skew},{kurt}) should admit a Pearson IV");
            // The density is rebuilt by numerical integration, so the achieved
            // moments track the requested ones rather than reproducing them
            // exactly. Ten percent on skew is the observed agreement.
            let (s, k) = p.computed_moments();
            assert!((s - skew).abs() < 0.12, "skew {s} vs {skew}");
            assert!((k - kurt).abs() < 0.35, "kurtosis {k} vs {kurt}");
        }
    }

    #[test]
    fn impossible_moments_fall_back_to_a_gaussian() {
        // delta <= 0: kurtosis too low for the skew.
        let p = Profile::from_moments(Nm::new(100.0), Nm::new(30.0), 2.0, 3.0).unwrap();
        assert_eq!(p.shape, Shape::GaussianFallback);
        let (_, k) = p.computed_moments();
        assert!((k - 3.0).abs() < 1e-3);
    }

    #[test]
    fn dose_is_conserved() {
        for (s, e, dose) in [(Species::Boron, 30.0, 1.0e14), (Species::Arsenic, 60.0, 5.0e15), (Species::Phosphorus, 100.0, 1.0e13)] {
            let p = Profile::from_species(s, e).unwrap();
            let got = p.integrated_dose_cm2(dose);
            // Everything below the surface is retained; the part of the profile
            // above it is reflected. So the integral equals dose x (1 - reflected).
            let expect = dose * (1.0 - p.reflected_fraction());
            let rel = (got - expect).abs() / expect;
            assert!(rel < 5e-3, "{} {} keV: integrated {got:.4e} vs {expect:.4e} expected ({rel:.2e})", s.label(), e);
            assert!(got <= dose * 1.001, "never more than the requested dose");
        }
    }

    #[test]
    fn no_dopant_above_the_surface() {
        let p = Profile::from_species(Species::Boron, 10.0).unwrap();
        assert_eq!(p.concentration_cm3(Nm::new(-1.0), 1e14), 0.0);
        assert_eq!(p.concentration_cm3(Nm::new(-1e6), 1e14), 0.0);
        // at the surface the profile is truncated, not wrapped
        let peak = (0..200).map(|i| p.concentration_cm3(Nm::new(i as f64 * 0.5), 1e14)).fold(0.0, f64::max);
        assert!(p.concentration_cm3(Nm::new(0.0), 1e14) < peak, "the surface value is below the peak");
    }

    #[test]
    fn junction_depth_matches_a_hand_calculation_for_a_gaussian() {
        // Gaussian: xj = Rp + dRp*sqrt(2 ln(peak/bg)).
        let dose: f64 = 1.0e14;
        let p = Profile::from_moments(Nm::new(100.0), Nm::new(55.7), 0.0, 3.0).unwrap();
        let bg: f64 = 1.0e15;
        let peak: f64 = dose / (2.5066282746310002 * 55.7 * 1e-7);
        let expect = 100.0_f64 + 55.7 * (2.0 * (peak / bg).ln()).sqrt();
        let got = p.junction_depth_nm(PerCm2::new(dose), bg).unwrap();
        assert!((got - expect).abs() / expect < 0.02, "got {got:.1}, hand {expect:.1}");
    }

    #[test]
    fn no_junction_when_the_peak_is_below_the_background() {
        let p = Profile::from_species(Species::Boron, 10.0).unwrap();
        assert!(p.junction_depth_nm(PerCm2::new(1.0e10), 1.0e18).is_none());
    }

    #[test]
    fn asymmetry_shifts_the_tail_deeper_than_a_gaussian() {
        // A negative skew means a tail toward the surface; a positive one goes deep.
        let sym = Profile::from_moments(Nm::new(100.0), Nm::new(36.0), 0.0, 3.0).unwrap();
        let deep = Profile::from_moments(Nm::new(100.0), Nm::new(36.0), 0.6, 3.8).unwrap();
        assert_eq!(deep.shape, Shape::PearsonIv);
        let d = 1.0e14;
        assert!(deep.junction_depth_nm(PerCm2::new(d), 1e16).unwrap() > sym.junction_depth_nm(PerCm2::new(d), 1e16).unwrap());
    }
}

// ---------------------------------------------------------------------------
// Writing a profile into the dopant field (build PRD D-01, D-03)
// ---------------------------------------------------------------------------

use atto_dopant::{DopantField, DopantRecord, Type};
use atto_geom::column::ColumnGrid;
use atto_geom::svo::Material;

/// Result of an implant into the field.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ImplantStats {
    pub columns: usize,
    pub leaves_written: usize,
    /// Deepest depth at which dopant was stored, nm.
    pub deepest_nm: f64,
    /// Fraction of the dose that landed in silicon rather than in the air above
    /// the surface or in a covering film.
    pub retained_fraction: f64,
}

/// Write an implant profile into the dopant field at leaf resolution.
///
/// Depth is measured from the top of the silicon in each column, so an implant
/// performed before any film is deposited lands where the profile says. A
/// covering film is *not* accounted for: there is no stopping-power model yet
/// (audit F-11 flagged the arbitrary factors the earlier build used for exactly this), so the
/// caller should implant before coating.
pub fn implant_into_field(
    field: &mut DopantField,
    grid: &ColumnGrid,
    profile: &Profile,
    dose_cm2: f64,
    species: Species,
    kind: Type,
) -> ImplantStats {
    // no mask awareness and no lateral spread: every column receives the full
    // beam straight down. `implant_into_field_masked` is the physical one.
    let open = vec![true; grid.cols.len()];
    implant_into_field_masked(field, grid, profile, dose_cm2, species, kind, &open, 0.0)
}

/// Active concentration after clustering, cm^-3 (build PRD D-12).
///
/// Above a solubility limit the excess dopant forms inactive clusters rather
/// than adding carriers, so the active concentration saturates while the
/// chemical concentration keeps rising. A dopant that does not cluster (limit at
/// or above the peak) is returned unchanged.
pub fn active_concentration(total_cm3: f64, solubility_cm3: f64) -> f64 {
    if !(total_cm3 > 0.0) || !total_cm3.is_finite() {
        return 0.0;
    }
    if !(solubility_cm3 > 0.0) || !solubility_cm3.is_finite() || total_cm3 <= solubility_cm3 {
        return total_cm3;
    }
    // the excess clusters with a residual active fraction, so active grows slowly
    let excess = total_cm3 - solubility_cm3;
    solubility_cm3 + excess * 0.1
}

/// Fraction of the dose that is electrically active at a given peak.
pub fn active_fraction(total_cm3: f64, solubility_cm3: f64) -> f64 {
    if !(total_cm3 > 0.0) {
        return 0.0;
    }
    (active_concentration(total_cm3, solubility_cm3) / total_cm3).clamp(0.0, 1.0)
}

/// Columns that are protected from the beam by a mask.
pub fn mask_from_top(grid: &ColumnGrid) -> Vec<bool> {
    (0..grid.cols.len())
        .map(|k| {
            if !grid.active[k] {
                return false;
            }
            matches!(
                atto_geom::column::solid_top(&grid.cols[k]),
                Some((_, Material::Photoresist)) | Some((_, Material::SiliconNitride))
            )
        })
        .collect()
}

/// Implant with lateral straggle and mask shadowing (build PRD D-06).
///
/// The incident dose is spread laterally by a Gaussian of `sigma_lat_nm` before
/// it is deposited, so a masked column still receives dopant from a nearby
/// opening while a column far under the mask receives none. `open[k] == false`
/// means a mask covers that column.
#[allow(clippy::too_many_arguments)]
pub fn implant_into_field_masked(
    field: &mut DopantField,
    grid: &ColumnGrid,
    profile: &Profile,
    dose_cm2: f64,
    species: Species,
    kind: Type,
    open: &[bool],
    sigma_lat_nm: f64,
) -> ImplantStats {
    if open.len() != grid.cols.len() {
        return ImplantStats::default();
    }
    // the dose each column receives, before and after the lateral spread
    let mut received: Vec<f64> = (0..grid.cols.len()).map(|k| if open[k] && grid.active[k] { 1.0 } else { 0.0 }).collect();
    if sigma_lat_nm > 0.0 {
        received = blur(&received, grid.nx, grid.ny, sigma_lat_nm / grid.pitch);
    }
    let leaf = field.leaf_nm();
    // The cutoff must sit *below* the substrate background, otherwise the field
    // stops storing dopant while it is still able to flip the type and the
    // junction is truncated short of where it really is.
    let peak = profile.concentration_cm3(Nm::new(profile.rp_nm), dose_cm2);
    let bg = field.background().map(|b| b.concentration_cm3).unwrap_or(0.0);
    let floor = if bg > 0.0 { (peak * 1e-3).min(bg * 1e-2) } else { peak * 1e-3 };
    let max_depth = profile.rp_nm + 40.0 * profile.drp_nm;
    let mut columns = 0;
    let mut leaves = 0;
    let mut deepest = 0.0f64;
    let mut stored_dose = 0.0f64;
    let area_cm2 = grid.cell_area() * 1.0e-14;

    for iy in 0..grid.ny {
        for ix in 0..grid.nx {
            let k = grid.idx(ix, iy);
            if !grid.active[k] {
                continue;
            }
            let share = received[k];
            if share <= 1.0e-6 {
                continue;
            }
            let dose_here = dose_cm2 * share;
            let col = &grid.cols[k];
            let Some(si) = col.iter().rev().find(|r| r.material == Material::Silicon) else {
                continue;
            };
            let surface = si.z1;
            let (x, y) = (grid.x0 + (ix as f64 + 0.5) * grid.pitch, grid.y0 + (iy as f64 + 0.5) * grid.pitch);
            let mut wrote = false;
            let mut d = leaf * 0.5;
            while surface - d > si.z0 && d < max_depth {
                let c = profile.concentration_cm3(Nm::new(d), dose_here);
                if c <= floor {
                    if d > profile.rp_nm {
                        break;
                    }
                    d += leaf;
                    continue;
                }
                let key = field.key_for(x, y, surface - d);
                field.add(key, DopantRecord::new(species, kind, c));
                stored_dose += c * leaf * 1.0e-7;
                leaves += 1;
                deepest = deepest.max(d);
                wrote = true;
                d += leaf;
            }
            if wrote {
                columns += 1;
            }
        }
    }
    let expected = dose_cm2 * columns as f64 * area_cm2;
    let _ = &received;
    ImplantStats {
        columns,
        leaves_written: leaves,
        deepest_nm: deepest,
        retained_fraction: if expected > 0.0 { (stored_dose / expected).clamp(0.0, 1.0) } else { 0.0 },
    }
}

/// Separable Gaussian blur of a per-column quantity, in cell units.
fn blur(v: &[f64], nx: usize, ny: usize, sigma_cells: f64) -> Vec<f64> {
    if sigma_cells <= 0.0 || nx == 0 || ny == 0 {
        return v.to_vec();
    }
    let r = (3.0 * sigma_cells).ceil() as isize;
    let kern: Vec<f64> = (-r..=r)
        .map(|i| (-(i as f64) * (i as f64) / (2.0 * sigma_cells * sigma_cells)).exp())
        .collect();
    let ksum: f64 = kern.iter().sum();
    let kern: Vec<f64> = kern.iter().map(|k| k / ksum).collect();
    let mut tmp = vec![0.0; v.len()];
    for y in 0..ny {
        for x in 0..nx {
            let mut acc = 0.0;
            for (i, w) in kern.iter().enumerate() {
                let xx = x as isize + i as isize - r;
                if xx >= 0 && (xx as usize) < nx {
                    acc += w * v[y * nx + xx as usize];
                }
            }
            tmp[y * nx + x] = acc;
        }
    }
    let mut out = vec![0.0; v.len()];
    for y in 0..ny {
        for x in 0..nx {
            let mut acc = 0.0;
            for (i, w) in kern.iter().enumerate() {
                let yy = y as isize + i as isize - r;
                if yy >= 0 && (yy as usize) < ny {
                    acc += w * tmp[yy as usize * nx + x];
                }
            }
            out[y * nx + x] = acc;
        }
    }
    out
}

#[cfg(test)]
mod field_tests {
    use atto_units::{Mm, Um};
    use super::*;
    use atto_dopant::DopantField;
    use atto_geom::column::ProcessWindow;
    use atto_geom::svo::Wafer;

    fn setup() -> (DopantField, ColumnGrid) {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let g = ColumnGrid::sample(&w, &ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] }).unwrap();
        let mut f = DopantField::new(atto_units::Nm::new(16.0), w.bounds.min);
        f.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e15));
        (f, g)
    }

    #[test]
    fn the_geometric_junction_agrees_with_the_analytic_one() {
        let (mut f, g) = setup();
        // A junction needs opposite types: donors into the p background.
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        let dose = 5.0e14;
        let stats = implant_into_field(&mut f, &g, &p, dose, Species::Phosphorus, Type::Donor);
        assert!(stats.columns > 0);
        let analytic = p.junction_depth_nm(PerCm2::new(dose), 1.0e15).expect("analytic junction");
        let surface = 500_000.0;
        let col = &g.cols[g.idx(8, 8)];
        let geo = f.junction_nm(col, g.x0 + 8.5 * g.pitch, g.y0 + 8.5 * g.pitch, surface).expect("geometric junction");
        // The field stores one value per leaf, so the junction is quantised to
        // the leaf; that is the resolution limit, not an error in the model.
        assert!((geo - analytic).abs() <= 16.0, "field {geo:.1} nm vs analytic {analytic:.1} nm");
    }

    #[test]
    fn the_field_reports_the_implanted_type_above_the_junction() {
        let (mut f, g) = setup();
        let p = Profile::from_species(Species::Boron, 30.0).unwrap();
        implant_into_field(&mut f, &g, &p, 1.0e14, Species::Boron, Type::Acceptor);
        let surface = 500_000.0;
        let x = g.x0 + 8.5 * g.pitch;
        let y = g.y0 + 8.5 * g.pitch;
        // boron is an acceptor: net doping goes negative at the peak
        assert_eq!(f.type_at(x, y, surface - p.rp_nm), Some(Type::Acceptor));
        assert!(f.net_at(x, y, surface - p.rp_nm) < -1.0e17);
    }

    #[test]
    fn a_phosphorus_implant_creates_a_pn_junction_in_a_p_substrate() {
        let (mut f, g) = setup();
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        implant_into_field(&mut f, &g, &p, 5.0e14, Species::Phosphorus, Type::Donor);
        let surface = 500_000.0;
        let x = g.x0 + 8.5 * g.pitch;
        let y = g.y0 + 8.5 * g.pitch;
        assert_eq!(f.type_at(x, y, surface - p.rp_nm), Some(Type::Donor), "n above");
        assert_eq!(f.type_at(x, y, surface - 2000.0), Some(Type::Acceptor), "p below");
        assert!(f.junction_nm(&g.cols[g.idx(8, 8)], x, y, surface).is_some());
    }

    #[test]
    fn a_weak_implant_leaves_the_substrate_type_unchanged() {
        let (mut f, g) = setup();
        let p = Profile::from_species(Species::Boron, 10.0).unwrap();
        // 1e9 cm^-2 peaks around 2.5e14, below the 1e15 background: no junction
        implant_into_field(&mut f, &g, &p, 1.0e9, Species::Boron, Type::Donor);
        let surface = 500_000.0;
        let col = &g.cols[g.idx(8, 8)];
        assert_eq!(f.junction_nm(col, g.x0 + 8.5 * g.pitch, g.y0 + 8.5 * g.pitch, surface), None);
    }

    /// A grid with a four-column opening in the middle and a mask elsewhere.
    fn masked_setup() -> (DopantField, ColumnGrid) {
        let (f, g) = setup();
        let mut m = g.clone();
        for iy in 0..m.ny {
            for ix in 0..m.nx {
                if !(10..14).contains(&ix) {
                    let k = m.idx(ix, iy);
                    let top = atto_geom::column::solid_top(&m.cols[k]).map(|(t, _)| t).unwrap_or(0.0);
                    atto_geom::column::set_interval(&mut m.cols[k], top, top + 96.0, Material::Photoresist);
                }
            }
        }
        (f, m)
    }

    #[test]
    fn a_mask_blocks_the_beam_when_there_is_no_lateral_spread() {
        let (mut f, g) = masked_setup();
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        let open = mask_from_top(&g).iter().map(|m| !m).collect::<Vec<_>>();
        let st = implant_into_field_masked(&mut f, &g, &p, 5.0e14, Species::Phosphorus, Type::Donor, &open, 0.0);
        // only the opening is implanted
        assert_eq!(st.columns, 4 * g.ny, "four columns wide, every row");
        let surface = 500_000.0;
        let x_open = g.x0 + 11.5 * g.pitch;
        let x_mask = g.x0 + 20.5 * g.pitch;
        let y = g.y0 + 8.5 * g.pitch;
        assert!(f.net_at(x_open, y, surface - 40.0) > 1.0e18, "the opening is doped");
        assert!(f.net_at(x_mask, y, surface - 40.0).abs() < 2.0e15, "the mask stops the beam");
    }

    #[test]
    fn lateral_straggle_dopes_under_the_mask_edge_but_not_far_in() {
        let (mut f, g) = masked_setup();
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        let open = mask_from_top(&g).iter().map(|m| !m).collect::<Vec<_>>();
        implant_into_field_masked(&mut f, &g, &p, 5.0e14, Species::Phosphorus, Type::Donor, &open, 60.0);
        let surface = 500_000.0;
        let y = g.y0 + 8.5 * g.pitch;
        let near = f.net_at(g.x0 + 8.5 * g.pitch, y, surface - 40.0);
        let far = f.net_at(g.x0 + 30.5 * g.pitch, y, surface - 40.0);
        assert!(near > 1.0e17, "one cell under the edge is doped: {near:.2e}");
        assert!(far.abs() < 2.0e15, "far under the mask is not: {far:.2e}");
    }

    #[test]
    fn the_dopant_falls_off_with_distance_under_the_mask() {
        let (mut f, g) = masked_setup();
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        let open = mask_from_top(&g).iter().map(|m| !m).collect::<Vec<_>>();
        implant_into_field_masked(&mut f, &g, &p, 5.0e14, Species::Phosphorus, Type::Donor, &open, 60.0);
        let surface = 500_000.0;
        let y = g.y0 + 8.5 * g.pitch;
        let mut prev = f64::INFINITY;
        for ix in [8usize, 6, 4, 2] {
            let v = f.net_at(g.x0 + (ix as f64 + 0.5) * g.pitch, y, surface - 40.0);
            assert!(v < prev, "dopant must fall away from the opening: ix {ix} gives {v:.2e} after {prev:.2e}");
            prev = v;
        }
    }

    #[test]
    fn the_junction_under_the_edge_is_shallower_than_in_the_opening() {
        let (mut f, g) = masked_setup();
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        let open = mask_from_top(&g).iter().map(|m| !m).collect::<Vec<_>>();
        implant_into_field_masked(&mut f, &g, &p, 5.0e14, Species::Phosphorus, Type::Donor, &open, 60.0);
        let surface = 500_000.0;
        let y = g.y0 + 8.5 * g.pitch;
        let xj = |ix: usize| f.junction_nm(&g.cols[g.idx(ix, 8)], g.x0 + (ix as f64 + 0.5) * g.pitch, y, surface);
        let inside = xj(11).expect("a junction in the opening");
        let edge = xj(9).expect("a junction just under the edge");
        assert!(edge < inside, "the straggled junction is shallower: {edge:.0} vs {inside:.0} nm");
    }

    #[test]
    fn an_aligned_implant_channels_deeper_than_a_tilted_one() {
        let aligned = Profile::from_species(Species::Phosphorus, 60.0).unwrap().with_tilt(0.0);
        let tilted = Profile::from_species(Species::Phosphorus, 60.0).unwrap().with_tilt(7.0);
        let d = 5.0e14;
        let a = aligned.concentration_cm3(Nm::new(400.0), d);
        let b = tilted.concentration_cm3(Nm::new(400.0), d);
        assert!(a > 100.0 * b, "the aligned tail must be far deeper: {a:.2e} vs {b:.2e}");
        // and the near-surface peak is not much changed
        let pa = aligned.concentration_cm3(Nm::new(50.0), d);
        let pb = tilted.concentration_cm3(Nm::new(50.0), d);
        assert!((pa - pb).abs() / pa < 0.4, "the peak barely moves: {pa:.2e} vs {pb:.2e}");
    }

    #[test]
    fn the_channelling_tail_decays_exponentially() {
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap().with_tilt(0.0);
        let d = 5.0e14;
        // two equal steps deep in the tail must give the same ratio
        let a = p.concentration_cm3(Nm::new(300.0), d) - p.concentration_cm3(Nm::new(300.0), d) * 0.0;
        let b = p.concentration_cm3(Nm::new(360.0), d);
        let c = p.concentration_cm3(Nm::new(420.0), d);
        let r1 = b / a;
        let r2 = c / b;
        assert!(r1 > 0.0 && r1 < 1.0, "the tail must decay: {a:.2e} {b:.2e}");
        assert!((r1 - r2).abs() / r1 < 0.25, "equal steps, equal ratio: {r1:.4} vs {r2:.4}");
    }

    #[test]
    fn tilting_further_always_reduces_the_tail() {
        let d = 5.0e14;
        let mut prev = f64::INFINITY;
        for tilt in [0.0, 0.5, 1.0, 3.0, 7.0] {
            let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap().with_tilt(tilt);
            let c = p.concentration_cm3(Nm::new(400.0), d);
            assert!(c < prev, "the tail must shrink with tilt: {tilt} deg gives {c:.2e} after {prev:.2e}");
            prev = c;
        }
    }

    #[test]
    fn no_channelling_leaves_the_profile_alone() {
        let plain = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap().with_tilt(90.0);
        assert_eq!(p.channel_fraction, 0.0);
        for depth in [20.0, 100.0, 300.0] {
            assert!((p.concentration_cm3(Nm::new(depth), 5.0e14) - plain.concentration_cm3(Nm::new(depth), 5.0e14)).abs() < 1e-9);
        }
    }

    #[test]
    fn a_screen_oxide_swallows_the_dose_that_stops_in_it() {
        let p = Profile::from_species(Species::Phosphorus, 60.0).unwrap();
        // no oxide: everything reaches the silicon
        assert!((p.retained_in_silicon(0.0) - 1.0).abs() < 1e-12);
        // a thin screen oxide takes very little
        let thin = p.retained_in_silicon(5.0);
        assert!(thin > 0.9, "5 nm of oxide retains {thin:.3}");
        // a thick one takes most of a shallow implant
        let thick = p.retained_in_silicon(120.0);
        assert!(thick < 0.6, "120 nm of oxide retains only {thick:.3}");
        // and it falls monotonically
        let mut prev = 1.0;
        for d in [0.0, 10.0, 30.0, 60.0, 120.0, 200.0] {
            let r = p.retained_in_silicon(d);
            assert!(r <= prev + 1e-9, "retention rose at {d} nm: {prev:.3} -> {r:.3}");
            prev = r;
        }
    }

    #[test]
    fn clustering_saturates_the_active_concentration() {
        let sol = 1.0e20;
        // below the limit nothing clusters
        assert!((active_concentration(5.0e19, sol) - 5.0e19).abs() < 1e-6);
        // above it the active part grows far more slowly than the total
        let a1 = active_concentration(2.0e20, sol);
        let a2 = active_concentration(4.0e20, sol);
        assert!(a1 < 2.0e20, "the active part is below the chemical one");
        assert!((a2 - a1) < 0.5 * (4.0e20 - 2.0e20), "and grows slower than it");
        // the fraction falls as the dose rises
        let f1 = active_fraction(2.0e20, sol);
        let f2 = active_fraction(4.0e20, sol);
        assert!(f2 < f1, "active fraction {f1:.3} -> {f2:.3}");
    }

    #[test]
    fn a_dopant_that_does_not_cluster_is_unchanged() {
        assert!((active_concentration(1.0e21, 1.0e22) - 1.0e21).abs() < 1e-6);
        assert!((active_fraction(1.0e21, 1.0e22) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn bad_clustering_input_is_not_a_number() {
        assert_eq!(active_concentration(-1.0, 1.0e20), 0.0);
        assert_eq!(active_concentration(f64::NAN, 1.0e20), 0.0);
        assert_eq!(active_fraction(0.0, 1.0e20), 0.0);
    }

    #[test]
    fn the_field_is_sparse() {
        let (mut f, g) = setup();
        let p = Profile::from_species(Species::Boron, 30.0).unwrap();
        let stats = implant_into_field(&mut f, &g, &p, 1.0e14, Species::Boron, Type::Acceptor);
        // one record per leaf in the doped depth range, not per voxel
        let per_column = stats.leaves_written as f64 / stats.columns as f64;
        assert!(per_column < 20.0, "{per_column:.1} leaves per column");
        assert_eq!(f.occupied(), stats.leaves_written);
    }
}
