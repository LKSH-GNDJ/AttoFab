// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Chemical-mechanical planarisation (proposed requirements P-01, P-02).
//!
//! **P-01, Preston removal.** `dh/dt = K * p * v`, with the down-force `p` and the
//! relative velocity `v` as recipe parameters. Preston (1927) is the published
//! model, so this part is checked against the closed form.
//!
//! **P-02, pattern-density dependence.** A line scan (one dimension, along the
//! wafer) is polished in time steps. The pad carries the whole down-force, so
//! it presses harder on whatever is highest: the pressure on a cell is
//! `p * w_i / mean(w)`, where `w_i = exp(-step_i / h_c)` is the contact weight
//! of a cell that sits `step_i` below the highest cell within one
//! planarisation length, and the mean runs over that same window. Force is
//! conserved inside the window. Narrow soft features dish a little and stop,
//! because they fall out of contact; a feature wider than the window never
//! sees a higher neighbour and keeps dishing; a dense soft array leaves a
//! small area of hard dielectric to carry the load and erodes.
//!
//! **Maturity: L1 for P-02.** The mechanism is the published qualitative one
//! (Ouma-style planarisation length), but the contact length and the
//! planarisation length are registered as *assumed*: no measured CMP profile
//! exists in the tree, so the feature is capped at L1 (FEATURES.md section 6).

use atto_units::{KPa, Minutes, MetresPerSecond, Nm};

/// Which film a cell of the line scan is made of.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CmpMaterial {
    /// The soft, fast-polishing film (a damascene metal).
    Metal,
    /// The hard film that carries the load once the metal is recessed.
    Dielectric,
}

/// Why a polish could not run.
#[derive(Debug, Clone, PartialEq)]
pub enum CmpError {
    /// The scan has no cells, or the cell pitch is not positive.
    EmptyScan,
    /// A recipe value was negative, zero where it must not be, or not finite.
    BadRecipe(&'static str),
}

impl std::fmt::Display for CmpError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CmpError::EmptyScan => write!(f, "the line scan has no cells or a non-positive pitch"),
            CmpError::BadRecipe(what) => write!(f, "invalid polish recipe: {what}"),
        }
    }
}

impl std::error::Error for CmpError {}

/// Preston coefficient of the metal, nm/min per (kPa * m/s). Assumed.
pub const K_METAL: f64 = 35.0;
/// Preston coefficient of the dielectric, nm/min per (kPa * m/s). Assumed.
pub const K_DIELECTRIC: f64 = 10.0;
/// Pad contact length `h_c`, nm: the step height over which a cell leaves contact. Assumed.
pub const CONTACT_LENGTH_NM: f64 = 15.0;
/// Planarisation length, nm: the radius over which the pad shares load. Assumed.
pub const PLANARISATION_LENGTH_NM: f64 = 20_000.0;

/// Preston removal rate and the pattern-density parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CmpModel {
    /// Preston coefficient of the metal, nm/min per (kPa * m/s).
    pub k_metal: f64,
    /// Preston coefficient of the dielectric, nm/min per (kPa * m/s).
    pub k_dielectric: f64,
    /// Pad contact length.
    pub contact_length: Nm,
    /// Planarisation length.
    pub planarisation_length: Nm,
}

impl CmpModel {
    /// The bundled, assumed parameter set.
    pub fn bundled() -> Self {
        CmpModel {
            k_metal: K_METAL,
            k_dielectric: K_DIELECTRIC,
            contact_length: Nm(CONTACT_LENGTH_NM),
            planarisation_length: Nm(PLANARISATION_LENGTH_NM),
        }
    }

    /// Preston coefficient of a material.
    pub fn k(&self, m: CmpMaterial) -> f64 {
        match m {
            CmpMaterial::Metal => self.k_metal,
            CmpMaterial::Dielectric => self.k_dielectric,
        }
    }

    /// Material removed from a flat, uniformly loaded film: `K * p * v * t` (P-01).
    pub fn removal(&self, m: CmpMaterial, p: KPa, v: MetresPerSecond, t: Minutes) -> Nm {
        Nm(self.k(m) * p.0 * v.0 * t.0)
    }

    /// Removal of a two-film stack: the polish clears the top film, then meets the stop layer.
    ///
    /// Over-polish is declared as a fraction of the time it took to clear the top film. Past
    /// that point the pad works on the stop layer at its own Preston rate, which selectivity
    /// makes small: a 30% over-polish removes `0.3 x k_stop / k_top` of what it removed before.
    /// Nothing else limits it, because nothing else is in the model — this is a two-film stack
    /// with a declared over-polish, not a full stop-layer model with endpoint detection.
    pub fn removal_on_stop(
        &self,
        top: CmpMaterial,
        stop: CmpMaterial,
        p: KPa,
        v: MetresPerSecond,
        t_clear: Minutes,
        over_polish: f64,
    ) -> Result<(Nm, Nm), CmpError> {
        if !p.0.is_finite() || p.0 <= 0.0 {
            return Err(CmpError::BadRecipe("pressure must be finite and > 0"));
        }
        if !v.0.is_finite() || v.0 <= 0.0 {
            return Err(CmpError::BadRecipe("velocity must be finite and > 0"));
        }
        if !t_clear.0.is_finite() || t_clear.0 <= 0.0 {
            return Err(CmpError::BadRecipe("time to clear must be finite and > 0"));
        }
        if !over_polish.is_finite() || over_polish < 0.0 {
            return Err(CmpError::BadRecipe("over-polish must be finite and >= 0"));
        }
        let cleared = self.removal(top, p, v, t_clear);
        let over = self.removal(stop, p, v, Minutes(t_clear.0 * over_polish));
        Ok((cleared, over))
    }

    /// Polish a line scan for `t`, in `steps` equal explicit steps (P-02).
    pub fn polish(
        &self,
        scan: &LineScan,
        p: KPa,
        v: MetresPerSecond,
        t: Minutes,
        steps: usize,
    ) -> Result<LineScan, CmpError> {
        if scan.material.is_empty() || scan.material.len() != scan.height.len() || !(scan.dx.0 > 0.0) {
            return Err(CmpError::EmptyScan);
        }
        if !(p.0.is_finite() && p.0 >= 0.0) {
            return Err(CmpError::BadRecipe("pressure must be finite and not negative"));
        }
        if !(v.0.is_finite() && v.0 >= 0.0) {
            return Err(CmpError::BadRecipe("velocity must be finite and not negative"));
        }
        if !(t.0.is_finite() && t.0 >= 0.0) {
            return Err(CmpError::BadRecipe("time must be finite and not negative"));
        }
        if steps == 0 {
            return Err(CmpError::BadRecipe("steps must be at least one"));
        }
        if !(self.contact_length.0 > 0.0) || !(self.planarisation_length.0 > 0.0) {
            return Err(CmpError::BadRecipe("contact and planarisation lengths must be positive"));
        }
        let n = scan.material.len();
        let radius = ((self.planarisation_length.0 / scan.dx.0).round() as usize).max(1);
        let dt = t.0 / steps as f64;
        let mut h: Vec<f64> = scan.height.iter().map(|x| x.0).collect();
        let mut w = vec![0.0f64; n];
        for _ in 0..steps {
            for i in 0..n {
                let lo = i.saturating_sub(radius);
                let hi = (i + radius).min(n - 1);
                let mut top = f64::NEG_INFINITY;
                for j in lo..=hi {
                    if h[j] > top {
                        top = h[j];
                    }
                }
                w[i] = (-(top - h[i]) / self.contact_length.0).exp();
            }
            let mut next = h.clone();
            for i in 0..n {
                let lo = i.saturating_sub(radius);
                let hi = (i + radius).min(n - 1);
                let mut sum = 0.0;
                for j in lo..=hi {
                    sum += w[j];
                }
                let mean = sum / (hi - lo + 1) as f64;
                let local_p = p.0 * w[i] / mean;
                next[i] = h[i] - self.k(scan.material[i]) * local_p * v.0 * dt;
            }
            h = next;
        }
        Ok(LineScan { dx: scan.dx, material: scan.material.clone(), height: h.into_iter().map(Nm).collect() })
    }
}

/// A one-dimensional surface profile: a material and a height for each cell.
#[derive(Debug, Clone, PartialEq)]
pub struct LineScan {
    /// Cell pitch.
    pub dx: Nm,
    /// Material of each cell.
    pub material: Vec<CmpMaterial>,
    /// Surface height of each cell; a polish lowers it.
    pub height: Vec<Nm>,
}

impl LineScan {
    /// A flat scan built from consecutive segments of `(material, width)`.
    /// Each width is rounded to a whole number of cells, at least one.
    pub fn from_segments(dx: Nm, segments: &[(CmpMaterial, Nm)]) -> Self {
        let mut material = Vec::new();
        for (m, wdt) in segments {
            let cells = ((wdt.0 / dx.0).round() as usize).max(1);
            material.extend(std::iter::repeat(*m).take(cells));
        }
        let height = vec![Nm::ZERO; material.len()];
        LineScan { dx, material, height }
    }

    /// Mean surface height over cells `[i0, i1)`.
    pub fn mean_height(&self, i0: usize, i1: usize) -> Nm {
        let s: f64 = self.height[i0..i1].iter().map(|h| h.0).sum();
        Nm(s / (i1 - i0) as f64)
    }

    /// Dishing of a metal span `[i0, i1)`: the dielectric height just beyond the
    /// span minus the lowest metal height inside it. Positive when the metal is recessed.
    pub fn dishing(&self, i0: usize, i1: usize) -> Nm {
        let n = self.height.len();
        let mut shoulder = f64::NEG_INFINITY;
        if i0 > 0 {
            shoulder = shoulder.max(self.height[i0 - 1].0);
        }
        if i1 < n {
            shoulder = shoulder.max(self.height[i1].0);
        }
        let low = self.height[i0..i1].iter().map(|h| h.0).fold(f64::INFINITY, f64::min);
        Nm(shoulder - low)
    }

    /// Erosion of a dielectric region `[i0, i1)`: the height it started at (zero)
    /// minus its mean height now. Positive when the dielectric has thinned.
    pub fn erosion(&self, i0: usize, i1: usize) -> Nm {
        Nm(-self.mean_height(i0, i1).0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DX: Nm = Nm(500.0);
    const P: KPa = KPa(14.0);
    const V: MetresPerSecond = MetresPerSecond(1.0);

    #[test]
    fn removal_matches_preston_k_p_v_t() {
        let m = CmpModel::bundled();
        let got = m.removal(CmpMaterial::Metal, P, V, Minutes(0.5)).0;
        let want = K_METAL * 14.0 * 1.0 * 0.5;
        assert!((got - want).abs() / want < 0.02, "{got} vs {want}");
    }

    #[test]
    fn doubling_pressure_doubles_removal() {
        let m = CmpModel::bundled();
        let a = m.removal(CmpMaterial::Dielectric, KPa(7.0), V, Minutes(1.0)).0;
        let b = m.removal(CmpMaterial::Dielectric, KPa(14.0), V, Minutes(1.0)).0;
        assert!((b / a - 2.0).abs() < 1e-12);
    }

    #[test]
    fn a_flat_uniform_film_polishes_at_the_preston_rate() {
        let m = CmpModel::bundled();
        let scan = LineScan::from_segments(DX, &[(CmpMaterial::Dielectric, Nm(100_000.0))]);
        let out = m.polish(&scan, P, V, Minutes(0.2), 20).unwrap();
        let want = m.removal(CmpMaterial::Dielectric, P, V, Minutes(0.2)).0;
        for h in &out.height {
            assert!(((-h.0) - want).abs() / want < 0.02);
        }
    }

    fn pad_scan(width_nm: f64) -> (LineScan, usize, usize) {
        let field = Nm(60_000.0);
        let scan = LineScan::from_segments(
            DX,
            &[
                (CmpMaterial::Dielectric, field),
                (CmpMaterial::Metal, Nm(width_nm)),
                (CmpMaterial::Dielectric, field),
            ],
        );
        let i0 = (field.0 / DX.0) as usize;
        let i1 = i0 + (width_nm / DX.0).round().max(1.0) as usize;
        (scan, i0, i1)
    }

    #[test]
    fn a_wide_pad_dishes_more_than_a_narrow_line() {
        let m = CmpModel::bundled();
        let (narrow, a0, a1) = pad_scan(1_000.0);
        let (wide, b0, b1) = pad_scan(50_000.0);
        let n = m.polish(&narrow, P, V, Minutes(0.1), 200).unwrap().dishing(a0, a1).0;
        let w = m.polish(&wide, P, V, Minutes(0.1), 200).unwrap().dishing(b0, b1).0;
        assert!(w > 2.0 * n, "wide pad dished {w} nm, narrow line {n} nm");
    }

    fn array_scan(metal_fraction: f64) -> (LineScan, usize, usize) {
        let pitch = 2_000.0;
        let metal = pitch * metal_fraction;
        let mut seg = vec![(CmpMaterial::Dielectric, Nm(60_000.0))];
        for _ in 0..20 {
            seg.push((CmpMaterial::Metal, Nm(metal)));
            seg.push((CmpMaterial::Dielectric, Nm(pitch - metal)));
        }
        seg.push((CmpMaterial::Dielectric, Nm(60_000.0)));
        let scan = LineScan::from_segments(DX, &seg);
        let n = scan.material.len();
        // The dielectric cells in the middle of the array.
        let mid = n / 2;
        (scan, mid - 8, mid + 8)
    }

    fn mean_dielectric_height(s: &LineScan, i0: usize, i1: usize) -> f64 {
        let mut sum = 0.0;
        let mut c = 0;
        for i in i0..i1 {
            if s.material[i] == CmpMaterial::Dielectric {
                sum += s.height[i].0;
                c += 1;
            }
        }
        sum / c as f64
    }

    #[test]
    fn erosion_grows_with_pattern_density() {
        let m = CmpModel::bundled();
        let mut prev = -1.0;
        for frac in [0.25, 0.5, 0.75] {
            let (scan, i0, i1) = array_scan(frac);
            let out = m.polish(&scan, P, V, Minutes(0.15), 300).unwrap();
            let e = -mean_dielectric_height(&out, i0, i1);
            assert!(e > prev, "erosion {e} nm at density {frac} did not exceed {prev}");
            prev = e;
        }
    }

    #[test]
    fn polishing_is_deterministic() {
        let m = CmpModel::bundled();
        let (scan, _, _) = array_scan(0.5);
        let a = m.polish(&scan, P, V, Minutes(0.1), 50).unwrap();
        let b = m.polish(&scan, P, V, Minutes(0.1), 50).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn a_bad_recipe_is_refused_not_clamped() {
        let m = CmpModel::bundled();
        let (scan, _, _) = pad_scan(1_000.0);
        assert!(matches!(m.polish(&scan, KPa(-1.0), V, Minutes(1.0), 10), Err(CmpError::BadRecipe(_))));
        assert!(matches!(m.polish(&scan, P, V, Minutes(f64::NAN), 10), Err(CmpError::BadRecipe(_))));
        assert!(matches!(m.polish(&scan, P, V, Minutes(1.0), 0), Err(CmpError::BadRecipe(_))));
        let empty = LineScan { dx: DX, material: vec![], height: vec![] };
        assert_eq!(m.polish(&empty, P, V, Minutes(1.0), 10), Err(CmpError::EmptyScan));
    }
}

#[cfg(test)]
mod stop_layer_tests {
    use super::*;

    fn model() -> CmpModel { CmpModel::bundled() }

    #[test]
    fn polish_stops_at_a_dissimilar_film() {
        // the stop layer polishes far slower than the film above it
        let m = model();
        let (top, stop) = m
            .removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, KPa::new(20.0), MetresPerSecond::new(0.5), Minutes::new(1.0), 0.3)
            .expect("a two-film polish");
        assert!(top.0 > 0.0, "the top film is removed: {}", top.0);
        assert!(stop.0 < top.0 * 0.5, "the stop layer barely moves: {} against {}", stop.0, top.0);
    }

    #[test]
    fn thirty_percent_over_polish_removes_only_the_selectivity_allowance() {
        // P-03's acceptance: over-polish of 30% removes no more of the stop layer than the
        // selectivity allows, and exactly that much in this model
        let m = model();
        let (p, v, t) = (KPa::new(20.0), MetresPerSecond::new(0.5), Minutes::new(1.0));
        let (top, stop) = m
            .removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, p, v, t, 0.3)
            .expect("a two-film polish");
        let selectivity = m.k(CmpMaterial::Dielectric) / m.k(CmpMaterial::Metal);
        let allowed = top.0 * 0.3 * selectivity;
        assert!((stop.0 - allowed).abs() < 1e-9, "stop removal {} must equal the allowance {}", stop.0, allowed);
    }

    #[test]
    fn no_over_polish_touches_nothing_of_the_stop_layer() {
        let m = model();
        let (_, stop) = m
            .removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, KPa::new(20.0), MetresPerSecond::new(0.5), Minutes::new(1.0), 0.0)
            .expect("a two-film polish");
        assert_eq!(stop.0, 0.0, "zero over-polish leaves the stop layer intact");
    }

    #[test]
    fn a_bad_recipe_is_refused_not_clamped() {
        let m = model();
        let (p, v, t) = (KPa::new(20.0), MetresPerSecond::new(0.5), Minutes::new(1.0));
        assert!(m.removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, p, v, t, -0.1).is_err(), "negative over-polish");
        assert!(m.removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, KPa::new(0.0), v, t, 0.3).is_err(), "zero pressure");
        assert!(m.removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, p, MetresPerSecond::new(0.0), t, 0.3).is_err(), "zero velocity");
        assert!(m.removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, p, v, Minutes::new(0.0), 0.3).is_err(), "zero time");
        assert!(m.removal_on_stop(CmpMaterial::Metal, CmpMaterial::Dielectric, p, v, t, f64::NAN).is_err(), "not a number");
    }
}
