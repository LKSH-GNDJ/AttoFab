// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Plasma etch rate model (build PRD E-01, E-04).
//!
//! Vertical rate is a function of RF power and chamber pressure, interpolated
//! log-log from a calibration table. Aspect-ratio dependent etching reduces the
//! rate in deep features. Both are the two inputs the earlier build accepted and ignored.

use atto_units::{MTorr, Nm, Watts};
use atto_geom::column::EtchRates;
use atto_geom::Material;

/// One calibration point.
#[derive(Debug, Clone, Copy)]
pub struct PlasmaPoint {
    pub pressure_mtorr: f64,
    pub power_w: f64,
    pub vert_nm_min: f64,
    pub lat_nm_min: f64,
}

#[derive(Debug, Clone)]
/// The etch model: rates as a function of the operating point.
pub struct EtchModel {
    /// Chemistry name, for the report.
    pub chemistry: String,
    pub points: Vec<PlasmaPoint>,
    /// Selectivity of the target material over each mask material.
    pub selectivity: Vec<(Material, f64)>,
    /// Aspect ratio above which the rate starts to fall.
    pub critical_aspect_ratio: f64,
    /// Fraction of the rate retained at twice the critical aspect ratio.
    pub arde_factor: f64,
    /// Loading: how strongly an open area depletes the reactant. 0 is no effect.
    pub loading_coeff: f64,
    /// Open fraction the tabulated rates were measured at. The loading factor is 1.0 here.
    pub reference_open_fraction: f64,
    /// Residual angular spread of the ion flux at zero pressure, degrees.
    pub sigma_zero_deg: f64,
    /// Additional angular spread per mTorr of pressure, degrees.
    pub sigma_per_mtorr: f64,
}

impl EtchModel {
    /// The bundled CF4/CHF3 oxide pack. Unverified provenance, as the audit found.
    pub fn bundled_oxide() -> Self {
        EtchModel {
            chemistry: "CF4/CHF3 (SiO2)".into(),
            points: vec![
                PlasmaPoint { pressure_mtorr: 10.0, power_w: 100.0, vert_nm_min: 45.0, lat_nm_min: 8.0 },
                PlasmaPoint { pressure_mtorr: 20.0, power_w: 200.0, vert_nm_min: 100.0, lat_nm_min: 12.0 },
                PlasmaPoint { pressure_mtorr: 30.0, power_w: 300.0, vert_nm_min: 160.0, lat_nm_min: 18.0 },
            ],
            selectivity: vec![(Material::Silicon, 5.0), (Material::Photoresist, 10.0), (Material::SiliconNitride, 1.0)],
            critical_aspect_ratio: 8.0,
            arde_factor: 0.4,
            loading_coeff: 0.5,
            reference_open_fraction: 0.0,
            sigma_zero_deg: 0.5,
            sigma_per_mtorr: 0.13,
        }
    }

    /// Rate in nm/min as `f(power) * g(pressure)`.
    ///
    /// The bundled table has power and pressure collinear (100 W/10 mT,
    /// 200 W/20 mT, 300 W/30 mT), so the two effects cannot be separated by the
    /// data alone. The pressure dependence is therefore declared: `g` is a
    /// square-root law normalised at the median table pressure, and `f` is the
    /// log-log interpolation of the table after dividing by `g`. That
    /// reproduces every calibration point exactly, is monotone in both inputs,
    /// and is honest about which part is measured and which is assumed. A
    /// non-collinear table would replace this with a bilinear grid.
    pub fn rate(&self, power: Watts, pressure: MTorr) -> (f64, f64) {
        // A-03: the operating point is dimensioned. A caller that passes mTorr
        // where watts belong now fails to compile instead of shifting the curve.
        let (power_w, pressure_mtorr) = (power.value(), pressure.value());
        if self.points.is_empty() {
            return (0.0, 0.0);
        }
        let mut pts = self.points.clone();
        pts.sort_by(|a, b| a.power_w.partial_cmp(&b.power_w).unwrap_or(std::cmp::Ordering::Equal));
        let p_ref = pts[pts.len() / 2].pressure_mtorr;
        let g = |p: f64| (p.max(1e-6) / p_ref).sqrt();
        let scaled: Vec<(f64, f64, f64)> = pts.iter().map(|q| (q.power_w, q.vert_nm_min / g(q.pressure_mtorr), q.lat_nm_min / g(q.pressure_mtorr))).collect();
        let interp = |vals: &dyn Fn(usize) -> f64, p: f64| -> f64 {
            let p = p.max(1e-6);
            if p <= scaled[0].0 { return vals(0); }
            if p >= scaled[scaled.len() - 1].0 { return vals(scaled.len() - 1); }
            let mut i = 0;
            for k in 1..scaled.len() {
                if scaled[k].0 >= p { i = k - 1; break; }
            }
            let (a, b) = (scaled[i].0, scaled[i + 1].0);
            let f = (p.ln() - a.ln()) / (b.ln() - a.ln());
            let (va, vb) = (vals(i), vals(i + 1));
            (va.ln() + f * (vb.ln() - va.ln())).exp()
        };
        let v = interp(&|i| scaled[i].1, power_w) * g(pressure_mtorr);
        let l = interp(&|i| scaled[i].2, power_w) * g(pressure_mtorr);
        (v, l)
    }

    /// Rate multiplier from reactant loading (build PRD E-05).
    ///
    /// A larger exposed area consumes the reactant faster, so the rate falls:
    /// `1 / (1 + k * open_fraction)`. `k` is an assumed constant, registered in
    /// the manifest, because the bundled pack does not measure it.
    pub fn loading_factor(&self, open_fraction: f64) -> f64 {
        let f = open_fraction.clamp(0.0, 1.0);
        if !(f.is_finite()) || !self.reference_open_fraction.is_finite() {
            return 1.0;
        }
        // `1 / (1 + c f)` penalises every open area including the one the rates came from:
        // at f = 0 it returns the full rate, so the table is effectively the rate at zero
        // open area, which nothing is measured at. Normalising against a reference makes the
        // tabulated rate apply where it was actually measured.
        let r = self.reference_open_fraction.clamp(0.0, 1.0);
        (1.0 + self.loading_coeff * r) / (1.0 + self.loading_coeff * f)
    }

    /// Rates after both loading and ARDE, nm/min.
    pub fn rates_loaded(&self, power: Watts, pressure: MTorr, open_fraction: f64, aspect_ratio: f64) -> (f64, f64) {
        let (v, l) = self.rate(power, pressure);
        let m = self.loading_factor(open_fraction) * self.arde_multiplier(aspect_ratio);
        (v * m, l * m)
    }

    /// Lateral etch fraction from the ion angular distribution (build PRD E-03).
    ///
    /// Ions cross the sheath with a Gaussian angular spread that widens with
    /// pressure. Only those arriving outside the sidewall acceptance angle of
    /// the feature etch laterally, so the lateral rate is the vertical rate
    /// times the tail probability:
    ///
    /// `tan(theta_max) = width / (2 * depth)`, `P(|theta| > theta_max)` for
    /// `theta ~ N(0, sigma(p))`.
    ///
    /// A *wider* opening has a larger acceptance angle and therefore a smaller
    /// lateral share; a narrow one lets more of the angular tail reach the
    /// sidewall. Sidewall polymer and mask shadowing are not modelled — those
    /// suppress lateral etch in a deep feature, and the depth dependence here
    /// comes from ARDE instead.
    pub fn lateral_fraction(&self, pressure: MTorr, width: Nm, depth_nm: f64) -> f64 {
        let width_nm = width.value();
        let pressure_mtorr = pressure.value();
        if !(width_nm > 0.0) || !(depth_nm > 0.0) || !pressure_mtorr.is_finite() {
            return 0.0;
        }
        let sigma_deg = self.angular_sigma_deg(pressure_mtorr);
        if sigma_deg <= 0.0 {
            return 0.0;
        }
        let theta_max = (width_nm / (2.0 * depth_nm)).atan().to_degrees();
        // two-sided tail of a normal with zero mean
        2.0 * 0.5 * erfc(theta_max / (sigma_deg * std::f64::consts::SQRT_2))
    }

    /// Angular spread of the ion flux, degrees, from sheath scattering.
    pub fn angular_sigma_deg(&self, pressure_mtorr: f64) -> f64 {
        if !(pressure_mtorr > 0.0) || !pressure_mtorr.is_finite() {
            return 0.0;
        }
        self.sigma_zero_deg + self.sigma_per_mtorr * pressure_mtorr
    }

    /// Rates with the lateral component derived from the ion angles rather than
    /// tabulated, nm/min.
    pub fn rates_anisotropic(&self, power: Watts, pressure: MTorr, width: Nm, depth_nm: f64) -> (f64, f64) {
        let (v, _) = self.rate(power, pressure);
        (v, v * self.lateral_fraction(pressure, width, depth_nm))
    }

    /// Rate multiplier for aspect-ratio dependent etching.
    pub fn arde_multiplier(&self, aspect_ratio: f64) -> f64 {
        if !aspect_ratio.is_finite() || aspect_ratio <= self.critical_aspect_ratio {
            return 1.0;
        }
        let excess = (aspect_ratio / self.critical_aspect_ratio - 1.0).min(1.0);
        1.0 - (1.0 - self.arde_factor) * excess
    }

    /// Per-material rates for the geometry operator, with the mask selectivity.
    pub fn rates(&self, power: Watts, pressure: MTorr) -> EtchRates {
        // the operating point goes straight to `rate`; it is not read here
        let (vert, _lat) = self.rate(power, pressure);
        let sel = |m: Material| self.selectivity.iter().find(|(mm, _)| *mm == m).map(|(_, s)| *s).unwrap_or(1.0);
        EtchRates {
            si: vert / sel(Material::Silicon).max(1e-6),
            sio2: vert,
            si3n4: vert / sel(Material::SiliconNitride).max(1e-6),
            resist: vert / sel(Material::Photoresist).max(1e-6),
            metal: vert / 20.0,
        }
    }
}

/// Complementary error function, Abramowitz & Stegun 7.1.26.
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * z);
    let poly = t * (0.254829592
        + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let e = poly * (-z * z).exp();
    if x >= 0.0 { e } else { 2.0 - e }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reproduces_the_calibration_points() {
        let m = EtchModel::bundled_oxide();
        for p in m.points.iter() {
            let (v, l) = m.rate(Watts::new(p.power_w), MTorr::new(p.pressure_mtorr));
            assert!((v - p.vert_nm_min).abs() < 1.0, "at {} W {} mT: {v:.1} vs {:.1}", p.power_w, p.pressure_mtorr, p.vert_nm_min);
            assert!((l - p.lat_nm_min).abs() < 1.0);
        }
    }

    #[test]
    fn rate_rises_monotonically_with_power_at_fixed_pressure() {
        let m = EtchModel::bundled_oxide();
        let mut prev = 0.0;
        for w in [50.0, 100.0, 150.0, 200.0, 300.0, 400.0] {
            let (v, _) = m.rate(Watts::new(w), MTorr::new(20.0));
            assert!(v >= prev - 1e-9, "power {w}: {v:.1} < {prev:.1}");
            prev = v;
        }
    }

    #[test]
    fn rate_rises_with_pressure_at_fixed_power() {
        let m = EtchModel::bundled_oxide();
        let mut prev = 0.0;
        for p in [5.0, 10.0, 20.0, 30.0, 40.0] {
            let (v, _) = m.rate(Watts::new(200.0), MTorr::new(p));
            assert!(v >= prev - 1e-9, "pressure {p}: {v:.1} < {prev:.1}");
            prev = v;
        }
    }

    #[test]
    fn arde_reduces_deep_features() {
        let m = EtchModel::bundled_oxide();
        assert_eq!(m.arde_multiplier(1.0), 1.0);
        assert_eq!(m.arde_multiplier(8.0), 1.0, "at the critical ratio there is no penalty");
        let deep = m.arde_multiplier(16.0);
        assert!((deep - 0.4).abs() < 1e-9, "at twice the critical ratio the factor applies: {deep}");
        assert!(m.arde_multiplier(40.0) >= 0.4 - 1e-9, "the penalty saturates");
    }

    #[test]
    fn loading_reduces_the_rate_as_more_area_is_exposed() {
        let m = EtchModel::bundled_oxide();
        let a = m.loading_factor(0.0);
        let b = m.loading_factor(0.25);
        let c = m.loading_factor(1.0);
        assert!((a - 1.0).abs() < 1e-12, "no open area, no depletion");
        assert!(b < a && c < b, "the rate must fall: {a:.3} {b:.3} {c:.3}");
        // 1/(1+0.5) at full open area
        assert!((c - 1.0 / 1.5).abs() < 1e-12);
    }

    #[test]
    fn loading_composes_with_arde() {
        let m = EtchModel::bundled_oxide();
        let (v0, _) = m.rate(Watts::new(200.0), MTorr::new(20.0));
        let (v1, _) = m.rates_loaded(Watts::new(200.0), MTorr::new(20.0), 0.0, 0.0);
        let (v2, _) = m.rates_loaded(Watts::new(200.0), MTorr::new(20.0), 0.5, 0.0);
        assert!((v1 - v0).abs() < 1e-9, "nothing open and no ARDE leaves the rate alone");
        assert!(v2 < v1, "loading reduces it");
    }

    #[test]
    fn the_rate_applies_at_the_reference_open_fraction() {
        // a zero reference reduces to the old formula exactly, so the default changes nothing;
        // a tool that knows the open fraction its rates came from sets it and gets the table back
        let mut m = EtchModel::bundled_oxide();
        m.reference_open_fraction = 0.5;
        let r = m.reference_open_fraction;
        assert!((m.loading_factor(r) - 1.0).abs() < 1e-12, "the tabulated rate applies where it was measured: {}", m.loading_factor(r));
        assert!(m.loading_factor(1.0) < 1.0, "a wider open area loads more");
        assert!(m.loading_factor(0.0) > 1.0, "a narrower one loads less");
        let mut d = EtchModel::bundled_oxide();
        d.reference_open_fraction = 0.0;
        assert!((d.loading_factor(0.0) - 1.0).abs() < 1e-12, "no reference keeps the old behaviour exactly");
    }

    #[test]
    fn a_zero_loading_coefficient_changes_nothing() {
        let mut m = EtchModel::bundled_oxide();
        m.loading_coeff = 0.0;
        for f in [0.0, 0.5, 1.0] {
            assert!((m.loading_factor(f) - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn loading_is_clamped_to_physical_areas() {
        let m = EtchModel::bundled_oxide();
        assert!((m.loading_factor(-1.0) - m.loading_factor(0.0)).abs() < 1e-12);
        assert!((m.loading_factor(5.0) - m.loading_factor(1.0)).abs() < 1e-12);
        assert!((m.loading_factor(f64::NAN) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn anisotropy_falls_as_pressure_rises() {
        let m = EtchModel::bundled_oxide();
        // a 100 nm wide, 400 nm deep feature: theta_max ~ 7 degrees
        let a = m.lateral_fraction(MTorr::new(5.0), Nm::new(100.0), 400.0);
        let b = m.lateral_fraction(MTorr::new(20.0), Nm::new(100.0), 400.0);
        let c = m.lateral_fraction(MTorr::new(40.0), Nm::new(100.0), 400.0);
        assert!(a < b && b < c, "the lateral share must rise with pressure: {a:.3} {b:.3} {c:.3}");
        assert!(a < 0.01, "5 mTorr is essentially vertical: {a:.2e}");
        assert!(b > 0.018 && b < 0.065, "20 mTorr is a few percent: {b:.3}");
        assert!(c > 0.15, "40 mTorr is clearly isotropic: {c:.3}");
    }

    #[test]
    fn a_wider_opening_etches_less_laterally() {
        // The acceptance angle widens with the opening, so less of the angular
        // tail reaches the sidewall. Shadowing is not modelled; ARDE carries the
        // depth dependence separately.
        let m = EtchModel::bundled_oxide();
        let narrow = m.lateral_fraction(MTorr::new(20.0), Nm::new(50.0), 400.0);
        let wide = m.lateral_fraction(MTorr::new(20.0), Nm::new(400.0), 400.0);
        assert!(wide < narrow, "a wide opening is more directional: {narrow:.4} -> {wide:.4}");
    }

    #[test]
    fn the_lateral_rate_is_derived_not_tabulated() {
        let m = EtchModel::bundled_oxide();
        let (v, l) = m.rates_anisotropic(Watts::new(200.0), MTorr::new(20.0), Nm::new(100.0), 400.0);
        assert!(v > 0.0);
        assert!(l > 0.0 && l < v, "lateral {l:.1} must sit below vertical {v:.1} nm/min");
        // zero pressure spread leaves no lateral component at all
        let (_, l0) = m.rates_anisotropic(Watts::new(200.0), MTorr::new(0.0), Nm::new(100.0), 400.0);
        assert!(l0 < 1e-3, "no scattering, no lateral etch: {l0:.4}");
    }

    #[test]
    fn the_angular_spread_is_linear_in_pressure() {
        let m = EtchModel::bundled_oxide();
        let a = m.angular_sigma_deg(10.0);
        let b = m.angular_sigma_deg(20.0);
        assert!((b - a - 10.0 * m.sigma_per_mtorr).abs() < 1e-9);
        assert!((m.angular_sigma_deg(0.0)).abs() < 1e-12);
    }

    #[test]
    fn selectivity_orders_the_rates() {
        let m = EtchModel::bundled_oxide();
        let r = m.rates(Watts::new(200.0), MTorr::new(20.0));
        assert!(r.sio2 > r.si, "oxide is etched faster than silicon: {} vs {}", r.sio2, r.si);
        assert!(r.sio2 > r.resist, "and faster than resist");
        assert!(r.si3n4 > 0.0);
    }
}
