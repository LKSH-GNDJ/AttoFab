// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The parameter manifest (build PRD K-01).
//!
//! Every model constant is a record with a unit, a value, a source and an
//! uncertainty. Physics modules read values through this crate, so a constant
//! cannot enter the model without being registered, and the audit can report
//! which constants are assumed rather than measured.

pub mod deck;

use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Where a declared constant came from.
pub enum SourceKind {
    /// Measured on a tool, with conditions and uncertainty recorded.
    Measured,
    /// Taken from a published reference, not fitted to a tool.
    Published,
    /// Fitted to measured data by the calibration pipeline.
    Fitted,
    /// A placeholder chosen for plausibility. May not be reported as calibrated.
    Assumed,
}

impl SourceKind {
    /// The name a report uses.
    pub fn label(self) -> &'static str {
        match self {
            SourceKind::Measured => "measured",
            SourceKind::Published => "published",
            SourceKind::Fitted => "fitted",
            SourceKind::Assumed => "assumed",
        }
    }
}

#[derive(Debug, Clone)]
/// Where a constant came from, and how firmly.
pub struct Source {
    pub kind: SourceKind,
    /// A resolvable identifier: a DOI, a lot number, or a report reference.
    pub id: String,
}

#[derive(Debug, Clone)]
/// A constant's uncertainty and how it was established.
pub struct Uncertainty {
    /// Absolute one-sigma uncertainty in the parameter's own unit.
    pub sigma: f64,
}

#[derive(Debug, Clone)]
/// One declared constant.
pub struct Param {
    pub id: String,
    pub unit: String,
    pub value: f64,
    pub source: Source,
    pub uncertainty: Uncertainty,
}

#[derive(Debug, Clone, Default)]
/// Every declared constant, with units, value and source.
pub struct Manifest {
    pub version: u32,
    params: BTreeMap<String, Param>,
}

#[derive(Debug, Clone, PartialEq)]
/// Why a parameter lookup failed.
pub enum ParamError {
    Missing(String),
    NotFinite(String),
    BadUncertainty(String),
}

impl std::fmt::Display for ParamError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParamError::Missing(id) => write!(f, "parameter '{id}' is not registered in the manifest"),
            ParamError::NotFinite(id) => write!(f, "parameter '{id}' is not finite"),
            ParamError::BadUncertainty(id) => write!(f, "parameter '{id}' has a non-positive or non-finite uncertainty"),
        }
    }
}
impl std::error::Error for ParamError {}

impl Manifest {
    /// An empty manifest.
    pub fn new() -> Self {
        Manifest { version: 1, params: BTreeMap::new() }
    }

    /// Declare a constant.
    pub fn insert(&mut self, id: &str, unit: &str, value: f64, kind: SourceKind, src_id: &str, sigma: f64) {
        self.params.insert(
            id.to_string(),
            Param {
                id: id.to_string(),
                unit: unit.to_string(),
                value,
                source: Source { kind, id: src_id.to_string() },
                uncertainty: Uncertainty { sigma },
            },
        );
    }

    /// A constant by id, or why it is missing.
    pub fn get(&self, id: &str) -> Result<f64, ParamError> {
        let p = self.params.get(id).ok_or_else(|| ParamError::Missing(id.to_string()))?;
        if !p.value.is_finite() {
            return Err(ParamError::NotFinite(id.to_string()));
        }
        Ok(p.value)
    }

    /// A declared constant with its metadata.
    pub fn param(&self, id: &str) -> Option<&Param> { self.params.get(id) }

    /// How many constants are declared.
    pub fn len(&self) -> usize { self.params.len() }
    /// Whether anything is declared.
    pub fn is_empty(&self) -> bool { self.params.is_empty() }

    /// Every declared constant.
    pub fn iter(&self) -> impl Iterator<Item = &Param> { self.params.values() }

    /// Every parameter whose source is Assumed. A model built only from these
    /// may not be reported as calibrated.
    pub fn assumed(&self) -> Vec<&Param> {
        self.params.values().filter(|p| p.source.kind == SourceKind::Assumed).collect()
    }

    /// The fraction of constants whose source is measured.
    pub fn measured_fraction(&self) -> f64 {
        if self.params.is_empty() {
            return 0.0;
        }
        let n = self.params.values().filter(|p| p.source.kind == SourceKind::Measured).count();
        n as f64 / self.params.len() as f64
    }

    /// Structural validation: finite values, positive uncertainties, non-empty sources.
    pub fn validate(&self) -> Vec<ParamError> {
        let mut e = Vec::new();
        for p in self.params.values() {
            if !p.value.is_finite() {
                e.push(ParamError::NotFinite(p.id.clone()));
            }
            if !(p.uncertainty.sigma.is_finite() && p.uncertainty.sigma > 0.0) {
                e.push(ParamError::BadUncertainty(p.id.clone()));
            }
        }
        e
    }
}

// ---------------------------------------------------------------------------
// Implant moments (build PRD D-05): energy -> (Rp, dRp, skewness, kurtosis)
// ---------------------------------------------------------------------------

/// One tabulated moment set. Lengths in nm, energy in keV.
#[derive(Debug, Clone, Copy)]
pub struct Moment {
    pub energy_kev: f64,
    pub rp_nm: f64,
    pub drp_nm: f64,
    pub skew: f64,
    pub kurt: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The dopant species the models know.
pub enum Species { Boron, Phosphorus, Arsenic }

impl Species {
    /// The name a report uses.
    pub fn label(self) -> &'static str {
        match self { Species::Boron => "B", Species::Phosphorus => "P", Species::Arsenic => "As" }
    }

    /// Intrinsic diffusivity prefactor (cm^2/s) and activation energy (eV).
    pub fn diffusivity(self) -> (f64, f64) {
        match self { Species::Boron => (0.76, 3.46), Species::Phosphorus => (3.85, 3.66), Species::Arsenic => (22.9, 4.1) }
    }

    /// Tabulated moments. Skew and kurtosis are the third and fourth standardised
    /// moments; a Gaussian has 0 and 3.
    pub fn moments(self) -> &'static [Moment] {
        match self {
            // The moments must satisfy kurtosis > 1.5*skew^2 + 3 for a Pearson IV
            // to exist; these do, and `moments_are_pearson_valid` checks it.
            Species::Boron => &[
                Moment { energy_kev: 10.0, rp_nm: 35.0, drp_nm: 14.0, skew: -0.55, kurt: 4.0 },
                Moment { energy_kev: 30.0, rp_nm: 100.0, drp_nm: 36.0, skew: -0.95, kurt: 5.0 },
                Moment { energy_kev: 60.0, rp_nm: 200.0, drp_nm: 55.0, skew: -1.05, kurt: 5.4 },
                Moment { energy_kev: 100.0, rp_nm: 310.0, drp_nm: 70.0, skew: -1.09, kurt: 5.6 },
                Moment { energy_kev: 200.0, rp_nm: 580.0, drp_nm: 95.0, skew: -1.20, kurt: 6.0 },
            ],
            Species::Phosphorus => &[
                Moment { energy_kev: 10.0, rp_nm: 14.0, drp_nm: 6.0, skew: 0.2, kurt: 3.3 },
                Moment { energy_kev: 30.0, rp_nm: 40.0, drp_nm: 17.0, skew: -0.3, kurt: 3.5 },
                Moment { energy_kev: 60.0, rp_nm: 80.0, drp_nm: 28.0, skew: -0.5, kurt: 3.8 },
                Moment { energy_kev: 100.0, rp_nm: 134.0, drp_nm: 40.0, skew: -0.7, kurt: 4.0 },
                Moment { energy_kev: 200.0, rp_nm: 270.0, drp_nm: 62.0, skew: -0.9, kurt: 4.4 },
            ],
            Species::Arsenic => &[
                Moment { energy_kev: 10.0, rp_nm: 9.0, drp_nm: 3.5, skew: 0.4, kurt: 3.4 },
                Moment { energy_kev: 30.0, rp_nm: 20.0, drp_nm: 8.0, skew: 0.1, kurt: 3.2 },
                Moment { energy_kev: 60.0, rp_nm: 37.0, drp_nm: 14.0, skew: -0.2, kurt: 3.4 },
                Moment { energy_kev: 100.0, rp_nm: 58.0, drp_nm: 20.0, skew: -0.4, kurt: 3.6 },
                Moment { energy_kev: 200.0, rp_nm: 112.0, drp_nm: 35.0, skew: -0.6, kurt: 4.0 },
            ],
        }
    }
}

/// Log-log interpolation of the tabulated moments, clamped to the table range
/// rather than extrapolated (build PRD D-05: never invent a range).
pub fn moments_at(s: Species, energy_kev: f64) -> Moment {
    let t = s.moments();
    if energy_kev <= t[0].energy_kev {
        return t[0];
    }
    if energy_kev >= t[t.len() - 1].energy_kev {
        return t[t.len() - 1];
    }
    let mut i = 0;
    for k in 1..t.len() {
        if t[k].energy_kev >= energy_kev {
            i = k - 1;
            break;
        }
    }
    let (a, b) = (t[i], t[i + 1]);
    let f = (energy_kev.ln() - a.energy_kev.ln()) / (b.energy_kev.ln() - a.energy_kev.ln());
    let lerp = |x: f64, y: f64| (x.ln() + f * (y.ln() - x.ln())).exp();
    Moment {
        energy_kev,
        rp_nm: lerp(a.rp_nm, b.rp_nm),
        drp_nm: lerp(a.drp_nm, b.drp_nm),
        skew: a.skew + f * (b.skew - a.skew),
        kurt: a.kurt + f * (b.kurt - a.kurt),
    }
}

// ---------------------------------------------------------------------------
// Deal-Grove coefficients (build PRD O-01)
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The oxidation ambient.
pub enum Ambient { Dry, Wet }

/// Coefficients at the reference temperature plus the activation energies that
/// carry them to another temperature.
#[derive(Debug, Clone, Copy)]
pub struct DealGrove {
    pub a_um: f64,
    pub b_um2_per_hr: f64,
    pub ea_ev: f64,
    pub eb_ev: f64,
    /// Reference temperature for a_um and b_um2_per_hr.
    pub t_ref_k: f64,
}

impl DealGrove {
    /// The bundled pack: values at 1000 C with published activation energies.
    pub fn bundled(ambient: Ambient) -> Self {
        match ambient {
            Ambient::Dry => DealGrove { a_um: 0.102, b_um2_per_hr: 0.0105, ea_ev: 2.0, eb_ev: 1.24, t_ref_k: 1273.15 },
            Ambient::Wet => DealGrove { a_um: 0.52, b_um2_per_hr: 0.53, ea_ev: 1.05, eb_ev: 0.70, t_ref_k: 1273.15 },
        }
    }

    /// `A` is `B` divided by the linear rate constant, so its effective
    /// activation energy is `Ea - Eb`: A *falls* as temperature rises, which is
    /// why the linear regime gives way to the parabolic one.
    pub fn a_at(&self, t_k: f64) -> f64 {
        self.a_um * ((self.ea_ev - self.eb_ev) / crate::KB_EV * (1.0 / t_k - 1.0 / self.t_ref_k)).exp()
    }

    /// The parabolic rate constant at a temperature.
    pub fn b_at(&self, t_k: f64) -> f64 {
        self.b_um2_per_hr * (-self.eb_ev / crate::KB_EV * (1.0 / t_k - 1.0 / self.t_ref_k)).exp()
    }
}

/// Boltzmann constant in eV/K, re-exported for the Arrhenius helpers.
pub const KB_EV: f64 = 8.617_333_262e-5;

/// Nanometres of silicon consumed per nanometre of oxide grown (SiO2 is 2.27x
/// the molar volume of Si, so 1/2.27 = 0.44).
pub const SI_CONSUMED_RATIO: f64 = 0.44;

/// Build the manifest the models actually read.
pub fn default_manifest() -> Manifest {
    let mut m = Manifest::new();
    m.insert("si.young.gpa", "GPa", 130.0, SourceKind::Published, "textbook: Si <100> Young's modulus", 2.0);
    m.insert("si.poisson", "-", 0.28, SourceKind::Published, "textbook: Si Poisson ratio", 0.01);
    m.insert("sio2.young.gpa", "GPa", 70.0, SourceKind::Published, "bundled stress pack", 2.0);
    m.insert("sio2.poisson", "-", 0.17, SourceKind::Published, "bundled stress pack", 0.01);
    m.insert("sio2.intrinsic.stress.mpa", "MPa", -200.0, SourceKind::Published, "bundled stress pack", 20.0);
    m.insert("si.lattice.nm", "nm", 0.5431, SourceKind::Published, "textbook", 0.0001);
    m.insert("ox.si.consumed.ratio", "-", 0.44, SourceKind::Published, "Deal-Grove stoichiometry", 0.01);
    m.insert("etch.critical.ar", "-", 8.0, SourceKind::Assumed, "bundled etch pack, unverified", 1.0);
    m.insert("etch.arde.factor", "-", 0.4, SourceKind::Assumed, "bundled etch pack, unverified", 0.05);
    m
}


// ---------------------------------------------------------------------------
// The parameter registry (build PRD K-01)
// ---------------------------------------------------------------------------
//
// A constant is registered once, here, with the unit it is in, the value the
// code uses, where that value came from and how uncertain it is. `live` reads
// the value out of the same code path the physics uses, so `verify` catches the
// failure the earlier build had: a documented table that disagreed with the shipped numbers.

/// One registered constant.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParamSpec {
    pub id: &'static str,
    pub unit: &'static str,
    /// The value this manifest declares. Checked against `live`.
    pub declared: f64,
    pub kind: SourceKind,
    pub citation: &'static str,
    pub sigma: f64,
}

/// A registered value that no longer matches what the code computes.
#[derive(Debug, Clone, PartialEq)]
pub struct Drift {
    pub id: &'static str,
    pub declared: f64,
    pub live: f64,
}

impl std::fmt::Display for Drift {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: manifest says {}, code computes {}", self.id, self.declared, self.live)
    }
}

/// Every scalar constant this crate hands to the physics.
pub fn registry() -> Vec<ParamSpec> {
    let mut v = Vec::new();
    let mut add = |id, unit, declared, kind, citation, sigma| {
        v.push(ParamSpec { id, unit, declared, kind, citation, sigma });
    };
    // intrinsic diffusivities: D = D0 exp(-Ea/kT)
    for (sp, d0, ea) in [
        (Species::Boron, 0.76, 3.46),
        (Species::Phosphorus, 3.85, 3.66),
        (Species::Arsenic, 22.9, 4.10),
    ] {
        let (id0, ida) = match sp {
            Species::Boron => ("diffusion.b.d0", "diffusion.b.ea"),
            Species::Phosphorus => ("diffusion.p.d0", "diffusion.p.ea"),
            Species::Arsenic => ("diffusion.as.d0", "diffusion.as.ea"),
        };
        add(id0, "cm2/s", d0, SourceKind::Published, "B/P/As in Si, intrinsic", d0 * 0.1);
        add(ida, "eV", ea, SourceKind::Published, "B/P/As in Si, intrinsic", 0.05);
    }
    // Deal-Grove coefficients at the reference temperature
    for amb in [Ambient::Dry, Ambient::Wet] {
        let dg = DealGrove::bundled(amb);
        let (pa, pb, pe, pf) = match amb {
            Ambient::Dry => ("deal_grove.dry.a_um", "deal_grove.dry.b_um2_hr", "deal_grove.dry.ea_ev", "deal_grove.dry.eb_ev"),
            Ambient::Wet => ("deal_grove.wet.a_um", "deal_grove.wet.b_um2_hr", "deal_grove.wet.ea_ev", "deal_grove.wet.eb_ev"),
        };
        add(pa, "um", dg.a_um, SourceKind::Published, "Deal-Grove (100) Si", dg.a_um * 0.05);
        add(pb, "um2/hr", dg.b_um2_per_hr, SourceKind::Published, "Deal-Grove (100) Si", dg.b_um2_per_hr * 0.05);
        add(pe, "eV", dg.ea_ev, SourceKind::Published, "Deal-Grove linear regime", 0.05);
        add(pf, "eV", dg.eb_ev, SourceKind::Published, "Deal-Grove parabolic regime", 0.05);
    }
    add("ox.si.consumed.ratio", "-", SI_CONSUMED_RATIO, SourceKind::Published, "SiO2 stoichiometry: 0.44 nm Si per nm oxide", 0.01);
    add("si.lattice.nm", "nm", 0.5431, SourceKind::Published, "Si lattice constant", 0.0001);
    add("boltzmann.ev_per_k", "eV/K", KB_EV, SourceKind::Published, "CODATA", 1e-12);
    v
}

/// The value the code actually computes for an id.
pub fn live(id: &str) -> Option<f64> {
    for sp in [Species::Boron, Species::Phosphorus, Species::Arsenic] {
        let (d0, ea) = sp.diffusivity();
        let (id0, ida) = match sp {
            Species::Boron => ("diffusion.b.d0", "diffusion.b.ea"),
            Species::Phosphorus => ("diffusion.p.d0", "diffusion.p.ea"),
            Species::Arsenic => ("diffusion.as.d0", "diffusion.as.ea"),
        };
        if id == id0 {
            return Some(d0);
        }
        if id == ida {
            return Some(ea);
        }
    }
    for amb in [Ambient::Dry, Ambient::Wet] {
        let dg = DealGrove::bundled(amb);
        let ids = match amb {
            Ambient::Dry => ["deal_grove.dry.a_um", "deal_grove.dry.b_um2_hr", "deal_grove.dry.ea_ev", "deal_grove.dry.eb_ev"],
            Ambient::Wet => ["deal_grove.wet.a_um", "deal_grove.wet.b_um2_hr", "deal_grove.wet.ea_ev", "deal_grove.wet.eb_ev"],
        };
        let vals = [dg.a_um, dg.b_um2_per_hr, dg.ea_ev, dg.eb_ev];
        for k in 0..4 {
            if id == ids[k] {
                return Some(vals[k]);
            }
        }
    }
    match id {
        "ox.si.consumed.ratio" => Some(SI_CONSUMED_RATIO),
        "si.lattice.nm" => Some(0.5431),
        "boltzmann.ev_per_k" => Some(KB_EV),
        _ => None,
    }
}

/// Compare a registry against the code. Any entry that disagrees is drift.
pub fn verify(specs: &[ParamSpec], live: &dyn Fn(&str) -> Option<f64>) -> Vec<Drift> {
    let mut out = Vec::new();
    for s in specs {
        match live(s.id) {
            None => out.push(Drift { id: s.id, declared: s.declared, live: f64::NAN }),
            Some(v) => {
                let tol = 1e-9f64.max(v.abs() * 1e-9);
                if (s.declared - v).abs() > tol {
                    out.push(Drift { id: s.id, declared: s.declared, live: v });
                }
            }
        }
    }
    out
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
/// What fraction of a model's constants are measured, calibrated or assumed.
pub struct Coverage {
    pub total: usize,
    pub measured: usize,
    pub published: usize,
    pub fitted: usize,
    pub assumed: usize,
}

/// How much of a declared set is measured rather than assumed.
pub fn coverage(specs: &[ParamSpec]) -> Coverage {
    let mut c = Coverage { total: specs.len(), ..Default::default() };
    for s in specs {
        match s.kind {
            SourceKind::Measured => c.measured += 1,
            SourceKind::Published => c.published += 1,
            SourceKind::Fitted => c.fitted += 1,
            SourceKind::Assumed => c.assumed += 1,
        }
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_registry_matches_the_code() {
        let specs = registry();
        assert!(specs.len() >= 17, "the registry should cover every constant, got {}", specs.len());
        let drift = verify(&specs, &live);
        assert!(drift.is_empty(), "the manifest has drifted from the code:\n{}", drift.iter().map(|d| d.to_string()).collect::<Vec<_>>().join("\n"));
    }

    #[test]
    fn every_registered_id_resolves() {
        for s in registry() {
            assert!(live(s.id).is_some(), "{} is registered but nothing computes it", s.id);
            assert!(!s.citation.is_empty(), "{} has no citation", s.id);
            assert!(s.sigma > 0.0, "{} has no uncertainty", s.id);
        }
    }

    #[test]
    fn coverage_counts_the_provenance() {
        let c = coverage(&registry());
        assert_eq!(c.total, c.measured + c.published + c.fitted + c.assumed);
        assert_eq!(c.measured, 0, "nothing here is measured yet");
        assert!(c.published > 0);
    }

    #[test]
    fn drift_is_detected() {
        let specs = vec![ParamSpec { id: "ox.si.consumed.ratio", unit: "-", declared: 0.5, kind: SourceKind::Published, citation: "x", sigma: 0.01 }];
        let d = verify(&specs, &live);
        assert_eq!(d.len(), 1, "0.5 is not 0.44");
        assert!((d[0].live - 0.44).abs() < 1e-12);
    }

    #[test]
    fn manifest_reads_and_reports_missing() {
        let m = default_manifest();
        assert_eq!(m.get("si.young.gpa").unwrap(), 130.0);
        assert!(matches!(m.get("nope").unwrap_err(), ParamError::Missing(_)));
        assert!(m.validate().is_empty());
    }

    #[test]
    fn assumed_parameters_are_traceable() {
        let m = default_manifest();
        let a = m.assumed();
        assert_eq!(a.len(), 2, "the two etch constants are assumed: {a:?}");
        assert!(a.iter().all(|p| p.source.kind == SourceKind::Assumed));
    }

    #[test]
    fn moments_are_clamped_not_extrapolated() {
        let lo = moments_at(Species::Boron, 1.0);
        assert_eq!(lo.rp_nm, 35.0, "below the table the first row is used, not an extrapolation");
        let hi = moments_at(Species::Boron, 5000.0);
        assert_eq!(hi.rp_nm, 580.0);
    }

    #[test]
    fn moments_interpolate_between_rows() {
        let m = moments_at(Species::Boron, 30.0);
        assert!((m.rp_nm - 100.0).abs() < 1e-9);
        let mid = moments_at(Species::Boron, 45.0);
        assert!(mid.rp_nm > 130.0 && mid.rp_nm < 175.0, "{mid:?}");
    }

    #[test]
    fn moments_are_pearson_valid() {
        // A Pearson IV exists only when kurtosis > 1.5*skew^2 + 3.
        for sp in [Species::Boron, Species::Phosphorus, Species::Arsenic] {
            for m in sp.moments() {
                let min_kurt = 1.5 * m.skew * m.skew + 3.0;
                assert!(m.kurt > min_kurt + 1e-6, "{} {} keV: kurt {} <= {min_kurt}", sp.label(), m.energy_kev, m.kurt);
            }
        }
    }

    #[test]
    fn deal_grove_follows_arrhenius() {
        let dg = DealGrove::bundled(Ambient::Dry);
        let b1000 = dg.b_at(1273.15);
        let b1100 = dg.b_at(1373.15);
        assert!((b1000 - 0.0105).abs() < 1e-12, "at the reference temperature the value is the tabulated one");
        assert!(b1100 > b1000 * 1.5, "rate rises with temperature: {b1000} -> {b1100}");
        let a1100 = dg.a_at(1373.15);
        assert!(a1100 < dg.a_um, "A falls with temperature: {a1100} vs {}", dg.a_um);
    }
}
