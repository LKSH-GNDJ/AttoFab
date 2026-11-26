// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Process models. Each is a pure function of its inputs plus the parameter
//! manifest; nothing here reads the clock or touches the filesystem.

pub mod via;
pub mod anneal;
pub mod pair_diffusion;
pub mod implant;
pub mod litho;
pub mod psm;
pub mod mask3d;
pub mod massoud;
pub mod profile;
pub mod stochastic;
pub mod mask_erosion;
pub mod passivation;
pub mod wet_etch;
pub mod oed;
pub mod ald;
pub mod yield_model;
pub mod epitaxy;
pub mod opc;
pub mod oxidation;
pub mod dill;
pub mod flux;
pub mod field3d;
pub mod resist;
pub mod stress;
pub mod stress_coupling;
pub mod deposition;
pub mod etch;
pub mod redeposition;
pub mod flow;
pub mod cmp;
pub mod coupling;

pub use anneal::{AnnealConfig, ThermalHistory, ThermalStep};
pub use implant::{Profile, ProfileError};
pub use oxidation::{DealGroveError, oxidation_nm};
pub use etch::{EtchModel, PlasmaPoint};

// ---------------------------------------------------------------------------
// Parameter registry (build PRD K-01)
// ---------------------------------------------------------------------------

use atto_params::{ParamSpec, SourceKind};

/// Every scalar constant the process models read, with the value the code uses.
pub fn registry() -> Vec<ParamSpec> {
    vec![
        ParamSpec {
            id: "ox.beak.ratio",
            unit: "-",
            declared: oxidation::BEAK_RATIO,
            kind: SourceKind::Assumed,
            citation: "bird's beak length as a multiple of the grown oxide; not measured here",
            sigma: 0.3,
        },
        ParamSpec {
            id: "etch.arde.factor",
            unit: "-",
            declared: etch::EtchModel::bundled_oxide().arde_factor,
            kind: SourceKind::Assumed,
            citation: "bundled etch pack, unverified (audit F-02 class)",
            sigma: 0.1,
        },
        ParamSpec {
            id: "etch.sigma.per.mtorr",
            unit: "deg/mTorr",
            declared: etch::EtchModel::bundled_oxide().sigma_per_mtorr,
            kind: SourceKind::Assumed,
            citation: "ion angular spread per mTorr of sheath scattering; not measured",
            sigma: 0.15,
        },
        ParamSpec {
            id: "etch.loading.coeff",
            unit: "-",
            declared: etch::EtchModel::bundled_oxide().loading_coeff,
            kind: SourceKind::Assumed,
            citation: "reactant loading coefficient; not measured in the bundled pack",
            sigma: 0.2,
        },
        ParamSpec {
            id: "etch.critical.ar",
            unit: "-",
            declared: etch::EtchModel::bundled_oxide().critical_aspect_ratio,
            kind: SourceKind::Assumed,
            citation: "bundled etch pack, unverified",
            sigma: 1.0,
        },
        ParamSpec {
            id: "defect.threshold.ar",
            unit: "-",
            declared: profile::DefectCoefficients::silicon().threshold_ar,
            kind: SourceKind::Assumed,
            citation: "aspect ratio below which no profile defect appears; assumed",
            sigma: 1.0,
        },
        ParamSpec {
            id: "defect.bow.per.ar",
            unit: "-",
            declared: profile::DefectCoefficients::silicon().bow_per_ar,
            kind: SourceKind::Assumed,
            citation: "sidewall bow fraction per unit aspect ratio above threshold; assumed",
            sigma: 0.002,
        },
        ParamSpec {
            id: "defect.notch.per.ar",
            unit: "-",
            declared: profile::DefectCoefficients::silicon().notch_per_ar,
            kind: SourceKind::Assumed,
            citation: "foot notch fraction per unit aspect ratio above threshold; assumed",
            sigma: 0.002,
        },
        ParamSpec {
            id: "defect.microtrench.per.ar",
            unit: "-",
            declared: profile::DefectCoefficients::silicon().microtrench_per_ar,
            kind: SourceKind::Assumed,
            citation: "microtrench fraction per unit aspect ratio above threshold; assumed",
            sigma: 0.001,
        },
        ParamSpec {
            id: "euv.photon.j",
            unit: "J",
            declared: stochastic::EUV_PHOTON_J,
            kind: SourceKind::Published,
            citation: "photon energy at 13.5 nm, hc/lambda",
            sigma: 0.0,
        },
        ParamSpec {
            id: "euv.lwr.response",
            unit: "nm",
            declared: stochastic::LWR_RESPONSE_NM,
            kind: SourceKind::Assumed,
            citation: "resist edge response per unit fractional dose fluctuation; assumed, not fitted",
            sigma: 1.0,
        },
        ParamSpec {
            id: "seg.m.boron",
            unit: "-",
            declared: atto_dopant::segregation::BORON_SEGREGATION_M,
            kind: SourceKind::Assumed,
            citation: "boron segregation coefficient at Si/SiO2; the oxide is the sink, magnitude assumed",
            sigma: 0.1,
        },
        ParamSpec {
            id: "oed.k",
            unit: "1/(nm/min)",
            declared: oed::OED_K,
            kind: SourceKind::Assumed,
            citation: "interstitial supersaturation per unit oxidation rate; chosen to make boron enhancement visible, not fitted",
            sigma: 1.0,
        },
        ParamSpec {
            id: "ald.nm.per.cycle",
            unit: "nm",
            declared: ald::AldProcess::oxide().nm_per_cycle,
            kind: SourceKind::Assumed,
            citation: "saturated growth per ALD cycle for a silica process; typical magnitude, not measured here",
            sigma: 0.02,
        },
        ParamSpec {
            id: "ald.saturation.dose",
            unit: "-",
            declared: ald::AldProcess::oxide().saturation_dose,
            kind: SourceKind::Assumed,
            citation: "precursor dose at which the surface saturates; process units, chosen for a visible knee",
            sigma: 0.2,
        },
        ParamSpec {
            id: "yield.threshold.ev",
            unit: "eV",
            declared: yield_model::YieldParams::silicon().threshold_ev,
            kind: SourceKind::Assumed,
            citation: "sputter threshold energy for Si; assumed",
            sigma: 5.0,
        },
        ParamSpec {
            id: "yield.reference.ev",
            unit: "eV",
            declared: yield_model::YieldParams::silicon().reference_ev,
            kind: SourceKind::Assumed,
            citation: "reference ion energy for the yield normalisation; assumed",
            sigma: 20.0,
        },
        ParamSpec {
            id: "yield.ref",
            unit: "atoms/ion",
            declared: yield_model::YieldParams::silicon().y_ref,
            kind: SourceKind::Assumed,
            citation: "yield at the reference energy and normal incidence; assumed",
            sigma: 0.3,
        },
        ParamSpec {
            id: "yield.peak.angle",
            unit: "deg",
            declared: yield_model::YieldParams::silicon().peak_angle_deg,
            kind: SourceKind::Assumed,
            citation: "oblique angle at which the sputter yield peaks; typical, not fitted",
            sigma: 5.0,
        },
        ParamSpec {
            id: "stress.v.diffusion",
            unit: "cm^3/mol",
            declared: stress_coupling::ActivationVolumes::silicon().diffusion,
            kind: SourceKind::Assumed,
            citation: "activation volume for diffusion under stress; assumed",
            sigma: 2.0,
        },
        ParamSpec {
            id: "stress.v.oxidation",
            unit: "cm^3/mol",
            declared: stress_coupling::ActivationVolumes::silicon().oxidation,
            kind: SourceKind::Assumed,
            citation: "activation volume for the oxidation rate under stress; assumed",
            sigma: 3.0,
        },
        ParamSpec {
            id: "stress.v.mobility",
            unit: "cm^3/mol",
            declared: stress_coupling::ActivationVolumes::silicon().mobility,
            kind: SourceKind::Assumed,
            citation: "activation volume for mobility under stress; assumed",
            sigma: 1.5,
        },
        ParamSpec {
            id: "epi.critical.a",
            unit: "nm",
            declared: epitaxy::CRITICAL_A_NM,
            kind: SourceKind::Assumed,
            citation: "coefficient in h_c = A/x^2 for SiGe; chosen so 20% Ge gives about 10 nm, not a Matthews-Blakeslee calculation",
            sigma: 0.2,
        },
        ParamSpec {
            id: "mask3d.k",
            unit: "nm/nm",
            declared: mask3d::MASK3D_K,
            kind: SourceKind::Assumed,
            citation: "thick-mask focus shift per unit absorber thickness at NA 1; chosen for a tens-of-nm shift, not fitted to a published EMF curve",
            sigma: 0.1,
        },
        ParamSpec {
            id: "pair.initial.supersaturation",
            unit: "-",
            declared: pair_diffusion::PairConfig::implant().initial_supersaturation,
            kind: SourceKind::Assumed,
            citation: "interstitial supersaturation after an implant; the observed order, not fitted",
            sigma: 5.0,
        },
        ParamSpec {
            id: "pair.k.recombine",
            unit: "1/min",
            declared: pair_diffusion::PairConfig::implant().k_recombine,
            kind: SourceKind::Assumed,
            citation: "second-order interstitial-vacancy recombination rate; chosen for a few-minute transient",
            sigma: 0.2,
        },
        ParamSpec {
            id: "ox.massoud.c",
            unit: "um/h",
            declared: massoud::Massoud::bundled().c_um_per_h,
            kind: SourceKind::Assumed,
            citation: "thin-regime acceleration for dry O2 on (100) Si; shape and magnitude assumed, not fitted to a published Massoud curve",
            sigma: 0.006,
        },
        ParamSpec {
            id: "ox.massoud.l",
            unit: "um",
            declared: massoud::Massoud::bundled().l_um,
            kind: SourceKind::Assumed,
            citation: "thin-regime characteristic length; assumed",
            sigma: 0.003,
        },
        ParamSpec {
            id: "sol.n0.b",
            unit: "cm^-3",
            declared: atto_dopant::activation::Solubility::of(atto_params::Species::Boron).n0_cm3,
            kind: SourceKind::Assumed,
            citation: "boron solubility pre-exponential; chosen for ordering and capping behaviour, not fitted to a published table",
            sigma: 1.0e22,
        },
        ParamSpec {
            id: "sol.ea.b",
            unit: "eV",
            declared: atto_dopant::activation::Solubility::of(atto_params::Species::Boron).ea_ev,
            kind: SourceKind::Assumed,
            citation: "boron solubility activation energy; assumed",
            sigma: 0.2,
        },
        ParamSpec {
            id: "sol.n0.as",
            unit: "cm^-3",
            declared: atto_dopant::activation::Solubility::of(atto_params::Species::Arsenic).n0_cm3,
            kind: SourceKind::Assumed,
            citation: "arsenic solubility pre-exponential; assumed",
            sigma: 1.0e22,
        },
        ParamSpec {
            id: "sol.ea.as",
            unit: "eV",
            declared: atto_dopant::activation::Solubility::of(atto_params::Species::Arsenic).ea_ev,
            kind: SourceKind::Assumed,
            citation: "arsenic solubility activation energy; assumed",
            sigma: 0.2,
        },
        ParamSpec {
            id: "sol.n0.p",
            unit: "cm^-3",
            declared: atto_dopant::activation::Solubility::of(atto_params::Species::Phosphorus).n0_cm3,
            kind: SourceKind::Assumed,
            citation: "phosphorus solubility pre-exponential; assumed",
            sigma: 1.0e22,
        },
        ParamSpec {
            id: "sol.ea.p",
            unit: "eV",
            declared: atto_dopant::activation::Solubility::of(atto_params::Species::Phosphorus).ea_ev,
            kind: SourceKind::Assumed,
            citation: "phosphorus solubility activation energy; assumed",
            sigma: 0.2,
        },
        ParamSpec {
            id: "cmp.k.metal",
            unit: "nm/min/(kPa*m/s)",
            declared: cmp::K_METAL,
            kind: SourceKind::Assumed,
            citation: "Preston coefficient of a damascene metal; typical magnitude, not measured here",
            sigma: 8.0,
        },
        ParamSpec {
            id: "cmp.k.dielectric",
            unit: "nm/min/(kPa*m/s)",
            declared: cmp::K_DIELECTRIC,
            kind: SourceKind::Assumed,
            citation: "Preston coefficient of an oxide dielectric; typical magnitude, not measured here",
            sigma: 3.0,
        },
        ParamSpec {
            id: "cmp.contact.length",
            unit: "nm",
            declared: cmp::CONTACT_LENGTH_NM,
            kind: SourceKind::Assumed,
            citation: "pad contact length h_c in the pattern-density model; not measured",
            sigma: 5.0,
        },
        ParamSpec {
            id: "cmp.planarisation.length",
            unit: "nm",
            declared: cmp::PLANARISATION_LENGTH_NM,
            kind: SourceKind::Assumed,
            citation: "planarisation length of the pad-load window; not measured",
            sigma: 10_000.0,
        },
        ParamSpec { id: "passivation.k.dep", unit: "1/min", declared: passivation::K_DEP_PER_MIN, kind: SourceKind::Assumed, citation: "polymer deposition constant at full passivant fraction; not fitted to a measured profile", sigma: 2.0 },
        ParamSpec { id: "passivation.k.rem", unit: "1/min", declared: passivation::K_REM_PER_MIN, kind: SourceKind::Assumed, citation: "ion-assisted polymer removal constant; not fitted to a measured profile", sigma: 0.1 },
    ]
}

/// The value the code computes for an id.
pub fn live(id: &str) -> Option<f64> {
    match id {
        "ox.beak.ratio" => Some(oxidation::BEAK_RATIO),
        "etch.arde.factor" => Some(etch::EtchModel::bundled_oxide().arde_factor),
        "etch.loading.coeff" => Some(etch::EtchModel::bundled_oxide().loading_coeff),
        "etch.sigma.per.mtorr" => Some(etch::EtchModel::bundled_oxide().sigma_per_mtorr),
        "etch.critical.ar" => Some(etch::EtchModel::bundled_oxide().critical_aspect_ratio),
        "defect.threshold.ar" => Some(profile::DefectCoefficients::silicon().threshold_ar),
        "defect.bow.per.ar" => Some(profile::DefectCoefficients::silicon().bow_per_ar),
        "defect.notch.per.ar" => Some(profile::DefectCoefficients::silicon().notch_per_ar),
        "defect.microtrench.per.ar" => Some(profile::DefectCoefficients::silicon().microtrench_per_ar),
        "euv.photon.j" => Some(stochastic::EUV_PHOTON_J),
        "euv.lwr.response" => Some(stochastic::LWR_RESPONSE_NM),
        "seg.m.boron" => Some(atto_dopant::segregation::BORON_SEGREGATION_M),
        "oed.k" => Some(oed::OED_K),
        "ald.nm.per.cycle" => Some(ald::AldProcess::oxide().nm_per_cycle),
        "ald.saturation.dose" => Some(ald::AldProcess::oxide().saturation_dose),
        "yield.threshold.ev" => Some(yield_model::YieldParams::silicon().threshold_ev),
        "yield.reference.ev" => Some(yield_model::YieldParams::silicon().reference_ev),
        "yield.ref" => Some(yield_model::YieldParams::silicon().y_ref),
        "yield.peak.angle" => Some(yield_model::YieldParams::silicon().peak_angle_deg),
        "stress.v.diffusion" => Some(stress_coupling::ActivationVolumes::silicon().diffusion),
        "stress.v.oxidation" => Some(stress_coupling::ActivationVolumes::silicon().oxidation),
        "stress.v.mobility" => Some(stress_coupling::ActivationVolumes::silicon().mobility),
        "epi.critical.a" => Some(epitaxy::CRITICAL_A_NM),
        "mask3d.k" => Some(mask3d::MASK3D_K),
        "pair.initial.supersaturation" => Some(pair_diffusion::PairConfig::implant().initial_supersaturation),
        "pair.k.recombine" => Some(pair_diffusion::PairConfig::implant().k_recombine),
        "ox.massoud.c" => Some(massoud::Massoud::bundled().c_um_per_h),
        "ox.massoud.l" => Some(massoud::Massoud::bundled().l_um),
        "sol.n0.b" => Some(atto_dopant::activation::Solubility::of(atto_params::Species::Boron).n0_cm3),
        "sol.ea.b" => Some(atto_dopant::activation::Solubility::of(atto_params::Species::Boron).ea_ev),
        "sol.n0.as" => Some(atto_dopant::activation::Solubility::of(atto_params::Species::Arsenic).n0_cm3),
        "sol.ea.as" => Some(atto_dopant::activation::Solubility::of(atto_params::Species::Arsenic).ea_ev),
        "sol.n0.p" => Some(atto_dopant::activation::Solubility::of(atto_params::Species::Phosphorus).n0_cm3),
        "sol.ea.p" => Some(atto_dopant::activation::Solubility::of(atto_params::Species::Phosphorus).ea_ev),
        "cmp.k.metal" => Some(cmp::K_METAL),
        "cmp.k.dielectric" => Some(cmp::K_DIELECTRIC),
        "cmp.contact.length" => Some(cmp::CONTACT_LENGTH_NM),
        "cmp.planarisation.length" => Some(cmp::PLANARISATION_LENGTH_NM),
        "passivation.k.dep" => Some(passivation::K_DEP_PER_MIN),
        "passivation.k.rem" => Some(passivation::K_REM_PER_MIN),
        _ => None,
    }
}

#[cfg(test)]
mod registry_tests {
    use super::*;
    use atto_params::verify;

    #[test]
    fn the_process_registry_matches_the_code() {
        let specs = registry();
        assert!(!specs.is_empty());
        let drift = verify(&specs, &live);
        assert!(drift.is_empty(), "{}", drift.iter().map(|d| d.to_string()).collect::<Vec<_>>().join("\n"));
    }

    #[test]
    fn every_process_constant_is_declared() {
        for s in registry() {
            assert!(live(s.id).is_some());
            assert!(!s.citation.is_empty(), "{} has no citation", s.id);
            // An exact constant has no uncertainty: the EUV photon energy is hc/lambda from
            // defined values, so its sigma is zero. The rule is "no uncertainty is silently
            // absent", not "every uncertainty is positive".
            assert!(s.sigma >= 0.0, "{} has a negative uncertainty", s.id);
        }
    }
}
