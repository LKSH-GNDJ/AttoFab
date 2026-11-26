// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The device layer's constants, in the same manifest form as the process layer (K-01).
//!
//! K-01 says no constant may be unregistered. Until now that rule could not be satisfied for
//! anything in this crate: the registry lived in `atto_process`, and `atto_device` depends on
//! `atto_process`, not the other way round — so a constant here had nowhere to be declared.
//!
//! The manifest types (`ParamSpec`, `SourceKind`) live in `atto_params`, which is upstream of
//! both. This module uses them directly, so a device constant is now registrable and a check
//! can aggregate the two layers without a dependency inversion.

use atto_params::{ParamSpec, SourceKind};
use crate::heavy::{E_REF_EV, N_REF_CM3};
use crate::mos::QE;

/// Every constant this crate introduces, with its provenance.
pub fn registry() -> Vec<ParamSpec> {
    vec![
        ParamSpec {
            id: "si.qe",
            unit: "C",
            declared: QE,
            kind: SourceKind::Published,
            citation: "elementary charge, CODATA 2018 exact value",
            sigma: 0.0,
        },
        ParamSpec {
            id: "sheet.mu.default",
            unit: "cm^2/Vs",
            declared: crate::sheet::DEFAULT_MOBILITY_CM2_VS,
            kind: SourceKind::Assumed,
            citation: "a typical drift mobility for a moderately doped silicon layer; not measured here and not species-dependent",
            sigma: 150.0,
        },
        ParamSpec {
            id: "meas.series.r",
            unit: "ohm",
            declared: crate::small_signal::DEFAULT_SERIES_R_OHM,
            kind: SourceKind::Assumed,
            citation: "representative series resistance for a C-V measurement; a measurement condition, not a material property",
            sigma: 50.0,
        },
        ParamSpec {
            id: "aval.a",
            unit: "cm^-1",
            declared: crate::avalanche::Ionisation::silicon().a_per_cm,
            kind: SourceKind::Assumed,
            citation: "Chynoweth ionisation pre-exponential for electrons in Si; commonly quoted, not fitted here",
            sigma: 1.0e5,
        },
        ParamSpec {
            id: "aval.b",
            unit: "V/cm",
            declared: crate::avalanche::Ionisation::silicon().b_v_per_cm,
            kind: SourceKind::Assumed,
            citation: "Chynoweth critical field for electrons in Si; commonly quoted, not fitted here",
            sigma: 0.2e6,
        },
        ParamSpec {
            id: "heavy.nref",
            unit: "cm^-3",
            declared: N_REF_CM3,
            kind: SourceKind::Published,
            citation: "Slotboom-de Graaff reference doping for bandgap narrowing",
            sigma: 0.0,
        },
        ParamSpec {
            id: "heavy.eref",
            unit: "eV",
            declared: E_REF_EV,
            kind: SourceKind::Published,
            citation: "Slotboom-de Graaff reference energy; the form of the expression is theirs, the fit is not checked here",
            sigma: 1.0e-3,
        },
    ]
}

/// Every registered constant in the tree, from both layers.
///
/// K-01 says no constant may be unregistered, and until this existed the rule could only be
/// checked one crate at a time: a reader had to know to look in both `atto_process::registry()`
/// and here. This is the list the rule actually means.
pub fn all() -> Vec<ParamSpec> {
    let mut v = atto_process::registry();
    v.extend(registry());
    v
}

/// The declared value of a registered constant, by id, from either layer.
pub fn value_anywhere(id: &str) -> Option<f64> {
    value_of(id).or_else(|| {
        // the process layer has its own accessor; probe it by looking the id up in its registry
        atto_process::registry()
            .into_iter()
            .find(|s| s.id == id)
            .map(|s| s.declared)
    })
}

/// The declared value of a registered constant, by id.
pub fn value_of(id: &str) -> Option<f64> {
    match id {
        "si.qe" => Some(QE),
        "sheet.mu.default" => Some(crate::sheet::DEFAULT_MOBILITY_CM2_VS),
        "meas.series.r" => Some(crate::small_signal::DEFAULT_SERIES_R_OHM),
        "aval.a" => Some(crate::avalanche::Ionisation::silicon().a_per_cm),
        "aval.b" => Some(crate::avalanche::Ionisation::silicon().b_v_per_cm),
        "heavy.nref" => Some(N_REF_CM3),
        "heavy.eref" => Some(E_REF_EV),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_declared_constant_resolves_to_its_value() {
        for spec in registry() {
            let v = value_of(spec.id);
            assert!(v.is_some(), "{} is declared but has no value", spec.id);
            assert!((v.unwrap() - spec.declared).abs() <= 1e-12 * spec.declared.abs().max(1.0), "{} resolves to a different number", spec.id);
        }
    }

    #[test]
    fn every_value_has_a_declaration() {
        let ids: Vec<&str> = registry().iter().map(|s| s.id).collect();
        for id in ["si.qe", "heavy.nref", "heavy.eref", "aval.a", "aval.b", "sheet.mu.default"] {
            assert!(ids.contains(&id), "{id} has a value but no declaration");
        }
    }

    #[test]
    fn an_unknown_id_is_not_invented() {
        assert!(value_of("heavy.not.a.constant").is_none());
    }

    #[test]
    fn no_constant_is_declared_twice() {
        let mut seen = std::collections::BTreeSet::new();
        for spec in registry() {
            assert!(seen.insert(spec.id), "{} is declared twice", spec.id);
        }
    }
}

#[cfg(test)]
mod aggregate_tests {
    use super::*;

    #[test]
    fn every_constant_in_the_tree_is_in_one_list() {
        let all = all();
        assert!(all.len() >= registry().len(), "the aggregate includes this crate");
        assert!(all.iter().any(|s| s.id == "si.qe"), "and the device layer's own");
        assert!(all.iter().any(|s| s.id == "ox.massoud.c"), "and the process layer's");
    }

    #[test]
    fn no_constant_is_declared_twice_across_layers() {
        // the invariant K-01 needs and could not check while the lists were separate
        let all = all();
        let mut seen = std::collections::BTreeSet::new();
        for s in &all {
            assert!(seen.insert(s.id), "{} is declared in more than one place", s.id);
        }
    }

    #[test]
    fn every_aggregated_constant_carries_a_citation() {
        for s in all() {
            assert!(!s.citation.is_empty(), "{} has no citation", s.id);
        }
    }

    #[test]
    fn a_value_can_be_read_from_either_layer() {
        assert!(value_anywhere("si.qe").is_some(), "a device constant");
        assert!(value_anywhere("ox.massoud.c").is_some(), "a process constant");
        assert!(value_anywhere("not.a.constant").is_none(), "and an unknown one is not invented");
    }
}
