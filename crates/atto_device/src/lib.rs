// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The device layer (build PRD Track H, V-01 to V-05).
//!
//! This crate reads the three things the process layer produces — geometry, net
//! doping and interface charge — onto a single cross-section mesh, and turns
//! them into electrical quantities. The 2-D drift-diffusion solve is the
//! destination; the depletion approximation is the first station, and it is
//! enough to prove that the mesh carries what a device needs.

pub mod block;
pub mod cross_section;
pub mod diode;
pub mod fastpath;
pub mod avalanche;
pub mod heavy;
pub mod continuity;
pub mod small_signal;
pub mod interface_effect;
pub mod registry;
pub mod mos;
pub mod sheet;
pub mod newton;
pub mod poisson;
pub mod volume3d;

pub use cross_section::{CrossSection, CrossSectionError};
pub use diode::{Diode, DiodePoint};
pub use mos::{MosParams, ThresholdResult};
pub use poisson::{cv_curve, CvPoint, PoissonColumn, PoissonError};

// ---------------------------------------------------------------------------
// Parameter registry (build PRD K-01)
// ---------------------------------------------------------------------------

use atto_params::{ParamSpec, SourceKind};

/// The physical constants the device layer reads.
pub fn registry() -> Vec<ParamSpec> {
    vec![
        ParamSpec { id: "si.eps.r", unit: "-", declared: mos::EPS_SI_R, kind: SourceKind::Published, citation: "Si relative permittivity", sigma: 0.05 },
        ParamSpec { id: "sio2.eps.r", unit: "-", declared: mos::EPS_OX_R, kind: SourceKind::Published, citation: "thermal SiO2 relative permittivity", sigma: 0.05 },
        ParamSpec { id: "si.ni.cm3", unit: "cm^-3", declared: mos::NI, kind: SourceKind::Published, citation: "intrinsic carrier concentration at 300 K", sigma: 1.0e9 },
        ParamSpec { id: "elementary.charge.c", unit: "C", declared: mos::QE, kind: SourceKind::Published, citation: "CODATA", sigma: 1.0e-28 },
        ParamSpec { id: "boltzmann.j_per_k", unit: "J/K", declared: mos::KB, kind: SourceKind::Published, citation: "CODATA", sigma: 1.0e-28 },
        ParamSpec { id: "device.ni.prefactor", unit: "cm^-3 K^-3/2", declared: newton::NI_PREFACTOR, kind: SourceKind::Assumed, citation: "calibrated so ni(300 K) = 1.0e10; the textbook 3.87e16 gives 8.6e10", sigma: 4.0e14 },
        ParamSpec { id: "mos.phi_ms.v", unit: "V", declared: MosParams::default().phi_ms, kind: SourceKind::Assumed, citation: "n+ polysilicon on p-silicon work function difference", sigma: 0.1 },
    ]
}

/// The value the code computes for an id.
pub fn live(id: &str) -> Option<f64> {
    match id {
        "si.eps.r" => Some(mos::EPS_SI_R),
        "sio2.eps.r" => Some(mos::EPS_OX_R),
        "si.ni.cm3" => Some(mos::NI),
        "elementary.charge.c" => Some(mos::QE),
        "boltzmann.j_per_k" => Some(mos::KB),
        "mos.phi_ms.v" => Some(MosParams::default().phi_ms),
        "device.ni.prefactor" => Some(newton::NI_PREFACTOR),
        _ => None,
    }
}

#[cfg(test)]
mod registry_tests {
    use super::*;
    use atto_params::verify;

    #[test]
    fn the_device_registry_matches_the_code() {
        let specs = registry();
        let drift = verify(&specs, &live);
        assert!(drift.is_empty(), "{}", drift.iter().map(|d| d.to_string()).collect::<Vec<_>>().join("\n"));
    }

    #[test]
    fn every_device_constant_is_declared() {
        for s in registry() {
            assert!(live(s.id).is_some(), "{}", s.id);
            assert!(!s.citation.is_empty());
            assert!(s.sigma > 0.0);
        }
    }
}
