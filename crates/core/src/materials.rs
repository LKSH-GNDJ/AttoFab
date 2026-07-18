//! Material Properties Database.
//!
//! Public-domain, textbook-typical physical constants (Deal & Grove 1965;
//! Plummer/Deen/Griffin "Silicon VLSI Technology"; Sze "Physics of
//! Semiconductor Devices"). These are illustrative reference values for an
//! EDUCATIONAL engine, not measured/fitted foundry process data. Swap in
//! your own published constants for higher fidelity, or to benchmark
//! against an open PDK (e.g. SkyWater SKY130 documentation).

use serde::{Deserialize, Serialize};

/// Boltzmann constant, eV/K.
pub const K_EV: f64 = 8.617_333_262e-5;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Material {
    Silicon,
    Oxide,
    Photoresist,
    PhotoresistExposed,
    Metal,
    Void,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Ambient {
    Dry,
    Wet,
}

/// Deal-Grove Arrhenius kinetics: B = C_b * exp(-Ea_b / kT) [um^2/hr];
/// B/A = C_ba * exp(-Ea_ba/kT) [um/hr]. Values are for the (100) plane.
pub struct OxidationKinetics {
    pub c_b: f64,
    pub ea_b: f64,
    pub c_ba: f64,
    pub ea_ba: f64,
}

/// Representative textbook constants for (100) single-crystal Si under
/// atmospheric pressure. See docs/math_references.md for provenance -
/// these were corrected from an earlier internal draft against a more
/// carefully sourced aggregated table.
pub fn oxidation_kinetics(ambient: Ambient) -> OxidationKinetics {
    match ambient {
        Ambient::Dry => OxidationKinetics {
            c_b: 7.72e2,  // um^2/hr
            ea_b: 1.23,   // eV
            c_ba: 3.71e6, // um/hr
            ea_ba: 2.00,  // eV
        },
        Ambient::Wet => OxidationKinetics {
            c_b: 3.86e2,  // um^2/hr
            ea_b: 0.78,   // eV
            c_ba: 9.70e7, // um/hr
            ea_ba: 2.05,  // eV
        },
    }
}

/// (111)-plane linear rate constant (B/A) scaling factor relative to
/// (100). The linear rate is reaction-limited (interface chemistry), which
/// scales with available Si bond density - (111) is more densely packed,
/// so it oxidizes faster in the linear regime. The parabolic constant B is
/// diffusion-limited (bulk transport through amorphous SiO2) and is
/// orientation-independent - this factor applies to B/A only, never to B.
pub const ORIENTATION_111_TO_100_LINEAR_RATIO: f64 = 1.68;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CrystalOrientation {
    Si100,
    Si111,
}

/// Native oxide thickness typically present before a "clean" thermal step (um).
pub const NATIVE_OXIDE_UM: f64 = 0.002;

/// Fraction of a unit of oxide growth that consumes underlying silicon
/// (classic Deal-Grove fact: ~0.44 units Si per unit SiO2 grown).
pub const SI_CONSUMPTION_RATIO: f64 = 0.44;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Hash)]
pub enum Dopant {
    Boron,
    Phosphorus,
    Arsenic,
    Antimony,
}

impl Dopant {
    /// P-type (acceptor) vs N-type (donor) convention used for junction finding.
    pub fn is_p_type(&self) -> bool {
        matches!(self, Dopant::Boron)
    }

    pub fn name(&self) -> &'static str {
        match self {
            Dopant::Boron => "Boron",
            Dopant::Phosphorus => "Phosphorus",
            Dopant::Arsenic => "Arsenic",
            Dopant::Antimony => "Antimony",
        }
    }
}

/// Diffusion coefficient Arrhenius parameters: D(T) = D0 * exp(-Ea/kT) [cm^2/s]
pub struct DiffusionParams {
    pub d0: f64,
    pub ea: f64,
}

pub fn diffusion_params(dopant: Dopant) -> DiffusionParams {
    match dopant {
        Dopant::Boron => DiffusionParams { d0: 0.76, ea: 3.46 },
        Dopant::Phosphorus => DiffusionParams { d0: 3.85, ea: 3.66 },
        Dopant::Arsenic => DiffusionParams { d0: 0.0724, ea: 3.44 },
        Dopant::Antimony => DiffusionParams { d0: 0.214, ea: 3.65 },
    }
}

/// Simple illustrative energy-scaling for implant range/straggle.
/// NOT a substitute for a SRIM/TRIM Monte-Carlo range table - pass
/// measured Rp/dRp directly when you have them.
pub struct RangeScaling {
    pub rp_per_kev_um: f64,
    pub straggle_fraction: f64,
}

pub fn range_scaling(dopant: Dopant) -> RangeScaling {
    match dopant {
        Dopant::Boron => RangeScaling { rp_per_kev_um: 0.0028, straggle_fraction: 0.42 },
        Dopant::Phosphorus => RangeScaling { rp_per_kev_um: 0.0016, straggle_fraction: 0.35 },
        Dopant::Arsenic => RangeScaling { rp_per_kev_um: 0.0009, straggle_fraction: 0.30 },
        Dopant::Antimony => RangeScaling { rp_per_kev_um: 0.0007, straggle_fraction: 0.28 },
    }
}

/// Arrhenius helper: rate = prefactor * exp(-Ea / (k * T))
pub fn arrhenius(prefactor: f64, ea_ev: f64, t_kelvin: f64) -> f64 {
    prefactor * (-ea_ev / (K_EV * t_kelvin)).exp()
}
