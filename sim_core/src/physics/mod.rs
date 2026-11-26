//! The 6 Fabrication Modules — deterministic, closed-form/phenomenological
//! physics for each stage of the manufacturing pipeline. Each submodule
//! exposes an `Inputs` struct, an `Outputs` struct, and a `run()` function
//! so the engine's dependency graph (see `crate::engine`) can invoke stages
//! uniformly.

pub mod deposition;
pub mod diffusion;
pub mod doping;
pub mod etching;
pub mod lithography;
pub mod mask_prep;
pub mod metrology;
pub mod oxidation;
pub mod prep;

/// Enumerates the seven pipeline stages in canonical order, for use by the UI
/// process-flow sidebar and the engine's dependency graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum Stage {
    SubstratePrep,
    Deposition,
    Lithography,
    Etching,
    Doping,
    FinalClean,
    Metrology,
    /// Chemical-mechanical planarisation: polish the top film flat.
    Cmp,
}

impl Stage {
    pub const ALL: [Stage; 8] = [
        Stage::SubstratePrep,
        Stage::Deposition,
        Stage::Lithography,
        Stage::Etching,
        Stage::Doping,
        Stage::FinalClean,
        Stage::Metrology,
        Stage::Cmp,
    ];

    /// The stages a new recipe starts with, in order. Planarisation is an optional extra: it is
    /// on the stage menu (`ALL`) but is not one of the default seven, and a recipe that wants it
    /// adds it after doping.
    pub const DEFAULT_RECIPE: [Stage; 7] = [
        Stage::SubstratePrep,
        Stage::Deposition,
        Stage::Lithography,
        Stage::Etching,
        Stage::Doping,
        Stage::FinalClean,
        Stage::Metrology,
    ];

    pub fn label(&self) -> &'static str {
        match self {
            Stage::SubstratePrep => "Substrate Prep",
            Stage::Deposition => "Deposition",
            Stage::Lithography => "Lithography",
            Stage::Etching => "Etching",
            Stage::Doping => "Doping & Anneal",
            Stage::FinalClean => "Final Clean",
            Stage::Metrology => "Metrology", Stage::Cmp => "Planarisation",
        }
    }
}
