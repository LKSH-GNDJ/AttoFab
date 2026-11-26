//! `sim_core`: AttoFab's platform-agnostic deterministic physics & I/O core.
//!
//! This crate has no UI or windowing dependencies, so it compiles identically
//! for desktop, WebAssembly, and Android targets. It owns:
//!   - the Sparse Voxel Octree wafer state (`wafer`),
//!   - the reactive dependency-graph engine and worker-thread plumbing (`engine`),
//!   - the six closed-form/phenomenological physics stages (`physics`),
//!   - GDSII layout ingestion (`io`),
//!   - and SVO -> GPU buffer flattening for the raymarching renderer (`render_wgpu`).
//!
//! Licensed under the GNU General Public License v3.0 or later. See
//! /README.md and /LICENSE for the full license text.

pub mod calibration;
pub mod engine;
pub mod fem;
pub mod geom;
pub mod io;
pub mod legend;
pub mod measure;
pub mod glue;
pub mod gds_sweep;
pub mod params;
pub mod physics;
pub mod pipeline;
pub mod playback;
pub mod render_wgpu;
pub mod report;
pub mod scale;
pub mod session;
pub mod signoff;
pub mod transform;
pub mod validation;
pub mod vv;
pub mod wafer;

/// Re-exports of the most commonly used types, for `use sim_core::prelude::*;`.
pub mod prelude {
    pub use crate::engine::{Engine, StageJob, StageResult};
    pub use crate::physics::Stage;
    pub use crate::pipeline::{evaluate_stage, StageEval};
    pub use crate::wafer::{Material, WaferState};
}
