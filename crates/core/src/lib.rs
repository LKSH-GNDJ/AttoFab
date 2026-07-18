//! attofab-core: AttoFab's physics engine core.
//!
//! Pure Rust, no FFI bindings (those live in `crates/core-wasm` /
//! `crates/core-py`). Ports the validated Phase-1 JS engine's math
//! (Deal-Grove oxidation, Gaussian implant, Fick's-Law FDM diffusion)
//! with regression tests pinned to that JS engine's verified output.
//!
//! Phase 2 adds `grid2d` (2D wafer cross-section), `lithography` (mask
//! expose/develop), and `etch` (anisotropic/isotropic etch + deposition),
//! built on the same validated Phase 1 solvers.

pub mod diffusion;
pub mod etch;
pub mod grid2d;
pub mod implant;
pub mod lithography;
pub mod materials;
pub mod oxidation;
pub mod wafer1d;

pub use grid2d::Wafer2d;
pub use materials::{Ambient, Dopant, Material};
pub use wafer1d::{AnnealStep, ImplantStep, OxidizeStep, Substrate, Wafer1d};
