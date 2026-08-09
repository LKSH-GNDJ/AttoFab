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
//!
//! Phase 3 adds `gds` (GDSII layout parsing) and `grid3d` (dense 3D wafer
//! mesh - RLE/chunked sparsity deferred per the "start dense, prove
//! correctness" plan, same as Phase 1/2 were).

pub mod diffusion;
pub mod etch;
pub mod gds;
pub mod grid2d;
pub mod grid3d;
pub mod implant;
pub mod lithography;
pub mod materials;
pub mod oxidation;
pub mod wafer1d;

pub use gds::{parse_gds, GdsLibrary};
pub use grid2d::Wafer2d;
pub use grid3d::Wafer3d;
pub use materials::{Ambient, Dopant, Material};
pub use wafer1d::{AnnealStep, ImplantStep, OxidizeStep, Substrate, Wafer1d};
