// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Geometry: the sparse voxel octree and the column operators (build PRD A-06).
//!
//! The SVO is unchanged in spirit from the earlier build: homogeneous regions collapse to a
//! single leaf, edits subdivide only where they must, and the tree re-collapses
//! afterwards. The operators are the the earlier build column model with the audit's defects
//! removed: the pitch is the true leaf size, conformal growth is separable
//! rather than O(r^2) per column, and etch depth excludes consumed resist.

pub mod levelset;
pub mod svo;
pub mod column;
pub mod gds;
pub mod oasis;
pub mod mask;
pub mod refine;
pub mod tiling;

pub use svo::{Bounds, Material, Run, SvoNode, Wafer};
pub use column::{ColumnGrid, ProcessWindow};
