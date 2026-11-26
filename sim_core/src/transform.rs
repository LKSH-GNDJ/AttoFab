// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Wafer environment transform (PRD §1).
//!
//! The wafer's SVO geometry always lives in *wafer-local* coordinates. Moving
//! or rotating the wafer in the environment only changes this transform; no
//! SVO node is rewritten. Rendering, picking and placement all convert through
//! it: camera rays are taken into wafer-local space for picking, and hit
//! points are brought back to world space for display.
//!
//! Orientation is a rotation about world Z (`yaw_rad`) only, so the wafer
//! always stays circular and horizontal in its local XY plane.

use serde::{Deserialize, Serialize};

pub type Vec3 = [f64; 3];

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct WaferTransform {
    /// World-space position of the wafer-local origin, nm.
    pub position: Vec3,
    /// Rotation about world +Z, radians (counter-clockwise seen from +Z).
    pub yaw_rad: f64,
}

impl Default for WaferTransform {
    fn default() -> Self { Self { position: [0.0; 3], yaw_rad: 0.0 } }
}

impl WaferTransform {
    pub const IDENTITY: WaferTransform = WaferTransform { position: [0.0; 3], yaw_rad: 0.0 };

    pub fn is_finite(&self) -> bool { self.position.iter().all(|v| v.is_finite()) && self.yaw_rad.is_finite() }

    /// Rotate a vector from wafer-local to world orientation (no translation).
    pub fn dir_to_world(&self, d: Vec3) -> Vec3 {
        let (s, c) = self.yaw_rad.sin_cos();
        [c * d[0] - s * d[1], s * d[0] + c * d[1], d[2]]
    }

    /// Rotate a vector from world to wafer-local orientation (no translation).
    pub fn dir_to_local(&self, d: Vec3) -> Vec3 {
        let (s, c) = self.yaw_rad.sin_cos();
        [c * d[0] + s * d[1], -s * d[0] + c * d[1], d[2]]
    }

    pub fn local_to_world(&self, p: Vec3) -> Vec3 {
        let r = self.dir_to_world(p);
        [r[0] + self.position[0], r[1] + self.position[1], r[2] + self.position[2]]
    }

    pub fn world_to_local(&self, p: Vec3) -> Vec3 {
        self.dir_to_local([p[0] - self.position[0], p[1] - self.position[1], p[2] - self.position[2]])
    }

    /// Transform a world-space ray into wafer-local space (direction stays unit length).
    pub fn ray_to_local(&self, origin: Vec3, dir: Vec3) -> (Vec3, Vec3) {
        (self.world_to_local(origin), self.dir_to_local(dir))
    }

    /// Translate in world space. Non-finite deltas are rejected (returns false).
    pub fn translate_world(&mut self, delta: Vec3) -> bool {
        if delta.iter().any(|v| !v.is_finite()) { return false; }
        for i in 0..3 { self.position[i] += delta[i]; }
        true
    }
}

/// Intersect a ray with the plane through `p0` with normal `n`.
/// Returns the ray parameter `t > 0`, or `None` when the ray is (nearly)
/// parallel to the plane, points away from it, or any input is non-finite.
pub fn ray_plane_t(origin: Vec3, dir: Vec3, n: Vec3, p0: Vec3) -> Option<f64> {
    if origin.iter().chain(dir.iter()).chain(n.iter()).chain(p0.iter()).any(|v| !v.is_finite()) { return None; }
    let denom = dir[0] * n[0] + dir[1] * n[1] + dir[2] * n[2];
    let dir_len = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
    let n_len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
    if dir_len < 1e-12 || n_len < 1e-12 || denom.abs() < 1e-9 * dir_len * n_len { return None; }
    let t = ((p0[0] - origin[0]) * n[0] + (p0[1] - origin[1]) * n[1] + (p0[2] - origin[2]) * n[2]) / denom;
    if t.is_finite() && t > 0.0 { Some(t) } else { None }
}

pub fn ray_point(origin: Vec3, dir: Vec3, t: f64) -> Vec3 {
    [origin[0] + dir[0] * t, origin[1] + dir[1] * t, origin[2] + dir[2] * t]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trip_is_identity_for_translation_and_yaw() {
        let t = WaferTransform { position: [1.0e6, -2.5e6, 300.0], yaw_rad: 0.7 };
        let p = [12_345.0, -9_876.0, 42.0];
        let q = t.world_to_local(t.local_to_world(p));
        for i in 0..3 { assert!((p[i] - q[i]).abs() < 1e-6, "axis {i}: {p:?} vs {q:?}"); }
    }

    #[test]
    fn yaw_keeps_wafer_horizontal() {
        let t = WaferTransform { position: [0.0; 3], yaw_rad: 1.1 };
        assert!((t.dir_to_world([0.0, 0.0, 1.0])[2] - 1.0).abs() < 1e-12);
        let up = t.dir_to_world([0.0, 0.0, 1.0]);
        assert!(up[0].abs() < 1e-12 && up[1].abs() < 1e-12);
    }

    #[test]
    fn ray_plane_rejects_parallel_backward_and_non_finite() {
        let o = [0.0, 0.0, 10.0];
        assert!(ray_plane_t(o, [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0; 3]).is_none(), "parallel");
        assert!(ray_plane_t(o, [0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0; 3]).is_none(), "pointing away");
        assert!(ray_plane_t(o, [0.0, 0.0, f64::NAN], [0.0, 0.0, 1.0], [0.0; 3]).is_none(), "non-finite");
        assert_eq!(ray_plane_t(o, [0.0, 0.0, -1.0], [0.0, 0.0, 1.0], [0.0; 3]), Some(10.0));
    }

    #[test]
    fn non_finite_translation_is_rejected_without_change() {
        let mut t = WaferTransform::default();
        assert!(!t.translate_world([f64::INFINITY, 0.0, 0.0]));
        assert_eq!(t, WaferTransform::default());
        assert!(t.translate_world([5.0, -3.0, 0.0]));
        assert_eq!(t.position, [5.0, -3.0, 0.0]);
    }
}
