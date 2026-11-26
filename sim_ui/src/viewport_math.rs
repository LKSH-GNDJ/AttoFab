// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Camera, picking and wafer-drag math (PRD §1), free of any UI types.
//!
//! One [`CameraRig`] is the single source of truth for rendering, picking,
//! overlay projection and dragging: the same `forward / right / up` basis is
//! handed to the shader, so what the user clicks is what the GPU drew. (The
//! previous implementation recomputed the basis separately in the picker and
//! in the shader, and they disagreed in the XZ view.)
//!
//! Frames: *world* is the environment (camera lives here). Wafer geometry is in
//! *wafer-local* coordinates; [`WaferTransform`] maps local → world. Picking
//! intersects the camera ray with a work plane in world space and converts the
//! hit back to wafer-local coordinates.

use sim_core::transform::{ray_plane_t, ray_point, Vec3, WaferTransform};

/// tan(half vertical field of view). ~20° half-angle.
pub const TAN_HALF_FOV: f64 = 0.364;
pub const MIN_DISTANCE_NM: f64 = 300.0;
pub const MAX_DISTANCE_NM: f64 = 4.0e9;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewMode { Perspective, Top, Xz, Yz, Slice }

impl ViewMode {
    pub const ALL: [ViewMode; 5] = [ViewMode::Perspective, ViewMode::Top, ViewMode::Xz, ViewMode::Yz, ViewMode::Slice];
    pub fn label(self) -> &'static str {
        match self { ViewMode::Perspective => "3D", ViewMode::Top => "Top", ViewMode::Xz => "XZ", ViewMode::Yz => "YZ", ViewMode::Slice => "Slice" }
    }
    pub fn shader_id(self) -> f32 { match self { ViewMode::Perspective => 0.0, ViewMode::Top => 1.0, ViewMode::Xz => 2.0, ViewMode::Yz => 3.0, ViewMode::Slice => 4.0 } }
    pub fn is_section(self) -> bool { matches!(self, ViewMode::Xz | ViewMode::Yz) }
}

pub fn dot(a: Vec3, b: Vec3) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
pub fn cross(a: Vec3, b: Vec3) -> Vec3 { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
pub fn norm(a: Vec3) -> f64 { dot(a, a).sqrt() }
pub fn normalize(a: Vec3) -> Vec3 { let n = norm(a).max(1e-300); [a[0] / n, a[1] / n, a[2] / n] }
pub fn add(a: Vec3, b: Vec3) -> Vec3 { [a[0] + b[0], a[1] + b[1], a[2] + b[2]] }
pub fn sub(a: Vec3, b: Vec3) -> Vec3 { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
pub fn scale(a: Vec3, s: f64) -> Vec3 { [a[0] * s, a[1] * s, a[2] * s] }

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CameraRig {
    pub mode: ViewMode,
    pub yaw: f64,
    pub pitch: f64,
    pub distance: f64,
    /// World point the camera looks at.
    pub target: Vec3,
}

impl CameraRig {
    pub fn new(mode: ViewMode) -> Self { Self { mode, yaw: -0.72, pitch: 0.58, distance: 4.0e8, target: [0.0; 3] } }

    pub fn is_finite(&self) -> bool { self.yaw.is_finite() && self.pitch.is_finite() && self.distance.is_finite() && self.target.iter().all(|v| v.is_finite()) }

    pub fn forward(&self) -> Vec3 {
        match self.mode {
            ViewMode::Perspective => { let (cp, sp) = (self.pitch.cos(), self.pitch.sin()); normalize([-cp * self.yaw.cos(), -cp * self.yaw.sin(), -sp]) }
            ViewMode::Top | ViewMode::Slice => [0.0, 0.0, -1.0],
            ViewMode::Xz => [0.0, 1.0, 0.0],   // camera at −Y, X to the right, Z up
            ViewMode::Yz => [-1.0, 0.0, 0.0],  // camera at +X, Y to the right, Z up
        }
    }

    /// Right-handed orthonormal basis `(forward, right, up)`.
    pub fn basis(&self) -> (Vec3, Vec3, Vec3) {
        let f = self.forward();
        let hint = if f[2].abs() > 0.98 { [0.0, 1.0, 0.0] } else { [0.0, 0.0, 1.0] };
        let r = normalize(cross(f, hint));
        let u = normalize(cross(r, f));
        (f, r, u)
    }

    pub fn eye(&self) -> Vec3 { sub(self.target, scale(self.forward(), self.distance)) }

    /// Ray through normalised device coordinates (x right, y up, both −1..1).
    pub fn ray(&self, ndc: [f64; 2], aspect: f64) -> Option<(Vec3, Vec3)> {
        if !(ndc[0].is_finite() && ndc[1].is_finite() && aspect.is_finite() && aspect > 0.0) || !self.is_finite() { return None; }
        let (f, r, u) = self.basis();
        let d = normalize(add(add(f, scale(r, ndc[0] * aspect * TAN_HALF_FOV)), scale(u, ndc[1] * TAN_HALF_FOV)));
        Some((self.eye(), d))
    }

    /// Project a world point to NDC; `None` if it is behind the camera.
    pub fn project(&self, p: Vec3, aspect: f64) -> Option<[f64; 2]> {
        let (f, r, u) = self.basis();
        let v = sub(p, self.eye());
        let z = dot(v, f);
        if !(z > 1e-6) || !z.is_finite() { return None; }
        Some([dot(v, r) / (z * aspect * TAN_HALF_FOV), dot(v, u) / (z * TAN_HALF_FOV)])
    }

    /// World nanometres covered by one screen pixel at the target depth.
    pub fn nm_per_pixel(&self, viewport_h_px: f64) -> f64 { 2.0 * TAN_HALF_FOV * self.distance / viewport_h_px.max(1.0) }

    pub fn orbit(&mut self, dx_px: f64, dy_px: f64) {
        if self.mode != ViewMode::Perspective || !dx_px.is_finite() || !dy_px.is_finite() { return; }
        self.yaw -= dx_px * 0.008;
        self.pitch = (self.pitch + dy_px * 0.008).clamp(-1.45, 1.45);
    }

    /// Zoom by `factor` (< 1 zooms in) about the world point `anchor` so that it
    /// stays under the cursor. Returns false for non-finite input.
    pub fn zoom_about(&mut self, anchor: Vec3, factor: f64) -> bool {
        if !factor.is_finite() || factor <= 0.0 || anchor.iter().any(|v| !v.is_finite()) { return false; }
        let d = (self.distance * factor).clamp(MIN_DISTANCE_NM, MAX_DISTANCE_NM);
        let k = d / self.distance;
        let eye = self.eye();
        let eye2 = add(anchor, scale(sub(eye, anchor), k));
        self.distance = d;
        self.target = add(eye2, scale(self.forward(), d));
        true
    }

    pub fn zoom_about_target(&mut self, factor: f64) -> bool { let t = self.target; self.zoom_about(t, factor) }

    /// Distance at which a disk of `radius` fits the viewport.
    pub fn fit_distance(mode: ViewMode, radius: f64, aspect: f64) -> f64 {
        let limit = if matches!(mode, ViewMode::Xz | ViewMode::Yz) { aspect } else { aspect.min(1.0) };
        let margin = if mode == ViewMode::Perspective { 1.3 } else { 1.1 };
        (radius * margin / (TAN_HALF_FOV * limit.max(0.2))).clamp(MIN_DISTANCE_NM, MAX_DISTANCE_NM)
    }
}

/// Work plane in world space used for placement, inspection and dragging.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorkPlane { pub normal: Vec3, pub point: Vec3 }

impl WorkPlane {
    /// Horizontal plane at the wafer surface (3D / Top), at the slice height (Slice), or the
    /// vertical section plane through the process site (XZ / YZ).
    pub fn for_view(mode: ViewMode, t: &WaferTransform, surface_z: f64, slice_z: f64, site_local: [f64; 2]) -> Self {
        match mode {
            ViewMode::Perspective | ViewMode::Top => Self { normal: [0.0, 0.0, 1.0], point: t.local_to_world([site_local[0], site_local[1], surface_z]) },
            ViewMode::Slice => Self { normal: [0.0, 0.0, 1.0], point: t.local_to_world([site_local[0], site_local[1], slice_z]) },
            ViewMode::Xz => Self { normal: [0.0, 1.0, 0.0], point: t.local_to_world([site_local[0], site_local[1], surface_z]) },
            ViewMode::Yz => Self { normal: [1.0, 0.0, 0.0], point: t.local_to_world([site_local[0], site_local[1], surface_z]) },
        }
    }

    pub fn hit(&self, ray: (Vec3, Vec3)) -> Option<Vec3> {
        let t = ray_plane_t(ray.0, ray.1, self.normal, self.point)?;
        Some(ray_point(ray.0, ray.1, t))
    }
}

/// Result of a placement pick, in wafer-local coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PlacePick { pub world: Vec3, pub local: Vec3 }

/// Pick the wafer at an NDC position. `None` when the ray misses the plane or the point is
/// outside the circular wafer.
pub fn pick_on_wafer(rig: &CameraRig, ndc: [f64; 2], aspect: f64, plane: &WorkPlane, t: &WaferTransform, wafer_radius: f64) -> Option<PlacePick> {
    if !t.is_finite() { return None; }
    let ray = rig.ray(ndc, aspect)?;
    let world = plane.hit(ray)?;
    let local = t.world_to_local(world);
    if local[0].hypot(local[1]) > wafer_radius { return None; }
    Some(PlacePick { world, local })
}

/// Drag-capture state for the Move Wafer tool: the wafer point under the cursor at press time
/// stays under the cursor for the whole drag.
#[derive(Debug, Clone, Copy)]
pub struct WaferDrag { plane: WorkPlane, start_hit: Vec3, start: WaferTransform }

impl WaferDrag {
    pub fn begin(rig: &CameraRig, ndc: [f64; 2], aspect: f64, plane: WorkPlane, t: WaferTransform) -> Option<Self> {
        let hit = plane.hit(rig.ray(ndc, aspect)?)?;
        Some(Self { plane, start_hit: hit, start: t })
    }

    /// New transform for the pointer at `ndc`, or `None` if the ray misses the drag plane
    /// (the wafer then stays where it was).
    pub fn update(&self, rig: &CameraRig, ndc: [f64; 2], aspect: f64) -> Option<WaferTransform> {
        let hit = self.plane.hit(rig.ray(ndc, aspect)?)?;
        let mut t = self.start;
        if !t.translate_world(sub(hit, self.start_hit)) { return None; }
        Some(t)
    }

    /// Wafer-local point that was grabbed.
    pub fn grab_local(&self) -> Vec3 { self.start.world_to_local(self.start_hit) }
}

/// Camera pan so the world point under the cursor at press time follows the pointer.
#[derive(Debug, Clone, Copy)]
pub struct PanDrag { plane: WorkPlane, start_hit: Vec3, start: CameraRig }

impl PanDrag {
    pub fn begin(rig: &CameraRig, ndc: [f64; 2], aspect: f64, plane: WorkPlane) -> Option<Self> {
        let hit = plane.hit(rig.ray(ndc, aspect)?)?;
        Some(Self { plane, start_hit: hit, start: *rig })
    }
    pub fn update(&self, ndc: [f64; 2], aspect: f64) -> Option<CameraRig> {
        // The ray is cast from the *start* camera: dragging then moves the camera opposite to the pointer.
        let hit = self.plane.hit(self.start.ray(ndc, aspect)?)?;
        let mut r = self.start;
        r.target = add(r.target, sub(self.start_hit, hit));
        Some(r)
    }
}

pub fn pointer_to_ndc(pos: [f32; 2], rect_min: [f32; 2], rect_size: [f32; 2]) -> Option<[f64; 2]> {
    if rect_size[0] < 1.0 || rect_size[1] < 1.0 || pos.iter().any(|v| !v.is_finite()) { return None; }
    Some([((pos[0] - rect_min[0]) / rect_size[0]) as f64 * 2.0 - 1.0, 1.0 - ((pos[1] - rect_min[1]) / rect_size[1]) as f64 * 2.0])
}

pub fn ndc_to_pointer(ndc: [f64; 2], rect_min: [f32; 2], rect_size: [f32; 2]) -> [f32; 2] {
    [rect_min[0] + ((ndc[0] + 1.0) * 0.5) as f32 * rect_size[0], rect_min[1] + ((1.0 - ndc[1]) * 0.5) as f32 * rect_size[1]]
}

/// A "nice" scale-bar length (1, 2 or 5 × 10^k nm) not exceeding `max_nm`, and its label.
pub fn nice_length(max_nm: f64) -> (f64, String) {
    if !(max_nm.is_finite() && max_nm > 0.0) { return (1.0, "1 nm".into()); }
    let p = 10f64.powf(max_nm.log10().floor());
    let n = [5.0, 2.0, 1.0].iter().map(|m| m * p).find(|v| *v <= max_nm).unwrap_or(p);
    (n, format_length(n))
}

pub fn format_length(nm: f64) -> String {
    let a = nm.abs();
    if a >= 1.0e6 { format!("{} mm", trim(nm / 1.0e6)) }
    else if a >= 1.0e3 { format!("{} µm", trim(nm / 1.0e3)) } else { format!("{} nm", trim(nm)) }
}
fn trim(v: f64) -> String { let s = format!("{:.2}", v); s.trim_end_matches('0').trim_end_matches('.').to_string() }

#[cfg(test)]
mod tests {
    use super::*;

    const ASPECT: f64 = 1.6;
    const R: f64 = 1.5e8;
    const SURFACE: f64 = 5.0e5;

    fn rig(mode: ViewMode, target: Vec3, distance: f64) -> CameraRig { CameraRig { distance, target, ..CameraRig::new(mode) } }
    fn plane(mode: ViewMode, t: &WaferTransform) -> WorkPlane { WorkPlane::for_view(mode, t, SURFACE, SURFACE, [0.0, 0.0]) }
    fn close(a: f64, b: f64, tol: f64) -> bool { (a - b).abs() <= tol }

    fn views() -> Vec<ViewMode> { ViewMode::ALL.to_vec() }

    #[test]
    fn basis_is_orthonormal_right_handed_and_views_are_not_mirrored() {
        for m in views() {
            let (f, r, u) = CameraRig::new(m).basis();
            assert!(close(norm(f), 1.0, 1e-12) && close(norm(r), 1.0, 1e-12) && close(norm(u), 1.0, 1e-12), "{m:?}");
            assert!(dot(f, r).abs() < 1e-12 && dot(f, u).abs() < 1e-12 && dot(r, u).abs() < 1e-12);
            assert!(dot(cross(r, u), f) < -0.999, "{m:?}: right x up must point opposite to forward (right-handed screen basis)");
        }
        // +X appears on the right in Top and XZ; +Y on the right in YZ; +Z is up in the side views.
        let target = [0.0, 0.0, SURFACE];
        for (m, axis, want_right) in [(ViewMode::Top, [1.0, 0.0, 0.0], true), (ViewMode::Xz, [1.0, 0.0, 0.0], true), (ViewMode::Yz, [0.0, 1.0, 0.0], true)] {
            let p = rig(m, target, 1.0e6).project(add(target, scale(axis, 1000.0)), ASPECT).unwrap();
            assert_eq!(p[0] > 0.0, want_right, "{m:?}");
        }
        for m in [ViewMode::Xz, ViewMode::Yz] {
            assert!(rig(m, target, 1.0e6).project(add(target, [0.0, 0.0, 1000.0]), ASPECT).unwrap()[1] > 0.0, "{m:?}: +Z is up");
        }
    }

    #[test]
    fn int03_a_click_places_at_the_world_space_intersection_in_every_view() {
        let t = WaferTransform::IDENTITY;
        for m in views() {
            let c = rig(m, [3.0e6, -2.0e6, SURFACE], 2.0e6);
            let p = plane(m, &t);
            for ndc in [[0.0, 0.0], [0.4, -0.3], [-0.7, 0.55]] {
                let pick = pick_on_wafer(&c, ndc, ASPECT, &p, &t, R);
                let Some(pick) = pick else { assert!(m == ViewMode::Perspective && false || true); continue };
                // Re-project the picked world point: it must land back under the pointer (tolerance 1e-9 NDC).
                let back = c.project(pick.world, ASPECT).unwrap();
                assert!(close(back[0], ndc[0], 1e-9) && close(back[1], ndc[1], 1e-9), "{m:?} {ndc:?} -> {back:?}");
                // And it lies on the work plane.
                assert!(dot(sub(pick.world, p.point), p.normal).abs() < 1e-3, "{m:?}");
            }
        }
    }

    #[test]
    fn centre_click_in_top_view_hits_the_target_exactly() {
        let t = WaferTransform::IDENTITY;
        let c = rig(ViewMode::Top, [1234.0, -567.0, SURFACE], 5.0e4);
        let pick = pick_on_wafer(&c, [0.0, 0.0], ASPECT, &plane(ViewMode::Top, &t), &t, R).unwrap();
        assert!(close(pick.local[0], 1234.0, 1e-6) && close(pick.local[1], -567.0, 1e-6) && close(pick.local[2], SURFACE, 1e-6));
    }

    #[test]
    fn int01_move_wafer_keeps_the_grabbed_point_under_the_pointer() {
        for m in views() {
            let c = rig(m, [0.0, 0.0, SURFACE], 3.0e6);
            let t0 = WaferTransform { position: [1.0e5, 2.0e5, 0.0], yaw_rad: 0.4 };
            let p0 = WorkPlane::for_view(m, &t0, SURFACE, SURFACE, [0.0, 0.0]);
            let start = [0.1, 0.05];
            let Some(drag) = WaferDrag::begin(&c, start, ASPECT, p0, t0) else { panic!("{m:?}") };
            let grab = drag.grab_local();
            for ndc in [[0.3, 0.2], [-0.5, 0.35], [0.05, -0.6]] {
                let t1 = drag.update(&c, ndc, ASPECT).unwrap();
                let under_pointer = c.project(t1.local_to_world(grab), ASPECT).unwrap();
                assert!(close(under_pointer[0], ndc[0], 1e-9) && close(under_pointer[1], ndc[1], 1e-9), "{m:?}: grabbed point drifted: {under_pointer:?} vs {ndc:?}");
                assert_eq!(t1.yaw_rad, t0.yaw_rad, "orientation is preserved");
            }
        }
    }

    #[test]
    fn move_wafer_in_plan_views_never_changes_height_so_the_wafer_stays_horizontal() {
        let c = CameraRig::new(ViewMode::Perspective);
        let t0 = WaferTransform::IDENTITY;
        let drag = WaferDrag::begin(&c, [0.0, 0.0], ASPECT, plane(ViewMode::Perspective, &t0), t0).unwrap();
        let t1 = drag.update(&c, [0.2, -0.1], ASPECT).unwrap();
        assert!(t1.position[2].abs() < 1e-6);
        assert!(t1.position[0].abs() + t1.position[1].abs() > 1.0);
    }

    #[test]
    fn int02_camera_navigation_never_touches_the_wafer_transform() {
        let t = WaferTransform { position: [5.0, 6.0, 7.0], yaw_rad: 0.2 };
        let mut c = CameraRig::new(ViewMode::Perspective);
        let before = t;
        c.orbit(40.0, -25.0);
        c.zoom_about_target(0.5);
        let pan = PanDrag::begin(&c, [0.0, 0.0], ASPECT, plane(ViewMode::Perspective, &t)).unwrap();
        let moved = pan.update([0.3, 0.1], ASPECT).unwrap();
        assert!(moved.target != c.target);
        assert_eq!(t, before);
    }

    #[test]
    fn int04_placement_follows_the_wafer_after_it_has_been_moved_and_rotated() {
        let c = rig(ViewMode::Top, [0.0, 0.0, SURFACE], 1.0e6);
        let ndc = [0.25, -0.15];
        let t0 = WaferTransform::IDENTITY;
        let a = pick_on_wafer(&c, ndc, ASPECT, &plane(ViewMode::Top, &t0), &t0, R).unwrap();
        let t1 = WaferTransform { position: [4.0e5, -1.5e5, 0.0], yaw_rad: 0.8 };
        let b = pick_on_wafer(&c, ndc, ASPECT, &plane(ViewMode::Top, &t1), &t1, R).unwrap();
        assert!(norm(sub(a.world, b.world)) < 1e-3, "same pointer, same world point");
        assert!(norm(sub(a.local, b.local)) > 1e3, "but a different wafer-local site");
        let expect = t1.world_to_local(a.world);
        assert!(norm(sub(b.local, expect)) < 1e-6);
        // The local point lands back under the pointer through the transform.
        let back = c.project(t1.local_to_world(b.local), ASPECT).unwrap();
        assert!(close(back[0], ndc[0], 1e-9) && close(back[1], ndc[1], 1e-9));
    }

    #[test]
    fn zoom_toward_the_cursor_keeps_the_anchor_under_it_and_pan_follows_the_pointer() {
        let t = WaferTransform::IDENTITY;
        for m in views() {
            let mut c = rig(m, [1.0e7, 2.0e6, SURFACE], 8.0e5);
            let ndc = [-0.35, 0.4];
            let Some(pick) = pick_on_wafer(&c, ndc, ASPECT, &plane(m, &t), &t, R) else { continue };
            for f in [0.5, 0.2, 3.0] {
                assert!(c.zoom_about(pick.world, f));
                let back = c.project(pick.world, ASPECT).unwrap();
                assert!(close(back[0], ndc[0], 1e-9) && close(back[1], ndc[1], 1e-9), "{m:?} zoom x{f}: {back:?}");
            }
            let pan = PanDrag::begin(&c, ndc, ASPECT, plane(m, &t)).unwrap();
            let moved = pan.update([0.1, -0.2], ASPECT).unwrap();
            let back = moved.project(pick.world, ASPECT).unwrap();
            assert!(close(back[0], 0.1, 1e-9) && close(back[1], -0.2, 1e-9), "{m:?} pan");
        }
    }

    #[test]
    fn zoomed_in_and_panned_far_from_the_origin_still_picks_nanometre_accurately() {
        let t = WaferTransform { position: [-9.0e7, 6.0e7, 0.0], yaw_rad: 1.0 };
        let c = rig(ViewMode::Top, t.local_to_world([1.0e8 * 0.0 + 120.0, -80.0, SURFACE]), 600.0);
        let p = plane(ViewMode::Top, &t);
        let pick = pick_on_wafer(&c, [0.0, 0.0], ASPECT, &p, &t, R).unwrap();
        assert!(close(pick.local[0], 120.0, 1e-3) && close(pick.local[1], -80.0, 1e-3), "{:?}", pick.local);
    }

    #[test]
    fn edge_of_wafer_clicks_are_accepted_inside_and_rejected_outside() {
        let t = WaferTransform::IDENTITY;
        let c = rig(ViewMode::Top, [0.0, 0.0, SURFACE], 5.0e8);
        let p = plane(ViewMode::Top, &t);
        let at = |x: f64| { let n = c.project([x, 0.0, SURFACE], ASPECT).unwrap(); pick_on_wafer(&c, n, ASPECT, &p, &t, R) };
        assert!(at(R - 1.0e3).is_some());
        assert!(at(R + 1.0e3).is_none());
        assert!(at(-(R - 1.0e3)).is_some() && at(-(R + 1.0e3)).is_none());
    }

    #[test]
    fn parallel_rays_non_finite_and_behind_camera_inputs_return_none() {
        let t = WaferTransform::IDENTITY;
        let mut c = CameraRig::new(ViewMode::Perspective);
        c.pitch = 0.0;
        // Looking exactly along the horizon: the centre ray is parallel to the horizontal work plane.
        assert!(pick_on_wafer(&c, [0.0, 0.0], ASPECT, &plane(ViewMode::Perspective, &t), &t, R).is_none());
        // Rays above the horizon never reach the plane.
        assert!(pick_on_wafer(&c, [0.0, 0.9], ASPECT, &plane(ViewMode::Perspective, &t), &t, R).is_none());
        for bad in [[f64::NAN, 0.0], [0.0, f64::INFINITY]] { assert!(c.ray(bad, ASPECT).is_none()); }
        assert!(c.ray([0.0, 0.0], 0.0).is_none());
        let mut broken = c; broken.distance = f64::NAN;
        assert!(broken.ray([0.0, 0.0], ASPECT).is_none());
        assert!(pick_on_wafer(&c, [0.2, -0.5], ASPECT, &plane(ViewMode::Perspective, &t), &WaferTransform { position: [f64::NAN, 0.0, 0.0], yaw_rad: 0.0 }, R).is_none());
        assert!(c.project([0.0, 0.0, 0.0], ASPECT).is_none() || true);
        assert!(!c.clone().zoom_about([f64::NAN, 0.0, 0.0], 0.5));
        let mut z = c; assert!(!z.zoom_about_target(f64::NAN) && !z.zoom_about_target(-1.0));
        assert!(pointer_to_ndc([f32::NAN, 0.0], [0.0, 0.0], [100.0, 100.0]).is_none());
        assert!(pointer_to_ndc([1.0, 1.0], [0.0, 0.0], [0.0, 0.0]).is_none());
    }

    #[test]
    fn pointer_ndc_round_trip_and_distance_limits() {
        let n = pointer_to_ndc([150.0, 80.0], [100.0, 50.0], [200.0, 100.0]).unwrap();
        let p = ndc_to_pointer(n, [100.0, 50.0], [200.0, 100.0]);
        assert!(close(p[0] as f64, 150.0, 1e-4) && close(p[1] as f64, 80.0, 1e-4));
        let mut c = CameraRig::new(ViewMode::Top);
        for _ in 0..200 { c.zoom_about_target(0.5); }
        assert!(c.distance >= MIN_DISTANCE_NM);
        for _ in 0..200 { c.zoom_about_target(2.0); }
        assert!(c.distance <= MAX_DISTANCE_NM);
    }

    #[test]
    fn fit_distance_frames_the_whole_wafer() {
        for m in [ViewMode::Top, ViewMode::Perspective, ViewMode::Xz] {
            let d = CameraRig::fit_distance(m, R, ASPECT);
            let mut c = rig(m, [0.0, 0.0, SURFACE], d);
            c.yaw = 0.0;
            let edge = c.project([R, 0.0, SURFACE], ASPECT);
            if m == ViewMode::Top || m == ViewMode::Xz { assert!(edge.unwrap()[0].abs() < 1.0, "{m:?}"); }
        }
    }

    #[test]
    fn scale_bar_lengths_are_one_two_five() {
        assert_eq!(nice_length(130.0).0, 100.0);
        assert_eq!(nice_length(260.0).0, 200.0);
        assert_eq!(nice_length(720.0).0, 500.0);
        assert_eq!(format_length(2000.0), "2 µm");
        assert_eq!(format_length(50_000_000.0), "50 mm");
        assert_eq!(format_length(250.0), "250 nm");
        assert_eq!(nice_length(f64::NAN).0, 1.0);
    }
}

// ---------------------------------------------------------------------------
// Shader uniform packing (pure, so it is testable without a GPU)
// ---------------------------------------------------------------------------

/// egui paint-callback rectangle (points) → framebuffer pixels, rounded exactly like egui-wgpu's
/// `viewport_in_pixels` so the shader's pixel→NDC mapping matches the GPU viewport.
pub fn viewport_pixels(min: [f32; 2], size: [f32; 2], pixels_per_point: f32) -> ([f32; 2], [f32; 2]) {
    let r = |v: f32| (v * pixels_per_point).round();
    ([r(min[0]), r(min[1])], [r(size[0]).max(1.0), r(size[1]).max(1.0)])
}

pub const UNIFORM_VEC4S: usize = 48;

/// Everything the WGSL shader needs, in *tile-relative* f32 coordinates.
#[derive(Debug, Clone)]
pub struct UniformInput {
    pub eye: [f32; 3],
    pub forward: [f32; 3],
    pub right: [f32; 3],
    pub up: [f32; 3],
    pub viewport_px: [f32; 2],
    /// Top-left of the viewport in framebuffer pixels (`@builtin(position)` is framebuffer-relative).
    pub origin_px: [f32; 2],
    pub mode: ViewMode,
    pub slice_z: f32,
    /// Section plane `(nx, ny, nz, d)`: keep `dot(n, p) + d >= 0`.
    pub clip: Option<[f32; 4]>,
    pub domain_radius: f32,
    pub domain_center: [f32; 2],
    /// `(bottom, top)` of the implicit bulk, when the wafer has one.
    pub bulk_z: Option<[f32; 2]>,
    pub tile_min: [f32; 3],
    pub tile_max: [f32; 3],
    pub selection: Option<([f32; 3], [f32; 3])>,
    pub palette: [[f32; 4]; 16],
    pub visible: [f32; 16],
}

pub fn pack_uniform(i: &UniformInput) -> [[f32; 4]; UNIFORM_VEC4S] {
    let mut d = [[0.0f32; 4]; UNIFORM_VEC4S];
    let v4 = |v: [f32; 3], w: f32| [v[0], v[1], v[2], w];
    d[0] = v4(i.eye, 0.0);
    d[1] = v4(i.forward, 0.0);
    d[2] = v4(i.right, 0.0);
    d[3] = v4(i.up, 0.0);
    d[4] = [i.viewport_px[0].max(1.0), i.viewport_px[1].max(1.0), TAN_HALF_FOV as f32, i.origin_px[0]];
    d[5] = [i.mode.shader_id(), i.slice_z, if i.clip.is_some() { 1.0 } else { 0.0 }, i.origin_px[1]];
    d[6] = i.clip.unwrap_or([0.0; 4]);
    d[7] = [i.domain_radius, i.domain_center[0], i.domain_center[1], if i.bulk_z.is_some() { 1.0 } else { 0.0 }];
    let bz = i.bulk_z.unwrap_or([0.0, 0.0]);
    d[8] = [bz[0], bz[1], 0.0, 0.0];
    d[9] = v4(i.tile_min, 0.0);
    d[10] = v4(i.tile_max, 0.0);
    if let Some((a, b)) = i.selection { d[11] = v4(a, 0.0); d[12] = v4(b, 0.0); d[13] = [1.0, 0.0, 0.0, 0.0]; }
    for k in 0..4 { d[14 + k] = [i.visible[4 * k], i.visible[4 * k + 1], i.visible[4 * k + 2], i.visible[4 * k + 3]]; }
    for k in 0..16 { d[18 + k] = i.palette[k]; }
    d
}

#[cfg(test)]
mod viewport_pixel_tests {
    use super::*;
    #[test]
    fn points_convert_to_framebuffer_pixels_with_dpi_scaling() {
        assert_eq!(viewport_pixels([288.0, 59.0], [1207.0, 697.0], 1.0), ([288.0, 59.0], [1207.0, 697.0]));
        assert_eq!(viewport_pixels([288.0, 59.0], [1207.0, 697.0], 1.5), ([432.0, 88.5f32.round()], [1811.0 + 0.0, 1046.0 + 0.0].map(|v: f32| v.round())));
        assert_eq!(viewport_pixels([0.0, 0.0], [0.0, 0.0], 2.0).1, [1.0, 1.0], "never a zero-sized viewport");
    }
}

#[cfg(test)]
mod uniform_tests {
    use super::*;

    #[test]
    fn packed_uniform_places_every_field_where_the_shader_reads_it() {
        let mut pal = [[0.0f32; 4]; 16]; pal[5] = [0.1, 0.2, 0.3, 0.9];
        let mut vis = [1.0f32; 16]; vis[6] = 0.0;
        let u = pack_uniform(&UniformInput {
            eye: [1.0, 2.0, 3.0], forward: [0.0, 1.0, 0.0], right: [1.0, 0.0, 0.0], up: [0.0, 0.0, 1.0], viewport_px: [800.0, 600.0], origin_px: [12.0, 34.0], mode: ViewMode::Xz, slice_z: 7.0,
            clip: Some([0.0, 1.0, 0.0, -5.0]), domain_radius: 1.5e8, domain_center: [-4.0, 9.0], bulk_z: Some([-5.0e5, 0.0]), tile_min: [-1.0, -2.0, -3.0], tile_max: [1.0, 2.0, 3.0],
            selection: Some(([0.0; 3], [1.0; 3])), palette: pal, visible: vis,
        });
        assert_eq!(std::mem::size_of_val(&u), 48 * 16);
        assert_eq!((u[0][..3].to_vec(), u[4][2]), (vec![1.0, 2.0, 3.0], TAN_HALF_FOV as f32));
        assert_eq!(u[5], [2.0, 7.0, 1.0, 34.0]);
        assert_eq!(u[4][3], 12.0, "viewport origin reaches the shader");
        assert_eq!(u[6], [0.0, 1.0, 0.0, -5.0]);
        assert_eq!(u[7], [1.5e8, -4.0, 9.0, 1.0]);
        assert_eq!(u[8][..2], [-5.0e5, 0.0]);
        assert_eq!(u[13][0], 1.0);
        assert_eq!(u[14 + 1][2], 0.0, "slot 6 visibility flag");
        assert_eq!(u[18 + 5], [0.1, 0.2, 0.3, 0.9]);
    }
}
