// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Three-dimensional particle-flux solve: shadowing, sticking and re-emission (PRD roadmap 6).
//!
//! `flux.rs` answers one question in two dimensions: can this column see the source at all?
//! That is shadowing, and it is the dominant term, but it stops there. A real precursor that
//! does not stick on first contact bounces, and the bounce is what coats a sidewall and
//! lifts the bottom of a deep feature above its line-of-sight value. Bounces go sideways, so
//! the geometry that matters is a round hole or a corner, which a cross-section cannot hold.
//!
//! This module follows particles through a voxel grid in three dimensions:
//!
//! * **Source.** Rays are launched uniformly over the top plane. Either a cosine-law
//!   (Lambertian) gas-phase source, or a beam with a tilt, an azimuth and a Gaussian spread.
//! * **Transport.** Each ray is walked voxel by voxel (Amanatides and Woo) until it strikes
//!   solid or leaves through the top. The lateral boundaries are periodic, so a feature
//!   stands for an infinite array of itself.
//! * **Surface.** At every strike a fraction `sticking` of the arriving weight is absorbed;
//!   the rest is re-emitted diffusely (Lambert's law about the surface normal).
//! * **Output.** The impact flux on every exposed voxel face, relative to the flux on an open
//!   horizontal surface (1.0), and a weight budget that must balance.
//!
//! # What is and is not modelled
//!
//! Modelled: line-of-sight shadowing from any azimuth, finite angular spread, diffuse
//! re-emission, constant sticking, periodic arrays, closed substrate floor.
//!
//! Not modelled: surface diffusion, sticking that depends on angle or coverage, gas-phase
//! collisions (Knudsen number is taken as infinite), and surface motion. This solves the
//! flux field for a *given* geometry; it does not advance the surface. A deposition step that
//! uses it re-solves as the film grows.
//!
//! # Determinism (NFR-1)
//!
//! Every ray owns a random stream derived only from `(seed, ray index)`. Rays are split into a
//! fixed number of chunks and the chunk tallies are summed in chunk order, so the result is
//! bit-identical however the chunks were scheduled. This crate starts no threads (PRD A-02):
//! [`Plan`] hands out independent chunks, and a caller that wants parallelism runs
//! [`Plan::run_chunk`] on its own workers and gives the parts back to [`Plan::finish`].
//!
//! # Cost (NFR-2, NFR-3)
//!
//! The grid is capped at [`MAX_VOXELS`] and the ray count at [`MAX_RAYS`]; a run past either is
//! refused before allocating. A solve is a batch computation: run it on a worker thread, never
//! on the interface thread.
//!
//! This module is a library: the `field3d` binary in `atto_cli` is its only consumer so far.
//! It is not yet called by the deposition stage (PRD G-5).

use atto_units::Nm;
use std::fmt;
use std::ops::Range;

/// Largest grid accepted, in voxels. At this size the face table is about 50 MB.
pub const MAX_VOXELS: usize = 1 << 21;

/// Largest ray count accepted in one solve.
pub const MAX_RAYS: u64 = 200_000_000;

/// Number of fixed chunks the rays are divided into. Fixed so the result never depends on how
/// the chunks were scheduled.
const CHUNKS: u64 = 16;

/// A ray whose remaining weight falls below this is dropped and counted as truncated.
const MIN_WEIGHT: f64 = 1.0e-4;

/// A ray that bounces more than this many times is dropped and counted as truncated.
const MAX_BOUNCES: u32 = 2000;

/// The six voxel faces, in the order used by every table here.
/// 0 is -x, 1 is +x, 2 is -y, 3 is +y, 4 is -z, 5 is +z.
pub const FACE_TOP: u8 = 5;

/// Why a grid, a source or a run was refused.
#[derive(Debug, Clone, PartialEq)]
pub enum FieldError {
    /// The voxel grid is malformed or too large.
    BadGrid(String),
    /// The source description is outside its domain.
    BadSource(String),
    /// The run configuration is outside its domain.
    BadRun(String),
}

impl fmt::Display for FieldError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FieldError::BadGrid(m) => write!(f, "voxel grid refused: {m}"),
            FieldError::BadSource(m) => write!(f, "source refused: {m}"),
            FieldError::BadRun(m) => write!(f, "run refused: {m}"),
        }
    }
}

impl std::error::Error for FieldError {}

// ---------------------------------------------------------------------------
// Geometry
// ---------------------------------------------------------------------------

/// A box of cubic voxels, solid or empty. Axes are x, y, z with z up. The lateral boundaries are
/// periodic; the bottom layer is solid (a closed substrate); the top layer is empty (headroom
/// for the source plane).
#[derive(Debug, Clone, PartialEq)]
pub struct VoxelGrid {
    nx: usize,
    ny: usize,
    nz: usize,
    cell: Nm,
    solid: Vec<bool>,
}

impl VoxelGrid {
    /// Build a grid from an explicit solid map, laid out as `((iz * ny) + iy) * nx + ix`.
    ///
    /// Refused: any dimension zero, fewer than two layers, more than [`MAX_VOXELS`] voxels, a
    /// non-positive or non-finite cell, a map of the wrong length, a bottom layer that is not
    /// entirely solid, or a top layer that is not entirely empty.
    pub fn new(nx: usize, ny: usize, nz: usize, cell: Nm, solid: Vec<bool>) -> Result<Self, FieldError> {
        if nx == 0 || ny == 0 || nz < 2 {
            return Err(FieldError::BadGrid(format!("dimensions {nx} x {ny} x {nz} need nx, ny >= 1 and nz >= 2")));
        }
        let total = nx
            .checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .ok_or_else(|| FieldError::BadGrid("the voxel count overflows".into()))?;
        if total > MAX_VOXELS {
            return Err(FieldError::BadGrid(format!("{total} voxels exceeds the cap of {MAX_VOXELS}")));
        }
        if !cell.is_finite() || cell.value() <= 0.0 {
            return Err(FieldError::BadGrid("the cell size must be finite and positive".into()));
        }
        if solid.len() != total {
            return Err(FieldError::BadGrid(format!("the solid map has {} entries for {total} voxels", solid.len())));
        }
        let layer = nx * ny;
        if solid[..layer].iter().any(|s| !*s) {
            return Err(FieldError::BadGrid("the bottom layer must be entirely solid".into()));
        }
        if solid[total - layer..].iter().any(|s| *s) {
            return Err(FieldError::BadGrid("the top layer must be entirely empty: it is the source plane's headroom".into()));
        }
        Ok(VoxelGrid { nx, ny, nz, cell, solid })
    }

    /// Build a grid from one surface height per column, `ny` rows of `nx` heights.
    ///
    /// Each height is rounded to a whole number of cells and everything below it is solid. A
    /// height below one cell is raised to one (the substrate floor); a height that leaves no
    /// empty top layer is refused.
    pub fn from_heights(nx: usize, ny: usize, nz: usize, cell: Nm, heights: &[Nm]) -> Result<Self, FieldError> {
        if nx == 0 || ny == 0 || nz < 2 {
            return Err(FieldError::BadGrid(format!("dimensions {nx} x {ny} x {nz} need nx, ny >= 1 and nz >= 2")));
        }
        if !cell.is_finite() || cell.value() <= 0.0 {
            return Err(FieldError::BadGrid("the cell size must be finite and positive".into()));
        }
        if heights.len() != nx * ny {
            return Err(FieldError::BadGrid(format!("{} heights for {} columns", heights.len(), nx * ny)));
        }
        let total = nx
            .checked_mul(ny)
            .and_then(|v| v.checked_mul(nz))
            .ok_or_else(|| FieldError::BadGrid("the voxel count overflows".into()))?;
        if total > MAX_VOXELS {
            return Err(FieldError::BadGrid(format!("{total} voxels exceeds the cap of {MAX_VOXELS}")));
        }
        let mut solid = vec![false; total];
        for iy in 0..ny {
            for ix in 0..nx {
                let h = heights[iy * nx + ix];
                if !h.is_finite() || h.value() < 0.0 {
                    return Err(FieldError::BadGrid(format!("the height at column ({ix}, {iy}) is {h}")));
                }
                let layers = ((h.value() / cell.value()).round() as usize).max(1);
                if layers > nz - 1 {
                    return Err(FieldError::BadGrid(format!(
                        "the height at column ({ix}, {iy}) needs {layers} layers but only {} fit below the headroom layer",
                        nz - 1
                    )));
                }
                for iz in 0..layers {
                    solid[(iz * ny + iy) * nx + ix] = true;
                }
            }
        }
        VoxelGrid::new(nx, ny, nz, cell, solid)
    }

    /// A trench running the full length of y, centred in x. `floor` cells of solid lie below
    /// the trench bottom and `depth` cells deep it is cut; two or more empty layers are left on
    /// top. With periodic boundaries this is an infinite array of identical trenches.
    pub fn trench(nx: usize, ny: usize, nz: usize, cell: Nm, width: usize, depth: usize, floor: usize) -> Result<Self, FieldError> {
        Self::cut(nx, ny, nz, cell, width, width, depth, floor, false)
    }

    /// A square hole `width` cells on a side, centred in x and y, otherwise as [`VoxelGrid::trench`].
    pub fn hole(nx: usize, ny: usize, nz: usize, cell: Nm, width: usize, depth: usize, floor: usize) -> Result<Self, FieldError> {
        Self::cut(nx, ny, nz, cell, width, width, depth, floor, true)
    }

    fn cut(nx: usize, ny: usize, nz: usize, cell: Nm, wx: usize, wy: usize, depth: usize, floor: usize, bounded_y: bool) -> Result<Self, FieldError> {
        if wx == 0 || depth == 0 || floor == 0 || wx >= nx || (bounded_y && wy >= ny) {
            return Err(FieldError::BadGrid("a cut needs 0 < width < the grid, depth >= 1 and floor >= 1".into()));
        }
        if nz < floor + depth + 2 {
            return Err(FieldError::BadGrid(format!("nz = {nz} leaves no headroom above a surface at {} cells", floor + depth)));
        }
        let (x0, y0) = ((nx - wx) / 2, (ny.saturating_sub(wy)) / 2);
        let mut heights = vec![Nm::new((floor + depth) as f64 * cell.value()); nx * ny];
        for iy in 0..ny {
            for ix in x0..x0 + wx {
                let inside_y = !bounded_y || (iy >= y0 && iy < y0 + wy);
                if inside_y {
                    heights[iy * nx + ix] = Nm::new(floor as f64 * cell.value());
                }
            }
        }
        VoxelGrid::from_heights(nx, ny, nz, cell, &heights)
    }

    /// The grid dimensions in voxels: x, y, z.
    pub fn dims(&self) -> (usize, usize, usize) {
        (self.nx, self.ny, self.nz)
    }

    /// The edge length of one voxel.
    pub fn cell(&self) -> Nm {
        self.cell
    }

    /// Whether the voxel at `(ix, iy, iz)` is solid. Out-of-range indices are empty.
    pub fn is_solid(&self, ix: usize, iy: usize, iz: usize) -> bool {
        ix < self.nx && iy < self.ny && iz < self.nz && self.solid[(iz * self.ny + iy) * self.nx + ix]
    }

    /// Solid test with periodic wrapping in x and y. The caller guarantees `0 <= iz < nz`.
    fn solid_wrapped(&self, ix: i64, iy: i64, iz: i64) -> bool {
        let x = ix.rem_euclid(self.nx as i64) as usize;
        let y = iy.rem_euclid(self.ny as i64) as usize;
        self.solid[(iz as usize * self.ny + y) * self.nx + x]
    }

    /// Index of the highest solid voxel in each column, `ny` rows of `nx`.
    fn column_tops(&self) -> Vec<usize> {
        let mut tops = vec![0usize; self.nx * self.ny];
        for iy in 0..self.ny {
            for ix in 0..self.nx {
                let mut iz = self.nz - 1;
                while iz > 0 && !self.solid[(iz * self.ny + iy) * self.nx + ix] {
                    iz -= 1;
                }
                tops[iy * self.nx + ix] = iz;
            }
        }
        tops
    }

    /// Give every exposed face (a solid voxel face with an empty neighbour) a slot number.
    /// Returns the table, `u32::MAX` meaning "not exposed", and the number of slots.
    fn exposed_faces(&self) -> (Vec<u32>, usize) {
        let mut slot = vec![u32::MAX; self.solid.len() * 6];
        let mut n = 0u32;
        for iz in 0..self.nz {
            for iy in 0..self.ny {
                for ix in 0..self.nx {
                    let v = (iz * self.ny + iy) * self.nx + ix;
                    if !self.solid[v] {
                        continue;
                    }
                    for face in 0..6u8 {
                        let axis = (face / 2) as usize;
                        let delta: i64 = if face % 2 == 0 { -1 } else { 1 };
                        let mut nb = [ix as i64, iy as i64, iz as i64];
                        nb[axis] += delta;
                        // below the floor there is nothing to be exposed to; above the top is
                        // impossible because the top layer is empty
                        if nb[2] < 0 || nb[2] >= self.nz as i64 {
                            continue;
                        }
                        if !self.solid_wrapped(nb[0], nb[1], nb[2]) {
                            slot[v * 6 + face as usize] = n;
                            n += 1;
                        }
                    }
                }
            }
        }
        (slot, n as usize)
    }
}

// ---------------------------------------------------------------------------
// Source and run configuration
// ---------------------------------------------------------------------------

/// The angular distribution of the arriving species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Angular {
    /// Cosine-law (Lambertian) crossing of the source plane: a gas-phase CVD source with no
    /// preferred direction. Tilt, azimuth and spread are ignored.
    Cosine,
    /// A beam about a central direction with a Gaussian angular spread: a sputter or
    /// evaporation source.
    Beam,
}

/// The arriving species.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Source3d {
    /// The angular distribution.
    pub angular: Angular,
    /// Beam tilt from vertical, degrees, in `[0, 90)`. Beam only.
    pub tilt_deg: f64,
    /// Direction the beam leans toward, degrees from +x toward +y. Beam only.
    pub azimuth_deg: f64,
    /// One-sigma angular spread about the central direction, degrees, in `[0, 45]`. Beam only.
    pub spread_deg: f64,
}

impl Default for Source3d {
    fn default() -> Self {
        Source3d { angular: Angular::Cosine, tilt_deg: 0.0, azimuth_deg: 0.0, spread_deg: 0.0 }
    }
}

impl Source3d {
    /// Refuse a source outside its domain.
    pub fn validate(&self) -> Result<(), FieldError> {
        if !self.tilt_deg.is_finite() || !self.azimuth_deg.is_finite() || !self.spread_deg.is_finite() {
            return Err(FieldError::BadSource("tilt, azimuth and spread must be finite".into()));
        }
        if self.angular == Angular::Beam {
            if !(0.0..90.0).contains(&self.tilt_deg) {
                return Err(FieldError::BadSource(format!("a beam tilt of {} degrees is outside [0, 90)", self.tilt_deg)));
            }
            if !(0.0..=45.0).contains(&self.spread_deg) {
                return Err(FieldError::BadSource(format!("a spread of {} degrees is outside [0, 45]", self.spread_deg)));
            }
        }
        Ok(())
    }
}

/// How one solve is run.
#[derive(Debug, Clone, PartialEq)]
pub struct Config {
    /// Number of rays launched, `1..=MAX_RAYS`. The statistical error of a flux falls as the
    /// square root of the hits it collects.
    pub rays: u64,
    /// Seed of every random stream.
    pub seed: u64,
    /// Probability a striking particle is absorbed, in `(0, 1]`. 1 means no re-emission.
    pub sticking: f64,
    /// The arriving species.
    pub source: Source3d,
}

impl Default for Config {
    fn default() -> Self {
        Config { rays: 200_000, seed: 1, sticking: 1.0, source: Source3d::default() }
    }
}

impl Config {
    /// Refuse a configuration outside its domain.
    pub fn validate(&self) -> Result<(), FieldError> {
        if self.rays == 0 || self.rays > MAX_RAYS {
            return Err(FieldError::BadRun(format!("{} rays is outside 1..={MAX_RAYS}", self.rays)));
        }
        if !self.sticking.is_finite() || self.sticking <= 0.0 || self.sticking > 1.0 {
            return Err(FieldError::BadRun(format!("sticking {} is outside (0, 1]", self.sticking)));
        }
        self.source.validate()
    }
}

// ---------------------------------------------------------------------------
// Random numbers
// ---------------------------------------------------------------------------

/// A splitmix64 stream. Each ray derives its own from the seed and its index, so no ray's
/// path depends on any other's.
struct Rng(u64);

impl Rng {
    fn new(seed: u64, ray: u64) -> Self {
        let mut r = Rng(seed ^ ray.wrapping_mul(0x9E37_79B9_7F4A_7C15));
        r.next_u64(); // decorrelate neighbouring rays
        r
    }

    fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

fn normalised(a: [f64; 3]) -> [f64; 3] {
    let n = (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt();
    [a[0] / n, a[1] / n, a[2] / n]
}

/// A direction for a ray crossing the source plane downward.
fn sample_source(src: &Source3d, rng: &mut Rng) -> [f64; 3] {
    match src.angular {
        Angular::Cosine => {
            let cos_t = rng.unit().sqrt();
            let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
            let phi = std::f64::consts::TAU * rng.unit();
            [sin_t * phi.cos(), sin_t * phi.sin(), -cos_t]
        }
        Angular::Beam => {
            let (t, a) = (src.tilt_deg.to_radians(), src.azimuth_deg.to_radians());
            let centre = [t.sin() * a.cos(), t.sin() * a.sin(), -t.cos()];
            if src.spread_deg == 0.0 {
                return centre;
            }
            let reference = if centre[2].abs() > 0.999 { [1.0, 0.0, 0.0] } else { [0.0, 0.0, 1.0] };
            let e1 = normalised(cross(centre, reference));
            let e2 = cross(centre, e1);
            let sigma = src.spread_deg.to_radians();
            for _ in 0..8 {
                // a two-dimensional Gaussian in angle space: Rayleigh radius, uniform phase
                let delta = sigma * (-2.0 * (1.0 - rng.unit()).ln()).sqrt();
                let psi = std::f64::consts::TAU * rng.unit();
                let (sd, cd) = (delta.sin(), delta.cos());
                let d = [
                    centre[0] * cd + sd * (e1[0] * psi.cos() + e2[0] * psi.sin()),
                    centre[1] * cd + sd * (e1[1] * psi.cos() + e2[1] * psi.sin()),
                    centre[2] * cd + sd * (e1[2] * psi.cos() + e2[2] * psi.sin()),
                ];
                // a ray that would travel upward is not crossing the source plane: redraw
                if d[2] < -1.0e-3 {
                    return d;
                }
            }
            centre
        }
    }
}

/// A diffuse re-emission direction from a struck face. `face` is the face of the solid voxel that
/// was hit; the surface normal points out of the solid, toward the empty side.
fn lambert(face: u8, rng: &mut Rng) -> [f64; 3] {
    let axis = (face / 2) as usize;
    let sign = if face % 2 == 0 { -1.0 } else { 1.0 };
    let (b, c) = ((axis + 1) % 3, (axis + 2) % 3);
    let cos_t = rng.unit().sqrt();
    let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
    let phi = std::f64::consts::TAU * rng.unit();
    let mut d = [0.0; 3];
    d[axis] = sign * cos_t;
    d[b] = sin_t * phi.cos();
    d[c] = sin_t * phi.sin();
    d
}

// ---------------------------------------------------------------------------
// Transport
// ---------------------------------------------------------------------------

/// How a flight ended.
enum Flight {
    /// Struck solid. `voxel` is the solid voxel (unwrapped), `from` the empty voxel the ray was in,
    /// `face` the struck face of the solid voxel and `point` where on that face.
    Hit { voxel: [i64; 3], from: [i64; 3], face: u8, point: [f64; 3] },
    /// Left through the top.
    Escaped,
    /// Gave up: too many steps. Counted, never silently dropped.
    Lost,
}

/// The tally of one chunk of rays, produced by [`Plan::run_chunk`] and consumed by
/// [`Plan::finish`]. Opaque: it can be moved between threads but not inspected.
pub struct Partial {
    chunk: usize,
    impacts: Vec<f64>,
    launched: f64,
    absorbed: f64,
    escaped: f64,
    truncated: f64,
}

struct Tracer<'a> {
    grid: &'a VoxelGrid,
    slot: &'a [u32],
    max_steps: usize,
}

impl<'a> Tracer<'a> {
    /// Walk a ray from `pos`, which lies in or on the boundary of the empty voxel `start`, along
    /// `dir`, through the voxels until it strikes solid or leaves through the top.
    fn trace(&self, start: [i64; 3], pos: [f64; 3], dir: [f64; 3]) -> Flight {
        let nz = self.grid.nz as i64;
        let mut v = start;
        let mut step = [0i64; 3];
        let mut t_max = [f64::INFINITY; 3];
        let mut t_delta = [f64::INFINITY; 3];
        for a in 0..3 {
            if dir[a] > 0.0 {
                step[a] = 1;
                t_delta[a] = 1.0 / dir[a];
                t_max[a] = (((v[a] + 1) as f64 - pos[a]) / dir[a]).max(0.0);
            } else if dir[a] < 0.0 {
                step[a] = -1;
                t_delta[a] = -1.0 / dir[a];
                t_max[a] = ((pos[a] - v[a] as f64) / -dir[a]).max(0.0);
            }
        }
        for _ in 0..self.max_steps {
            let a = if t_max[0] <= t_max[1] && t_max[0] <= t_max[2] {
                0
            } else if t_max[1] <= t_max[2] {
                1
            } else {
                2
            };
            let t = t_max[a];
            let from = v;
            v[a] += step[a];
            t_max[a] += t_delta[a];
            if v[2] >= nz {
                return Flight::Escaped;
            }
            if v[2] < 0 {
                return Flight::Lost; // unreachable with a solid floor; kept honest anyway
            }
            if self.grid.solid_wrapped(v[0], v[1], v[2]) {
                // entering the solid voxel through the face that faces the direction of travel
                let face = if step[a] > 0 { 2 * a as u8 } else { 2 * a as u8 + 1 };
                let mut point = [pos[0] + dir[0] * t, pos[1] + dir[1] * t, pos[2] + dir[2] * t];
                // snap onto the face exactly so the restart never begins inside the solid
                point[a] = if step[a] > 0 { v[a] as f64 } else { (v[a] + 1) as f64 };
                return Flight::Hit { voxel: v, from, face, point };
            }
        }
        Flight::Lost
    }

    fn run_ray(&self, cfg: &Config, ray: u64, tally: &mut Partial) {
        let (nx, ny, nz) = self.grid.dims();
        let mut rng = Rng::new(cfg.seed, ray);
        let x = rng.unit() * nx as f64;
        let y = rng.unit() * ny as f64;
        let mut dir = sample_source(&cfg.source, &mut rng);
        let mut start = [(x.floor() as i64).min(nx as i64 - 1), (y.floor() as i64).min(ny as i64 - 1), nz as i64 - 1];
        let mut pos = [x, y, nz as f64];
        let mut weight = 1.0f64;
        let mut bounces = 0u32;
        tally.launched += 1.0;
        loop {
            match self.trace(start, pos, dir) {
                Flight::Escaped => {
                    tally.escaped += weight;
                    return;
                }
                Flight::Lost => {
                    tally.truncated += weight;
                    return;
                }
                Flight::Hit { voxel, from, face, point } => {
                    let vx = voxel[0].rem_euclid(nx as i64) as usize;
                    let vy = voxel[1].rem_euclid(ny as i64) as usize;
                    let vz = voxel[2] as usize;
                    let s = self.slot[((vz * ny + vy) * nx + vx) * 6 + face as usize];
                    // every face a ray can strike was entered from an empty neighbour, so it is exposed
                    debug_assert!(s != u32::MAX, "a ray struck a face that is not exposed");
                    if s != u32::MAX {
                        tally.impacts[s as usize] += weight;
                    }
                    let absorbed = weight * cfg.sticking;
                    tally.absorbed += absorbed;
                    weight -= absorbed;
                    bounces += 1;
                    if weight < MIN_WEIGHT || bounces > MAX_BOUNCES {
                        tally.truncated += weight;
                        return;
                    }
                    dir = lambert(face, &mut rng);
                    start = from;
                    pos = point;
                }
            }
        }
    }

    fn run_chunk(&self, cfg: &Config, chunk: u64, n_slots: usize) -> Partial {
        let first = cfg.rays * chunk / CHUNKS;
        let last = cfg.rays * (chunk + 1) / CHUNKS;
        let mut tally = Partial { chunk: chunk as usize, impacts: vec![0.0; n_slots], launched: 0.0, absorbed: 0.0, escaped: 0.0, truncated: 0.0 };
        for ray in first..last {
            self.run_ray(cfg, ray, &mut tally);
        }
        tally
    }
}

// ---------------------------------------------------------------------------
// Solution
// ---------------------------------------------------------------------------

/// Where every unit of launched weight went. The four add up to the launched weight; a solve in
/// which they do not is wrong, and [`Budget::imbalance`] says by how much.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Budget {
    /// Weight launched: one per ray.
    pub launched: f64,
    /// Weight absorbed by the surface.
    pub absorbed: f64,
    /// Weight that left through the top.
    pub escaped: f64,
    /// Weight dropped because a ray bounced too long or flew too far. Should be negligible.
    pub truncated: f64,
}

impl Budget {
    /// `|launched - absorbed - escaped - truncated| / launched`.
    pub fn imbalance(&self) -> f64 {
        if self.launched <= 0.0 {
            return 0.0;
        }
        (self.launched - self.absorbed - self.escaped - self.truncated).abs() / self.launched
    }

    /// The fraction of launched weight that was dropped.
    pub fn truncated_fraction(&self) -> f64 {
        if self.launched <= 0.0 {
            0.0
        } else {
            self.truncated / self.launched
        }
    }
}

/// The flux field of one solve.
#[derive(Debug, Clone, PartialEq)]
pub struct Solution {
    nx: usize,
    ny: usize,
    nz: usize,
    slot: Vec<u32>,
    tops: Vec<usize>,
    flux: Vec<f64>,
    sticking: f64,
    budget: Budget,
}

impl Solution {
    /// The weight budget of the solve.
    pub fn budget(&self) -> Budget {
        self.budget
    }

    /// The impact flux on every exposed face, in slot order, relative to an open horizontal surface.
    /// Two solves are identical exactly when these slices are.
    pub fn flux_values(&self) -> &[f64] {
        &self.flux
    }

    /// Relative impact flux on one face of one voxel, or `None` when that face is not exposed.
    /// Faces are numbered as in [`FACE_TOP`]'s documentation: 0 is -x ... 5 is +z.
    pub fn face_flux(&self, ix: usize, iy: usize, iz: usize, face: u8) -> Option<f64> {
        if ix >= self.nx || iy >= self.ny || iz >= self.nz || face > 5 {
            return None;
        }
        let s = self.slot[((iz * self.ny + iy) * self.nx + ix) * 6 + face as usize];
        if s == u32::MAX {
            None
        } else {
            Some(self.flux[s as usize])
        }
    }

    /// Relative impact flux on the upward face of the highest solid voxel in a column: what a film
    /// growing on that column's surface is fed by.
    pub fn top_flux(&self, ix: usize, iy: usize) -> Option<f64> {
        if ix >= self.nx || iy >= self.ny {
            return None;
        }
        self.face_flux(ix, iy, self.tops[iy * self.nx + ix], FACE_TOP)
    }

    /// Mean [`Solution::top_flux`] over a block of columns, or `None` when the block is empty or
    /// leaves the grid.
    pub fn mean_top_flux(&self, xs: Range<usize>, ys: Range<usize>) -> Option<f64> {
        if xs.is_empty() || ys.is_empty() || xs.end > self.nx || ys.end > self.ny {
            return None;
        }
        let mut sum = 0.0;
        let mut n = 0usize;
        for iy in ys {
            for ix in xs.clone() {
                sum += self.top_flux(ix, iy)?;
                n += 1;
            }
        }
        Some(sum / n as f64)
    }

    /// Relative deposition rate on a column's surface: the impact flux times the sticking
    /// probability. Re-emitted particles that arrive again are counted, which is what lifts a
    /// feature's bottom above its line-of-sight value.
    pub fn deposition_rate(&self, ix: usize, iy: usize) -> Option<f64> {
        self.top_flux(ix, iy).map(|f| f * self.sticking)
    }
}

/// A solve prepared for execution: the grid, its exposed-face table and the run configuration.
///
/// The rays are divided into [`Plan::chunk_count`] independent chunks. Run each with
/// [`Plan::run_chunk`], in any order and on any threads, then hand all the parts to
/// [`Plan::finish`]. The answer does not depend on the order or the schedule.
#[derive(Debug, Clone)]
pub struct Plan {
    grid: VoxelGrid,
    slot: Vec<u32>,
    n_slots: usize,
    tops: Vec<usize>,
    cfg: Config,
    max_steps: usize,
}

impl Plan {
    /// Validate `cfg` and prepare the face table. Refuses before any ray is traced.
    pub fn new(grid: &VoxelGrid, cfg: &Config) -> Result<Plan, FieldError> {
        cfg.validate()?;
        let (slot, n_slots) = grid.exposed_faces();
        let (nx, ny, nz) = grid.dims();
        Ok(Plan {
            grid: grid.clone(),
            slot,
            n_slots,
            tops: grid.column_tops(),
            cfg: cfg.clone(),
            max_steps: 64 * (nx + ny + nz) + 4096,
        })
    }

    /// How many chunks the rays are divided into. Fixed, whatever the configuration.
    pub fn chunk_count(&self) -> usize {
        CHUNKS as usize
    }

    /// Trace one chunk of rays. A `chunk` past [`Plan::chunk_count`] traces nothing.
    pub fn run_chunk(&self, chunk: usize) -> Partial {
        let tracer = Tracer { grid: &self.grid, slot: &self.slot, max_steps: self.max_steps };
        if chunk as u64 >= CHUNKS {
            return Partial { chunk, impacts: vec![0.0; self.n_slots], launched: 0.0, absorbed: 0.0, escaped: 0.0, truncated: 0.0 };
        }
        tracer.run_chunk(&self.cfg, chunk as u64, self.n_slots)
    }

    /// Combine every chunk's part into the flux field. Refused unless each chunk is present
    /// exactly once; the parts may arrive in any order.
    pub fn finish(self, mut parts: Vec<Partial>) -> Result<Solution, FieldError> {
        parts.sort_by_key(|p| p.chunk);
        if parts.len() != CHUNKS as usize || parts.iter().enumerate().any(|(i, p)| p.chunk != i) {
            return Err(FieldError::BadRun(format!("expected each of the {CHUNKS} chunks exactly once, got {} parts", parts.len())));
        }
        // summed in chunk order, whatever thread produced each chunk: the result is bit-identical
        let mut impacts = vec![0.0f64; self.n_slots];
        let mut budget = Budget { launched: 0.0, absorbed: 0.0, escaped: 0.0, truncated: 0.0 };
        for part in &parts {
            for (total, x) in impacts.iter_mut().zip(part.impacts.iter()) {
                *total += *x;
            }
            budget.launched += part.launched;
            budget.absorbed += part.absorbed;
            budget.escaped += part.escaped;
            budget.truncated += part.truncated;
        }
        let (nx, ny, nz) = self.grid.dims();
        // rays per unit area on the source plane is rays / (nx * ny) per face area, so this is the
        // arriving weight per face relative to the open-field value
        let scale = (nx * ny) as f64 / self.cfg.rays as f64;
        let flux: Vec<f64> = impacts.iter().map(|w| w * scale).collect();
        Ok(Solution { nx, ny, nz, slot: self.slot, tops: self.tops, flux, sticking: self.cfg.sticking, budget })
    }
}

/// Solve the flux field of `grid` under `cfg` on the calling thread.
///
/// Refuses an invalid configuration before doing any work. The result depends only on the grid
/// and `cfg`. To use several threads, drive a [`Plan`] directly.
pub fn solve(grid: &VoxelGrid, cfg: &Config) -> Result<Solution, FieldError> {
    let plan = Plan::new(grid, cfg)?;
    let parts: Vec<Partial> = (0..plan.chunk_count()).map(|c| plan.run_chunk(c)).collect();
    plan.finish(parts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn nm(v: f64) -> Nm {
        Nm::new(v)
    }

    fn cfg(rays: u64, sticking: f64) -> Config {
        Config { rays, sticking, ..Config::default() }
    }

    /// The flux on the bottom centre of a wide-open slot, from a Lambertian source, is the view
    /// factor of the slot opening: (w/2) / sqrt((w/2)^2 + d^2) for width w and depth d.
    fn trench_view_factor(w: f64, d: f64) -> f64 {
        (w / 2.0) / ((w / 2.0).powi(2) + d * d).sqrt()
    }

    /// View factor from a point on the axis to a square opening of side `w` at height `d`:
    /// four corner rectangles of half the side, each (1/2pi)[A/sqrt(1+A^2) atan(B/sqrt(1+A^2))
    /// + B/sqrt(1+B^2) atan(A/sqrt(1+B^2))] with A = B = w / (2d).
    fn square_view_factor(w: f64, d: f64) -> f64 {
        let a = w / (2.0 * d);
        let corner = (a / (1.0 + a * a).sqrt() * (a / (1.0 + a * a).sqrt()).atan()) * 2.0 / (2.0 * std::f64::consts::PI);
        4.0 * corner
    }

    #[test]
    fn an_open_flat_surface_receives_exactly_the_reference_flux() {
        let g = VoxelGrid::from_heights(16, 16, 16, nm(10.0), &vec![nm(40.0); 256]).unwrap();
        let s = solve(&g, &cfg(100_000, 0.5)).unwrap();
        let mean = s.mean_top_flux(0..16, 0..16).unwrap();
        assert!((mean - 1.0).abs() < 1e-9, "every ray lands on the open surface once: {mean}");
        let b = s.budget();
        assert!((b.absorbed - 0.5 * 100_000.0).abs() < 1e-6, "half of the arriving weight sticks: {}", b.absorbed);
    }

    #[test]
    fn weight_is_conserved_in_a_hole_with_re_emission() {
        let g = VoxelGrid::hole(24, 24, 40, nm(10.0), 10, 16, 4).unwrap();
        let s = solve(&g, &cfg(60_000, 0.3)).unwrap();
        let b = s.budget();
        assert!(b.imbalance() < 1e-9, "launched = absorbed + escaped + truncated, imbalance {}", b.imbalance());
        assert!(b.truncated_fraction() < 1e-2, "truncation must be negligible: {}", b.truncated_fraction());
        assert!(b.absorbed > 0.0 && b.escaped > 0.0);
    }

    #[test]
    fn a_trench_bottom_sees_its_view_factor_of_a_lambertian_source() {
        // 16 cells wide, 32 deep: aspect ratio 2, expected 8 / sqrt(64 + 1024) = 0.2425
        let g = VoxelGrid::trench(48, 4, 48, nm(10.0), 16, 32, 4).unwrap();
        let s = solve(&g, &cfg(400_000, 1.0)).unwrap();
        let got = s.mean_top_flux(23..25, 0..4).unwrap();
        let want = trench_view_factor(16.0, 32.0);
        assert!((got - want).abs() / want < 0.08, "bottom centre flux {got}, view factor {want}");
        let field = s.mean_top_flux(0..8, 0..4).unwrap();
        assert!((field - 1.0).abs() < 0.03, "the open field sees the whole source: {field}");
    }

    #[test]
    fn a_hole_bottom_is_darker_than_a_trench_bottom_of_the_same_aspect_ratio() {
        let hole = VoxelGrid::hole(32, 32, 40, nm(10.0), 16, 16, 4).unwrap();
        let trench = VoxelGrid::trench(32, 32, 40, nm(10.0), 16, 16, 4).unwrap();
        let h = solve(&hole, &cfg(600_000, 1.0)).unwrap().mean_top_flux(15..17, 15..17).unwrap();
        let t = solve(&trench, &cfg(600_000, 1.0)).unwrap().mean_top_flux(15..17, 0..32).unwrap();
        assert!(h < t, "a hole closes in on the bottom from four sides, a trench from two: {h} vs {t}");
        let want = square_view_factor(16.0, 16.0);
        assert!((h - want).abs() / want < 0.12, "hole bottom flux {h}, square-aperture view factor {want}");
    }

    #[test]
    fn re_emission_lifts_the_bottom_above_the_line_of_sight_value() {
        let g = VoxelGrid::trench(48, 4, 48, nm(10.0), 16, 16, 4).unwrap();
        let direct = solve(&g, &cfg(300_000, 1.0)).unwrap().mean_top_flux(23..25, 0..4).unwrap();
        let bouncing = solve(&g, &cfg(300_000, 0.05)).unwrap().mean_top_flux(23..25, 0..4).unwrap();
        assert!(bouncing > 1.15 * direct, "particles that bounce reach the bottom again: {bouncing} vs {direct}");
    }

    #[test]
    fn a_tilted_beam_leaves_the_bottom_of_a_deep_trench_in_shadow() {
        // 45 degrees toward +x crosses 32 cells of depth in 32 cells of x: wider than the 16-cell slot
        let g = VoxelGrid::trench(48, 4, 48, nm(10.0), 16, 32, 4).unwrap();
        let c = Config {
            source: Source3d { angular: Angular::Beam, tilt_deg: 45.0, azimuth_deg: 0.0, spread_deg: 0.0 },
            ..cfg(200_000, 1.0)
        };
        let s = solve(&g, &c).unwrap();
        let bottom = s.mean_top_flux(18..30, 0..4).unwrap();
        assert!(bottom < 1e-9, "a collimated beam at 45 degrees cannot reach a 1:2 trench bottom: {bottom}");
        let field = s.mean_top_flux(0..8, 0..4).unwrap();
        assert!((field - 1.0).abs() < 0.03, "the open field is unaffected by the tilt: {field}");
    }

    #[test]
    fn a_solve_is_bit_identical_across_runs_and_chunk_schedules() {
        let g = VoxelGrid::hole(20, 20, 32, nm(10.0), 8, 12, 4).unwrap();
        let one = solve(&g, &cfg(40_000, 0.4)).unwrap();
        let again = solve(&g, &cfg(40_000, 0.4)).unwrap();
        // the chunks run in reverse, as if scheduled on workers that finish in any order
        let plan = Plan::new(&g, &cfg(40_000, 0.4)).unwrap();
        let parts: Vec<Partial> = (0..plan.chunk_count()).rev().map(|c| plan.run_chunk(c)).collect();
        let scrambled = plan.finish(parts).unwrap();
        let bits = |s: &Solution| s.flux_values().iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&one), bits(&again), "the same inputs give the same output");
        assert_eq!(bits(&one), bits(&scrambled), "the schedule changes the time, never the answer");
        let other = solve(&g, &Config { seed: 2, ..cfg(40_000, 0.4) }).unwrap();
        assert_ne!(bits(&one), bits(&other), "a different seed is a different sample");
    }

    #[test]
    fn finishing_with_a_missing_or_repeated_chunk_is_refused() {
        let g = VoxelGrid::from_heights(8, 8, 8, nm(10.0), &vec![nm(20.0); 64]).unwrap();
        let plan = Plan::new(&g, &cfg(1_000, 1.0)).unwrap();
        let mut parts: Vec<Partial> = (0..plan.chunk_count()).map(|c| plan.run_chunk(c)).collect();
        parts.pop();
        assert!(plan.clone().finish(parts).is_err(), "a chunk is missing");
        let mut twice: Vec<Partial> = (0..plan.chunk_count()).map(|c| plan.run_chunk(c)).collect();
        twice[3] = plan.run_chunk(2);
        assert!(plan.finish(twice).is_err(), "chunk 2 appears twice and chunk 3 not at all");
    }

    #[test]
    fn a_grid_or_run_outside_its_domain_is_refused() {
        let flat = vec![true; 16 * 16 * 4];
        assert!(VoxelGrid::new(16, 16, 4, nm(10.0), flat).is_err(), "no empty top layer");
        let mut open = vec![false; 16 * 16 * 4];
        assert!(VoxelGrid::new(16, 16, 4, nm(10.0), open.clone()).is_err(), "no solid floor");
        for v in open.iter_mut().take(16 * 16) {
            *v = true;
        }
        assert!(VoxelGrid::new(16, 16, 4, nm(10.0), open.clone()).is_ok());
        assert!(VoxelGrid::new(16, 16, 4, nm(f64::NAN), open.clone()).is_err(), "cell is not a number");
        assert!(VoxelGrid::new(16, 16, 4, nm(-1.0), open.clone()).is_err(), "negative cell");
        assert!(VoxelGrid::new(0, 16, 4, nm(10.0), Vec::new()).is_err(), "empty axis");
        assert!(VoxelGrid::new(1024, 1024, 4, nm(10.0), Vec::new()).is_err(), "over the voxel cap");
        assert!(VoxelGrid::from_heights(4, 4, 8, nm(10.0), &vec![nm(1.0e6); 16]).is_err(), "taller than the grid");

        let g = VoxelGrid::from_heights(8, 8, 8, nm(10.0), &vec![nm(20.0); 64]).unwrap();
        assert!(solve(&g, &cfg(0, 1.0)).is_err(), "no rays");
        assert!(solve(&g, &cfg(MAX_RAYS + 1, 1.0)).is_err(), "too many rays");
        assert!(solve(&g, &cfg(100, 0.0)).is_err(), "nothing ever sticks");
        assert!(solve(&g, &cfg(100, 1.5)).is_err(), "sticking above one");
        assert!(solve(&g, &cfg(100, f64::NAN)).is_err(), "sticking is not a number");
        let bad_beam = Source3d { angular: Angular::Beam, tilt_deg: 90.0, ..Source3d::default() };
        assert!(solve(&g, &Config { source: bad_beam, ..cfg(100, 1.0) }).is_err(), "a grazing beam");
    }

    #[test]
    fn a_cut_that_does_not_fit_is_refused() {
        assert!(VoxelGrid::trench(16, 4, 20, nm(10.0), 16, 8, 4).is_err(), "trench as wide as the grid");
        assert!(VoxelGrid::trench(16, 4, 12, nm(10.0), 4, 8, 4).is_err(), "no headroom above the surface");
        assert!(VoxelGrid::hole(16, 8, 20, nm(10.0), 4, 0, 4).is_err(), "zero depth");
    }
}
