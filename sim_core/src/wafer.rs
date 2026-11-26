//! Sparse Voxel Octree (SVO) wafer memory state.
//!
//! A dense 3D grid of a 300mm wafer at nanometer voxel resolution would
//! require terabytes of memory. Instead, the wafer is represented as an
//! octree: any node whose entire volume is a single homogeneous material
//! (bulk silicon, or air above the surface) collapses into one leaf, and
//! only regions with actual nanoscale structure (trenches, deposited
//! layers, doped pockets) get subdivided down to fine-grained leaves.

use serde::{Deserialize, Serialize};

/// A material tag for a voxel/leaf. Kept as a small integer-backed enum so
/// leaves stay cheap; a side table (see `materials::MaterialDb`) maps these
/// to physical constants.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Material {
    Air,
    Silicon,
    SiliconDioxide,
    SiliconCarbide,
    GalliumNitride,
    Photoresist,
    Metal,
    Dopant(u8), // coarse dopant concentration bucket, for voxels with implant data
}

/// Axis-aligned bounding box, in nanometers, for one octree node.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Bounds {
    pub min: (f64, f64, f64),
    pub max: (f64, f64, f64),
}

impl Bounds {
    pub fn center(&self) -> (f64, f64, f64) {
        (
            (self.min.0 + self.max.0) / 2.0,
            (self.min.1 + self.max.1) / 2.0,
            (self.min.2 + self.max.2) / 2.0,
        )
    }

    pub fn octant(&self, index: usize) -> Bounds {
        let c = self.center();
        let (x_lo, x_hi) = if index & 1 == 0 { (self.min.0, c.0) } else { (c.0, self.max.0) };
        let (y_lo, y_hi) = if index & 2 == 0 { (self.min.1, c.1) } else { (c.1, self.max.1) };
        let (z_lo, z_hi) = if index & 4 == 0 { (self.min.2, c.2) } else { (c.2, self.max.2) };
        Bounds { min: (x_lo, y_lo, z_lo), max: (x_hi, y_hi, z_hi) }
    }
}

/// A node in the sparse voxel octree: either a homogeneous leaf, or an
/// internal node with (up to) 8 children.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum SvoNode {
    Leaf(Material),
    Internal(Box<[SvoNode; 8]>),
}

/// Circular physical wafer domain, independent of the cuboidal SVO storage root.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
pub struct WaferDomain {
    pub center_x_nm: f64,
    pub center_y_nm: f64,
    pub radius_nm: f64,
}

impl WaferDomain {
    pub fn contains_xy(&self, x_nm: f64, y_nm: f64) -> bool {
        if !x_nm.is_finite() || !y_nm.is_finite() || !self.radius_nm.is_finite() || self.radius_nm <= 0.0 { return false; }
        let dx = x_nm - self.center_x_nm;
        let dy = y_nm - self.center_y_nm;
        dx * dx + dy * dy <= self.radius_nm * self.radius_nm
    }

    pub fn intersects_xy_bounds(&self, b: Bounds) -> bool {
        let x = self.center_x_nm.clamp(b.min.0, b.max.0);
        let y = self.center_y_nm.clamp(b.min.1, b.max.1);
        let dx = x - self.center_x_nm;
        let dy = y - self.center_y_nm;
        dx * dx + dy * dy <= self.radius_nm * self.radius_nm
    }
}

/// The wafer's full SVO state plus its overall bounds and minimum leaf size
/// (the finest resolution the tree will subdivide to).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WaferState {
    pub root: SvoNode,
    pub bounds: Bounds,
    pub domain: WaferDomain,
    /// Smallest edge length (nm) a node may have before subdivision stops.
    pub min_leaf_size_nm: f64,
    /// Z of the original polished substrate surface (nm). Equals `bounds.max.2`
    /// for legacy wafers without headroom; lies below it when the root box
    /// carries an air headroom that films can grow into.
    #[serde(default)]
    pub surface_z_nm: f64,
    /// Physical substrate thickness (nm). The substrate spans z = 0 .. this value.
    #[serde(default)]
    pub substrate_thickness_nm: f64,
    /// When true, `bounds` is a cubic high-resolution *tile* around the
    /// process site and everything outside it is an implicit homogeneous bulk
    /// (silicon below the surface, air above, clipped to `domain`). A 300 mm
    /// × 0.5 mm wafer cannot be resolved to nanometre cells by a single
    /// octree without exploding the node count on every vertical sidewall, so
    /// only the tile is subdivided.
    #[serde(default)]
    pub implicit_bulk: bool,
    /// Process-state side data that is not geometry (surface state, implant
    /// profiles, film stress). Travels with the wafer through undo/redo.
    #[serde(default)]
    pub meta: WaferMeta,
}

/// Surface condition tracked by preparation / cleaning stages.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SurfaceState {
    pub native_oxide_nm: f64,
    pub particle_density_cm2: f64,
    pub residue_fraction: f64,
    pub contact_angle_deg: Option<f64>,
}

impl Default for SurfaceState {
    fn default() -> Self {
        Self { native_oxide_nm: 1.5, particle_density_cm2: 5.0e3, residue_fraction: 0.0, contact_angle_deg: None }
    }
}

/// Analytic description of one implant (+ anneal) so metrology can reference
/// the dopant profile without a dense concentration grid.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DopantProfile {
    pub species: String,
    pub dose_cm2: f64,
    pub energy_kev: f64,
    pub rp_nm: f64,
    /// Straggle after any anneal broadening, nm.
    pub sigma_nm: f64,
    pub anneal_temp_c: f64,
    pub anneal_min: f64,
    pub peak_cm3: f64,
    pub junction_depth_nm: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FilmStress {
    pub material: Material,
    pub thickness_nm: f64,
    pub stress_mpa: f64,
}

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
pub struct WaferMeta {
    pub surface: SurfaceState,
    pub dopant_profiles: Vec<DopantProfile>,
    pub film_stress: Vec<FilmStress>,
}

/// One vertical material run in a column, nm.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Run {
    pub z0: f64,
    pub z1: f64,
    pub material: Material,
}

#[derive(Debug, thiserror::Error, PartialEq)]
pub enum WaferError {
    #[error("wafer diameter must be finite and > 0 mm")]
    InvalidDiameter,
    #[error("wafer thickness must be finite and > 0 um")]
    InvalidThickness,
    #[error("minimum leaf size must be finite and > 0 nm")]
    InvalidLeafSize,
    #[error("tile must be centred inside the wafer and no thicker than the substrate allows")]
    InvalidTile,
}

impl WaferState {
    /// Create a fresh wafer: a single homogeneous silicon substrate leaf.
    /// Invalid dimensions are rejected immediately rather than allowing
    /// recursive SVO subdivision to run with NaN/zero resolution.
    pub fn new_substrate(diameter_mm: f64, thickness_um: f64, min_leaf_size_nm: f64) -> Self {
        Self::try_new_substrate(diameter_mm, thickness_um, min_leaf_size_nm)
            .expect("invalid wafer dimensions or SVO leaf size")
    }

    pub fn try_new_substrate(
        diameter_mm: f64,
        thickness_um: f64,
        min_leaf_size_nm: f64,
    ) -> Result<Self, WaferError> {
        if !diameter_mm.is_finite() || diameter_mm <= 0.0 {
            return Err(WaferError::InvalidDiameter);
        }
        if !thickness_um.is_finite() || thickness_um <= 0.0 {
            return Err(WaferError::InvalidThickness);
        }
        if !min_leaf_size_nm.is_finite() || min_leaf_size_nm <= 0.0 {
            return Err(WaferError::InvalidLeafSize);
        }
        let r_nm = diameter_mm * 1e6 / 2.0;
        let t_nm = thickness_um * 1e3;
        Ok(Self {
            root: SvoNode::Leaf(Material::Silicon),
            bounds: Bounds {
                min: (-r_nm, -r_nm, 0.0),
                max: (r_nm, r_nm, t_nm),
            },
            domain: WaferDomain { center_x_nm: 0.0, center_y_nm: 0.0, radius_nm: r_nm },
            min_leaf_size_nm,
            surface_z_nm: t_nm,
            substrate_thickness_nm: t_nm,
            implicit_bulk: false,
            meta: WaferMeta::default(),
        })
    }

    /// Create a wafer with a cubic high-resolution tile at `tile_center`.
    ///
    /// The tile edge is `min_leaf_size_nm * 2^levels` so the finest octree
    /// cell is exactly `min_leaf_size_nm` in X, Y *and* Z. The substrate
    /// surface sits on the tile's mid-height split plane, so the blank tile is
    /// nine nodes (four silicon octants, four air octants and the root).
    pub fn new_wafer_tile(diameter_mm: f64, thickness_um: f64, min_leaf_size_nm: f64, levels: u32, tile_center: (f64, f64)) -> Result<Self, WaferError> {
        let mut w = Self::try_new_substrate(diameter_mm, thickness_um, min_leaf_size_nm)?;
        if !(tile_center.0.is_finite() && tile_center.1.is_finite()) { return Err(WaferError::InvalidTile); }
        let levels = levels.clamp(1, 14);
        let edge = min_leaf_size_nm * (1u64 << levels) as f64;
        let t = w.substrate_thickness_nm;
        let h = edge * 0.5;
        if edge > t * 2.0 || !w.domain.contains_xy(tile_center.0, tile_center.1) { return Err(WaferError::InvalidTile); }
        w.bounds = Bounds { min: (tile_center.0 - h, tile_center.1 - h, t - h), max: (tile_center.0 + h, tile_center.1 + h, t + h) };
        w.implicit_bulk = true;
        let children: [SvoNode; 8] = std::array::from_fn(|i| {
            if i & 4 == 0 { SvoNode::Leaf(Material::Silicon) } else { SvoNode::Leaf(Material::Air) }
        });
        w.root = SvoNode::Internal(Box::new(children));
        Ok(w)
    }

    /// Centre of the high-resolution tile (wafer-local nm) at the surface.
    pub fn tile_origin(&self) -> (f64, f64, f64) {
        let c = self.bounds.center();
        (c.0, c.1, self.surface_z_nm)
    }

    pub fn tile_contains_xy(&self, x: f64, y: f64) -> bool {
        x >= self.bounds.min.0 && x <= self.bounds.max.0 && y >= self.bounds.min.1 && y <= self.bounds.max.1
    }

    /// Edge length (nm) of the finest octree cell in X/Y: the root is halved
    /// until its X edge is no larger than `min_leaf_size_nm`, exactly as the
    /// edit routines do. Process grids use this pitch so every grid cell maps
    /// onto exactly one finest-level column.
    pub fn leaf_pitch_nm(&self) -> f64 {
        let mut edge = self.bounds.max.0 - self.bounds.min.0;
        let mut guard = 0;
        while edge > self.min_leaf_size_nm && guard < 64 { edge *= 0.5; guard += 1; }
        edge
    }

    /// Material at a point, or `None` outside the root box or the circular domain.
    pub fn material_at(&self, x: f64, y: f64, z: f64) -> Option<Material> {
        if !(x.is_finite() && y.is_finite() && z.is_finite()) { return None; }
        let b = self.bounds;
        if !self.domain.contains_xy(x, y) { return None; }
        let inside_tile = x >= b.min.0 && x <= b.max.0 && y >= b.min.1 && y <= b.max.1 && z >= b.min.2 && z <= b.max.2;
        if !inside_tile {
            if self.implicit_bulk && z >= 0.0 {
                return Some(if z < self.surface_z_nm { Material::Silicon } else { Material::Air });
            }
            return None;
        }
        let mut node = &self.root;
        let mut nb = b;
        loop {
            match node {
                SvoNode::Leaf(m) => return Some(*m),
                SvoNode::Internal(children) => {
                    let c = nb.center();
                    let idx = (x >= c.0) as usize | (((y >= c.1) as usize) << 1) | (((z >= c.2) as usize) << 2);
                    nb = nb.octant(idx);
                    node = &children[idx];
                }
            }
        }
    }

    /// Vertical material stack at (x, y): contiguous runs from the root bottom
    /// to the root top, adjacent equal materials merged. Empty when (x, y) is
    /// outside the root box.
    pub fn column_runs(&self, x: f64, y: f64) -> Vec<Run> {
        let b = self.bounds;
        let mut out: Vec<Run> = Vec::new();
        if !(x.is_finite() && y.is_finite()) || x < b.min.0 || x > b.max.0 || y < b.min.1 || y > b.max.1 { return out; }
        fn walk(node: &SvoNode, nb: Bounds, x: f64, y: f64, out: &mut Vec<Run>) {
            match node {
                SvoNode::Leaf(m) => {
                    if let Some(last) = out.last_mut() {
                        if last.material == *m { last.z1 = nb.max.2; return; }
                    }
                    out.push(Run { z0: nb.min.2, z1: nb.max.2, material: *m });
                }
                SvoNode::Internal(children) => {
                    let c = nb.center();
                    let ix = (x >= c.0) as usize | (((y >= c.1) as usize) << 1);
                    for zb in 0..2 {
                        let idx = ix | (zb << 2);
                        walk(&children[idx], nb.octant(idx), x, y, out);
                    }
                }
            }
        }
        walk(&self.root, b, x, y, &mut out);
        out
    }

    /// Deterministic 64-bit structural fingerprint (FNV-1a) of the geometry,
    /// domain and process metadata that metrology depends on. Two wafers with
    /// the same fingerprint are the same committed state.
    pub fn fingerprint(&self) -> u64 {
        const PRIME: u64 = 0x0000_0100_0000_01B3;
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |bytes: &[u8]| { for b in bytes { h ^= *b as u64; h = h.wrapping_mul(PRIME); } };
        for v in [self.bounds.min.0, self.bounds.min.1, self.bounds.min.2, self.bounds.max.0, self.bounds.max.1, self.bounds.max.2,
                  self.domain.radius_nm, self.min_leaf_size_nm, self.surface_z_nm, self.substrate_thickness_nm] { feed(&v.to_bits().to_le_bytes()); }
        fn code(m: Material) -> [u8; 2] {
            match m {
                Material::Air => [0, 0], Material::Silicon => [1, 0], Material::SiliconDioxide => [2, 0],
                Material::SiliconCarbide => [3, 0], Material::GalliumNitride => [4, 0], Material::Photoresist => [5, 0],
                Material::Metal => [6, 0], Material::Dopant(b) => [7, b],
            }
        }
        fn walk(n: &SvoNode, feed: &mut impl FnMut(&[u8])) {
            match n {
                SvoNode::Leaf(m) => { feed(&[b'L']); feed(&code(*m)); }
                SvoNode::Internal(c) => { feed(&[b'I']); for ch in c.iter() { walk(ch, feed); } }
            }
        }
        walk(&self.root, &mut feed);
        for p in &self.meta.dopant_profiles {
            for v in [p.dose_cm2, p.rp_nm, p.sigma_nm, p.junction_depth_nm] { feed(&v.to_bits().to_le_bytes()); }
        }
        feed(&self.meta.surface.particle_density_cm2.to_bits().to_le_bytes());
        feed(&self.meta.surface.residue_fraction.to_bits().to_le_bytes());
        // surface state and film stress change the state, so they change the id
        feed(&self.meta.surface.native_oxide_nm.to_bits().to_le_bytes());
        for f in &self.meta.film_stress {
            // Material carries Dopant(u8) and has no repr, so `as u8` is not a
            // legal cast; a stable code is matched instead.
            let code: u8 = match f.material {
                Material::Air => 0,
                Material::Silicon => 1,
                Material::SiliconDioxide => 2,
                Material::SiliconCarbide => 3,
                Material::GalliumNitride => 4,
                Material::Photoresist => 5,
                Material::Metal => 6,
                Material::Dopant(_) => 7,
            };
            feed(&[code]);
            feed(&f.thickness_nm.to_bits().to_le_bytes());
            feed(&f.stress_mpa.to_bits().to_le_bytes());
        }
        h
    }

    /// Set the material within a region, subdividing as needed (and
    /// re-collapsing homogeneous subtrees where possible) — e.g. after a
    /// deposition, etch, or implant step mutates part of the wafer.
    pub fn set_region(&mut self, region: Bounds, material: Material) {
        Self::set_region_rec(&mut self.root, self.bounds, region, material, self.min_leaf_size_nm, self.domain);
    }

    /// Spatially clipped deposition: modifies only eligible solid/void leaves
    /// intersecting the operation volume. This is a localized geometry operator,
    /// not a calibrated film-growth or surface-normal model.
    pub fn deposit_region(&mut self, region: Bounds, material: Material) -> usize {
        self.set_region_filtered(region, material, 0)
    }

    /// Apply a photoresist mask volume to process-compatible exposed materials.
    /// Optical exposure, development, and mask transfer remain separate stages.
    pub fn photoresist_region(&mut self, region: Bounds) -> usize {
        self.set_region_filtered(region, Material::Photoresist, 1)
    }

    /// Remove eligible material in the selected etch volume, preserving voids.
    pub fn etch_region(&mut self, region: Bounds) -> usize {
        self.set_region_filtered(region, Material::Air, 2)
    }

    fn set_region_filtered(&mut self, region: Bounds, target: Material, operation: u8) -> usize {
        let mut eligible_leaves = 0usize;
        self.visit_leaves_in_region(region, |_, material| {
            if Self::material_eligible(material, operation) { eligible_leaves += 1; }
        });
        Self::set_region_filtered_rec(&mut self.root, self.bounds, region, target, self.min_leaf_size_nm, self.domain, operation);
        eligible_leaves
    }

    fn material_eligible(material: Material, operation: u8) -> bool {
        match operation {
            0 => !matches!(material, Material::Photoresist | Material::Dopant(_)),
            1 => matches!(material, Material::Silicon | Material::SiliconDioxide | Material::Metal | Material::GalliumNitride | Material::SiliconCarbide),
            2 => !matches!(material, Material::Air),
            _ => false,
        }
    }

    fn set_region_filtered_rec(node: &mut SvoNode, node_bounds: Bounds, region: Bounds, target: Material, min_leaf: f64, domain: WaferDomain, operation: u8) {
        if !overlaps(&node_bounds, &region) || !domain.intersects_xy_bounds(node_bounds) { return; }
        if let SvoNode::Leaf(existing) = node {
            if !Self::material_eligible(*existing, operation) { return; }
            if contains(&region, &node_bounds) { *node = SvoNode::Leaf(target); return; }
            let edge = node_bounds.max.0 - node_bounds.min.0;
            if edge <= min_leaf {
                let center = node_bounds.center();
                if point_in(&region, center) && domain.contains_xy(center.0, center.1) { *node = SvoNode::Leaf(target); }
                return;
            }
            let existing = *existing;
            let children: [SvoNode; 8] = std::array::from_fn(|_| SvoNode::Leaf(existing));
            *node = SvoNode::Internal(Box::new(children));
        }
        if let SvoNode::Internal(children) = node {
            for i in 0..8 { Self::set_region_filtered_rec(&mut children[i], node_bounds.octant(i), region, target, min_leaf, domain, operation); }
        }
        collapse_if_homogeneous(node);
    }

    fn set_region_rec(node: &mut SvoNode, node_bounds: Bounds, region: Bounds, material: Material, min_leaf: f64, domain: WaferDomain) {
        if !overlaps(&node_bounds, &region) || !domain.intersects_xy_bounds(node_bounds) {
            return;
        }

        // Store edits sparsely. The physical-domain mask clips rendering and
        // physical queries at the circular edge; forcing octree subdivision
        // all the way to nanometre leaves along the 300 mm circumference would
        // be computationally impossible. Regions completely outside the disk
        // are rejected above; boundary-crossing leaves remain clipped by domain.
        if contains(&region, &node_bounds) {
            *node = SvoNode::Leaf(material);
            return;
        }

        let edge = node_bounds.max.0 - node_bounds.min.0;
        if edge <= min_leaf {
            // Can't subdivide further; approximate by taking the material
            // if the region covers the node's center.
            let c = node_bounds.center();
            if point_in(&region, c) && domain.contains_xy(c.0, c.1) {
                *node = SvoNode::Leaf(material);
            }
            return;
        }

        // Ensure this node is subdivided, preserving its current material
        // as the starting point for all 8 children.
        if let SvoNode::Leaf(existing) = node {
            let existing = *existing;
            let children: [SvoNode; 8] = std::array::from_fn(|_| SvoNode::Leaf(existing));
            *node = SvoNode::Internal(Box::new(children));
        }

        if let SvoNode::Internal(children) = node {
            for i in 0..8 {
                let child_bounds = node_bounds.octant(i);
                Self::set_region_rec(&mut children[i], child_bounds, region, material, min_leaf, domain);
            }
        }

        collapse_if_homogeneous(node);
    }

    /// Visit homogeneous leaves without materializing a dense voxel array.
    /// Physics operators should use this sparse traversal (or a localized
    /// region query) as their geometry input rather than allocating Wafer2d/3d.
    pub fn visit_leaves(&self, mut visitor: impl FnMut(Bounds, Material)) {
        fn walk(node: &SvoNode, bounds: Bounds, visitor: &mut impl FnMut(Bounds, Material)) {
            match node {
                SvoNode::Leaf(material) => visitor(bounds, *material),
                SvoNode::Internal(children) => {
                    for (index, child) in children.iter().enumerate() {
                        walk(child, bounds.octant(index), visitor);
                    }
                }
            }
        }
        walk(&self.root, self.bounds, &mut visitor);
    }

    /// Visit only leaves intersecting `region`; useful for sparse/narrow-band
    /// process updates that should not scan unrelated wafer regions.
    pub fn visit_leaves_in_region(
        &self,
        region: Bounds,
        mut visitor: impl FnMut(Bounds, Material),
    ) {
        fn walk(
            node: &SvoNode,
            bounds: Bounds,
            region: Bounds,
            visitor: &mut impl FnMut(Bounds, Material),
        ) {
            if !overlaps(&bounds, &region) {
                return;
            }
            match node {
                SvoNode::Leaf(material) => visitor(bounds, *material),
                SvoNode::Internal(children) => {
                    for (index, child) in children.iter().enumerate() {
                        walk(child, bounds.octant(index), region, visitor);
                    }
                }
            }
        }
        walk(&self.root, self.bounds, region, &mut visitor);
    }

    /// Estimate physical material volume by integrating the circular XY mask
    /// over each leaf. This deterministic midpoint quadrature avoids counting
    /// square-root corner volume outside the wafer. Boundary accuracy increases
    /// with `samples_per_axis`; this remains an estimate, not exact CAD geometry.
    pub fn material_volume_nm3(&self, target: Material) -> f64 {
        self.material_volume_nm3_with_samples(target, 16)
    }

    pub fn material_volume_nm3_with_samples(&self, target: Material, samples_per_axis: usize) -> f64 {
        let n = samples_per_axis.max(1);
        let mut volume = 0.0;
        self.visit_leaves(|b, material| {
            if material != target { return; }
            let dx = (b.max.0 - b.min.0).max(0.0);
            let dy = (b.max.1 - b.min.1).max(0.0);
            let dz = (b.max.2 - b.min.2).max(0.0);
            if dx == 0.0 || dy == 0.0 || dz == 0.0 { return; }
            let mut inside = 0usize;
            for ix in 0..n {
                let x = b.min.0 + (ix as f64 + 0.5) * dx / n as f64;
                for iy in 0..n {
                    let y = b.min.1 + (iy as f64 + 0.5) * dy / n as f64;
                    if self.domain.contains_xy(x, y) { inside += 1; }
                }
            }
            volume += dx * dy * dz * inside as f64 / (n * n) as f64;
        });
        volume
    }

    /// CPU picking guard shared by all viewport modes. A hit is valid only if
    /// its world-space XY coordinate lies inside the physical wafer disk.
    pub fn accepts_pick_xy(&self, x_nm: f64, y_nm: f64) -> bool {
        self.domain.contains_xy(x_nm, y_nm)
    }

    /// Count leaves carrying a particular material tag.
    pub fn material_leaf_count(&self, target: Material) -> usize {
        let mut count = 0;
        self.visit_leaves(|_, material| if material == target { count += 1; });
        count
    }

    /// Count total nodes currently materialized (leaves + internal),
    /// primarily for diagnostics / memory-wall reporting in the UI status
    /// bar.
    pub fn node_count(&self) -> usize {
        fn count(n: &SvoNode) -> usize {
            match n {
                SvoNode::Leaf(_) => 1,
                SvoNode::Internal(children) => 1 + children.iter().map(count).sum::<usize>(),
            }
        }
        count(&self.root)
    }
}

fn overlaps(a: &Bounds, b: &Bounds) -> bool {
    a.min.0 < b.max.0 && a.max.0 > b.min.0 &&
    a.min.1 < b.max.1 && a.max.1 > b.min.1 &&
    a.min.2 < b.max.2 && a.max.2 > b.min.2
}

fn contains(outer: &Bounds, inner: &Bounds) -> bool {
    outer.min.0 <= inner.min.0 && outer.max.0 >= inner.max.0 &&
    outer.min.1 <= inner.min.1 && outer.max.1 >= inner.max.1 &&
    outer.min.2 <= inner.min.2 && outer.max.2 >= inner.max.2
}

fn point_in(b: &Bounds, p: (f64, f64, f64)) -> bool {
    p.0 >= b.min.0 && p.0 <= b.max.0 &&
    p.1 >= b.min.1 && p.1 <= b.max.1 &&
    p.2 >= b.min.2 && p.2 <= b.max.2
}

/// If every child of an internal node is a leaf of the same material,
/// collapse the node back into a single leaf — this is what keeps bulk
/// silicon/air regions cheap after many localized edits.
fn collapse_if_homogeneous(node: &mut SvoNode) {
    if let SvoNode::Internal(children) = node {
        if let SvoNode::Leaf(first) = &children[0] {
            let first = *first;
            let all_same = children.iter().all(|c| matches!(c, SvoNode::Leaf(m) if *m == first));
            if all_same {
                *node = SvoNode::Leaf(first);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_fingerprint_separates_states_that_differ() {
        // The fingerprint decides whether a committed stage is stale, and the earlier build had no
        // test for it: reverting F-17's film-material feed left every wafer test
        // passing. The material case is the one that isolates that feed — two films
        // of the same thickness and stress differ only by their material.
        let base = WaferState::new_substrate(1.0, 1.0, 100_000.0);
        let id = base.fingerprint();
        assert_eq!(base.fingerprint(), id, "an unchanged state keeps its id");
        let mut oxide = base.clone();
        oxide.meta.film_stress.push(FilmStress { material: Material::SiliconDioxide, thickness_nm: 100.0, stress_mpa: -250.0 });
        assert_ne!(oxide.fingerprint(), id, "a film must change the fingerprint");
        let mut metal = base.clone();
        metal.meta.film_stress.push(FilmStress { material: Material::Metal, thickness_nm: 100.0, stress_mpa: -250.0 });
        assert_ne!(metal.fingerprint(), oxide.fingerprint(), "the film material must change the fingerprint");
        let mut rough = base.clone();
        rough.meta.surface.native_oxide_nm += 1.0;
        assert_ne!(rough.fingerprint(), id, "the native oxide must change the fingerprint");
    }

    #[test]
    fn fresh_wafer_is_a_single_leaf() {
        let wafer = WaferState::new_substrate(300.0, 500.0, 1.0);
        assert_eq!(wafer.node_count(), 1);
    }

    #[test]
    fn invalid_dimensions_are_rejected_by_fallible_constructor() {
        assert!(matches!(WaferState::try_new_substrate(0.0, 500.0, 10.0), Err(WaferError::InvalidDiameter)));
        assert!(matches!(WaferState::try_new_substrate(300.0, f64::NAN, 10.0), Err(WaferError::InvalidThickness)));
        assert!(matches!(WaferState::try_new_substrate(300.0, 500.0, 0.0), Err(WaferError::InvalidLeafSize)));
    }

    #[test]
    fn material_volume_is_conserved_by_tag_reassignment() {
        let mut wafer = WaferState::new_substrate(1.0, 1.0, 100_000.0);
        let initial = wafer.material_volume_nm3(Material::Silicon);
        wafer.set_region(Bounds { min: (-100_000.0, -100_000.0, 0.0), max: (100_000.0, 100_000.0, 100.0) }, Material::SiliconDioxide);
        let silicon = wafer.material_volume_nm3(Material::Silicon);
        let oxide = wafer.material_volume_nm3(Material::SiliconDioxide);
        assert!((silicon + oxide - initial).abs() < initial * 0.05, "circular-boundary quadrature should conserve volume within 5%");
        assert!(oxide > 0.0);
    }

    #[test]
    fn localized_edit_subdivides_then_stays_small() {
        let mut wafer = WaferState::new_substrate(300.0, 500.0, 10.0);
        let region = Bounds { min: (0.0, 0.0, 0.0), max: (100.0, 100.0, 100.0) };
        wafer.set_region(region, Material::SiliconDioxide);
        // Should have subdivided (more than 1 node) but stayed many orders
        // of magnitude below a dense grid across the whole 300mm wafer at
        // 10nm voxels (which would be on the order of 10^19 voxels).
        assert!(wafer.node_count() > 1);
        assert!(wafer.node_count() < 1_000_000);
    }

    #[test]
    fn sparse_leaf_visitation_covers_the_root_once() {
        let wafer = WaferState::new_substrate(300.0, 500.0, 10.0);
        let mut seen = Vec::new();
        wafer.visit_leaves(|bounds, material| seen.push((bounds, material)));
        assert_eq!(seen.len(), 1);
        assert_eq!(seen[0].1, Material::Silicon);
    }

    #[test]
    fn localized_leaf_query_does_not_visit_disjoint_region() {
        let wafer = WaferState::new_substrate(300.0, 500.0, 10.0);
        let mut count = 0;
        wafer.visit_leaves_in_region(
            Bounds { min: (1.0e9, 1.0e9, 0.0), max: (1.1e9, 1.1e9, 100.0) },
            |_, _| count += 1,
        );
        assert_eq!(count, 0);
    }

    #[test]
    fn full_wafer_overwrite_recollapses_to_one_leaf() {
        let mut wafer = WaferState::new_substrate(300.0, 500.0, 10.0);
        wafer.set_region(wafer.bounds, Material::Air);
        assert_eq!(wafer.node_count(), 1);
    }
    #[test]
    fn domain_is_circular_and_uses_actual_diameter() {
        let wafer = WaferState::new_substrate(300.0, 500.0, 1_000_000.0);
        assert!((wafer.domain.radius_nm - 150_000_000.0).abs() < 1e-6);
        assert!(wafer.domain.contains_xy(150_000_000.0, 0.0));
        assert!(!wafer.domain.contains_xy(149_000_000.0, 149_000_000.0));
    }

    #[test]
    fn edits_outside_circle_do_not_subdivide_or_change_tree() {
        let mut wafer = WaferState::new_substrate(300.0, 500.0, 1_000_000.0);
        let before = wafer.node_count();
        wafer.set_region(
            Bounds { min: (140_000_000.0, 140_000_000.0, 0.0), max: (150_000_000.0, 150_000_000.0, 100_000.0) },
            Material::Metal,
        );
        assert_eq!(wafer.node_count(), before);
    }

    #[test]
    fn edge_region_edit_is_clipped_by_domain_center_sampling() {
        let mut wafer = WaferState::new_substrate(10.0, 1.0, 500_000.0);
        wafer.set_region(
            Bounds { min: (4_000_000.0, 4_000_000.0, 0.0), max: (5_000_000.0, 5_000_000.0, 500.0) },
            Material::Metal,
        );
        assert_eq!(wafer.material_leaf_count(Material::Metal), 0);
    }

    #[test]
    fn circular_volume_is_close_to_pi_over_four_of_bounding_cylinder() {
        let wafer = WaferState::new_substrate(2.0, 1.0, 10_000_000.0);
        let estimated = wafer.material_volume_nm3_with_samples(Material::Silicon, 128);
        let square_volume = 2_000_000.0_f64.powi(2) * 1_000.0;
        let ratio = estimated / square_volume;
        assert!((ratio - std::f64::consts::PI / 4.0).abs() < 0.015, "ratio={ratio}");
    }

    #[test]
    fn process_operators_are_spatial_and_material_aware() {
        let mut wafer = WaferState::new_substrate(1.0, 1.0, 100_000.0);
        let left = Bounds { min: (-300_000.0, -100_000.0, 500.0), max: (-100_000.0, 100_000.0, 1_000.0) };
        let right = Bounds { min: (100_000.0, -100_000.0, 500.0), max: (300_000.0, 100_000.0, 1_000.0) };
        assert!(wafer.deposit_region(left, Material::SiliconDioxide) > 0);
        assert!(wafer.material_leaf_count(Material::SiliconDioxide) > 0);
        // PR targets compatible solids only and does not overwrite the void tag.
        assert!(wafer.photoresist_region(left) > 0);
        assert!(wafer.material_leaf_count(Material::Photoresist) > 0);
        // Etching is local: it removes material in its region while leaving a
        // spatially disjoint oxide/PR feature untouched.
        wafer.set_region(right, Material::SiliconDioxide);
        assert!(wafer.etch_region(right) > 0);
        assert!(wafer.material_leaf_count(Material::Air) > 0);
        assert!(wafer.material_leaf_count(Material::Photoresist) > 0);
    }

    #[test]
    fn deposition_does_not_overwrite_photoresist() {
        let mut wafer = WaferState::new_substrate(1.0, 1.0, 100_000.0);
        let region = Bounds { min: (-100_000.0, -100_000.0, 500.0), max: (100_000.0, 100_000.0, 1_000.0) };
        wafer.set_region(region, Material::Photoresist);
        let pr_before = wafer.material_leaf_count(Material::Photoresist);
        wafer.deposit_region(region, Material::SiliconDioxide);
        assert_eq!(wafer.material_leaf_count(Material::Photoresist), pr_before);
    }

    #[test]
    fn cpu_pick_domain_rejects_square_corner_and_accepts_edge() {
        let wafer = WaferState::new_substrate(300.0, 500.0, 1_000_000.0);
        assert!(wafer.accepts_pick_xy(150_000_000.0, 0.0));
        assert!(wafer.accepts_pick_xy(0.0, 0.0));
        assert!(!wafer.accepts_pick_xy(149_000_000.0, 149_000_000.0));
        assert!(!wafer.accepts_pick_xy(f64::NAN, 0.0));
    }

}
