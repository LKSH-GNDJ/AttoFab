// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The sparse voxel octree (ported from the earlier build, with the audit's fixes).

use atto_units::{Mm, Nm, Um};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
/// The materials the geometry knows about.
pub enum Material {
    Air,
    Silicon,
    SiliconDioxide,
    SiliconNitride,
    Photoresist,
    Metal,
}

impl Material {
    /// Whether the material occupies space.
    pub fn is_solid(self) -> bool { self != Material::Air }
    /// The name used in a report.
    pub fn label(self) -> &'static str {
        match self {
            Material::Air => "Air",
            Material::Silicon => "Si",
            Material::SiliconDioxide => "SiO2",
            Material::SiliconNitride => "Si3N4",
            Material::Photoresist => "resist",
            Material::Metal => "metal",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// An axis-aligned box.
pub struct Bounds {
    pub min: (f64, f64, f64),
    pub max: (f64, f64, f64),
}

impl Bounds {
    /// A box from two corners.
    pub fn new(min: (f64, f64, f64), max: (f64, f64, f64)) -> Self { Bounds { min, max } }
    /// The centre of the box.
    pub fn center(&self) -> (f64, f64, f64) {
        ((self.min.0 + self.max.0) * 0.5, (self.min.1 + self.max.1) * 0.5, (self.min.2 + self.max.2) * 0.5)
    }
    /// The box edge length.
    pub fn edge(&self) -> f64 { self.max.0 - self.min.0 }
    /// The i-th child box.
    pub fn octant(&self, i: usize) -> Bounds {
        let c = self.center();
        let (x0, x1) = if i & 1 == 0 { (self.min.0, c.0) } else { (c.0, self.max.0) };
        let (y0, y1) = if i & 2 == 0 { (self.min.1, c.1) } else { (c.1, self.max.1) };
        let (z0, z1) = if i & 4 == 0 { (self.min.2, c.2) } else { (c.2, self.max.2) };
        Bounds { min: (x0, y0, z0), max: (x1, y1, z1) }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// A contiguous stretch of one material in a column.
pub struct Run {
    pub z0: f64,
    pub z1: f64,
    pub material: Material,
}

impl Run {
    /// The run's height, nm.
    pub fn thickness(&self) -> f64 { self.z1 - self.z0 }
}

#[derive(Debug, Clone)]
/// A sparse voxel node: a uniform leaf or eight children.
pub enum SvoNode {
    Leaf(Material),
    Internal(Box<[SvoNode; 8]>),
}

fn overlaps(a: &Bounds, b: &Bounds) -> bool {
    a.min.0 < b.max.0 && a.max.0 > b.min.0
        && a.min.1 < b.max.1 && a.max.1 > b.min.1
        && a.min.2 < b.max.2 && a.max.2 > b.min.2
}

fn contains(outer: &Bounds, inner: &Bounds) -> bool {
    outer.min.0 <= inner.min.0 && outer.max.0 >= inner.max.0
        && outer.min.1 <= inner.min.1 && outer.max.1 >= inner.max.1
        && outer.min.2 <= inner.min.2 && outer.max.2 >= inner.max.2
}

fn collapse_if_homogeneous(node: &mut SvoNode) {
    if let SvoNode::Internal(children) = node {
        if let SvoNode::Leaf(first) = &children[0] {
            let first = *first;
            if children.iter().all(|c| matches!(c, SvoNode::Leaf(m) if *m == first)) {
                *node = SvoNode::Leaf(first);
            }
        }
    }
}

#[derive(Debug, Clone)]
/// A sparse voxel wafer: the geometry everything else reads.
pub struct Wafer {
    pub root: SvoNode,
    pub bounds: Bounds,
    pub radius_nm: f64,
    pub leaf_nm: f64,
    /// Z of the polished substrate surface.
    pub surface_z_nm: f64,
    /// True when the octree is a high-resolution tile and everything outside it
    /// is implicit bulk (silicon below the surface, air above).
    pub implicit_bulk: bool,
}

impl Wafer {
    /// A cubic tile of `leaf_nm * 2^levels` centred on `site`, with the substrate
    /// surface on the tile's mid-height plane.
    pub fn new_tile(diameter: Mm, thickness: Um, leaf: Nm, levels: u32, site: (f64, f64)) -> Result<Self, String> {
        // A-03: three lengths, three units, made explicit rather than left to the
        // parameter names. The values are unchanged — only the type says which.
        let (diameter_mm, thickness_um, leaf_nm) = (diameter.value(), thickness.value(), leaf.value());
        if !(diameter_mm.is_finite() && diameter_mm > 0.0) { return Err("diameter must be finite and positive".into()); }
        if !(thickness_um.is_finite() && thickness_um > 0.0) { return Err("thickness must be finite and positive".into()); }
        if !(leaf_nm.is_finite() && leaf_nm > 0.0) { return Err("leaf size must be finite and positive".into()); }
        let levels = levels.clamp(1, 14);
        let radius_nm = diameter_mm * 1e6 * 0.5;
        let t = thickness_um * 1e3;
        let edge = leaf_nm * (1u64 << levels) as f64;
        let h = edge * 0.5;
        if edge > t * 2.0 { return Err("tile is thicker than the substrate".into()); }
        if site.0 * site.0 + site.1 * site.1 > radius_nm * radius_nm { return Err("tile centre is outside the wafer".into()); }
        let bounds = Bounds::new((site.0 - h, site.1 - h, t - h), (site.0 + h, site.1 + h, t + h));
        let children: [SvoNode; 8] = std::array::from_fn(|i| {
            if i & 4 == 0 { SvoNode::Leaf(Material::Silicon) } else { SvoNode::Leaf(Material::Air) }
        });
        Ok(Wafer { root: SvoNode::Internal(Box::new(children)), bounds, radius_nm, leaf_nm, surface_z_nm: t, implicit_bulk: true })
    }

    /// Whether a point lies on the wafer disk.
    pub fn contains_xy(&self, x: f64, y: f64) -> bool {
        if !(x.is_finite() && y.is_finite()) { return false; }
        x * x + y * y <= self.radius_nm * self.radius_nm
    }

    /// The true finest leaf edge: the root halved until it is no larger than
    /// `leaf_nm`, which is exactly what the subdivision test uses (audit F-15).
    pub fn leaf_pitch_nm(&self) -> f64 {
        let mut edge = self.bounds.edge();
        let mut guard = 0;
        while edge > self.leaf_nm && guard < 64 {
            edge *= 0.5;
            guard += 1;
        }
        edge
    }

    /// What occupies a point, or None outside the wafer.
    pub fn material_at(&self, x: f64, y: f64, z: f64) -> Option<Material> {
        if !(x.is_finite() && y.is_finite() && z.is_finite()) { return None; }
        if !self.contains_xy(x, y) { return None; }
        let b = self.bounds;
        let inside = x >= b.min.0 && x <= b.max.0 && y >= b.min.1 && y <= b.max.1 && z >= b.min.2 && z <= b.max.2;
        if !inside {
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
                SvoNode::Internal(ch) => {
                    let c = nb.center();
                    let idx = (x >= c.0) as usize | (((y >= c.1) as usize) << 1) | (((z >= c.2) as usize) << 2);
                    nb = nb.octant(idx);
                    node = &ch[idx];
                }
            }
        }
    }

    /// The vertical material stack at (x, y), bottom to top, adjacent equal
    /// materials merged.
    pub fn column_runs(&self, x: f64, y: f64) -> Vec<Run> {
        let b = self.bounds;
        let mut out: Vec<Run> = Vec::new();
        if !(x.is_finite() && y.is_finite()) { return out; }
        if x < b.min.0 || x > b.max.0 || y < b.min.1 || y > b.max.1 { return out; }
        fn walk(node: &SvoNode, nb: Bounds, x: f64, y: f64, out: &mut Vec<Run>) {
            match node {
                SvoNode::Leaf(m) => {
                    if let Some(last) = out.last_mut() {
                        if last.material == *m && (last.z1 - nb.min.2).abs() < 1e-9 {
                            last.z1 = nb.max.2;
                            return;
                        }
                    }
                    out.push(Run { z0: nb.min.2, z1: nb.max.2, material: *m });
                }
                SvoNode::Internal(ch) => {
                    let c = nb.center();
                    let ix = (x >= c.0) as usize | (((y >= c.1) as usize) << 1);
                    for zb in 0..2 {
                        let idx = ix | (zb << 2);
                        walk(&ch[idx], nb.octant(idx), x, y, out);
                    }
                }
            }
        }
        walk(&self.root, b, x, y, &mut out);
        out
    }

    /// Fill a region, subdividing only where it reaches.
    pub fn set_region(&mut self, region: Bounds, material: Material) {
        Self::rec(&mut self.root, self.bounds, region, material, self.leaf_nm, self.radius_nm);
    }

    fn rec(node: &mut SvoNode, nb: Bounds, region: Bounds, mat: Material, leaf: f64, radius: f64) {
        if !overlaps(&nb, &region) { return; }
        let c = nb.center();
        if c.0 * c.0 + c.1 * c.1 > radius * radius && nb.edge() < leaf * 64.0 { return; }
        if contains(&region, &nb) {
            *node = SvoNode::Leaf(mat);
            return;
        }
        let edge = nb.edge();
        if edge <= leaf {
            if c.0 >= region.min.0 && c.0 <= region.max.0
                && c.1 >= region.min.1 && c.1 <= region.max.1
                && c.2 >= region.min.2 && c.2 <= region.max.2
                && c.0 * c.0 + c.1 * c.1 <= radius * radius
            {
                *node = SvoNode::Leaf(mat);
            }
            return;
        }
        if let SvoNode::Leaf(existing) = node {
            let e = *existing;
            *node = SvoNode::Internal(Box::new(std::array::from_fn(|_| SvoNode::Leaf(e))));
        }
        if let SvoNode::Internal(ch) = node {
            for i in 0..8 {
                Self::rec(&mut ch[i], nb.octant(i), region, mat, leaf, radius);
            }
        }
        collapse_if_homogeneous(node);
    }

    /// Subdivide leaves where `test` says the geometry error is large (PRD W-06).
    ///
    /// The tree is already adaptive where geometry *exists*: `set_region` only
    /// descends into subtrees a region touches. This is the other half — refining
    /// where the estimate says the answer is uncertain, which geometry alone
    /// cannot know. A leaf that fails the test is left alone, and a subdivision
    /// that would produce one material throughout is collapsed again.
    ///
    /// Returns how many leaves were split.
    pub fn refine_where(&mut self, levels: u32, test: impl Fn(Bounds) -> bool) -> usize {
        let mut split = 0;
        Self::refine_rec(&mut self.root, self.bounds, levels, &test, &mut split);
        split
    }

    fn refine_rec(
        node: &mut SvoNode,
        nb: Bounds,
        levels: u32,
        test: &impl Fn(Bounds) -> bool,
        split: &mut usize,
    ) {
        if levels == 0 {
            return;
        }
        match node {
            SvoNode::Leaf(m) => {
                if !test(nb) {
                    return;
                }
                let m = *m;
                let mut kids: Box<[SvoNode; 8]> =
                    Box::new([SvoNode::Leaf(m), SvoNode::Leaf(m), SvoNode::Leaf(m), SvoNode::Leaf(m),
                              SvoNode::Leaf(m), SvoNode::Leaf(m), SvoNode::Leaf(m), SvoNode::Leaf(m)]);
                for (i, kid) in kids.iter_mut().enumerate() {
                    Self::refine_rec(kid, nb.octant(i), levels - 1, test, split);
                }
                *node = SvoNode::Internal(kids);
                *split += 1;
            }
            SvoNode::Internal(kids) => {
                for (i, kid) in kids.iter_mut().enumerate() {
                    Self::refine_rec(kid, nb.octant(i), levels, test, split);
                }
            }
        }
    }

    /// Walk every leaf with its bounds and material.
    pub fn visit_leaves(&self, mut f: impl FnMut(Bounds, Material)) {
        fn walk(n: &SvoNode, b: Bounds, f: &mut impl FnMut(Bounds, Material)) {
            match n {
                SvoNode::Leaf(m) => f(b, *m),
                SvoNode::Internal(ch) => {
                    for (i, c) in ch.iter().enumerate() {
                        walk(c, b.octant(i), f);
                    }
                }
            }
        }
        walk(&self.root, self.bounds, &mut f);
    }

    /// How many nodes the tree holds, for a refinement report.
    pub fn node_count(&self) -> usize {
        fn walk(n: &SvoNode) -> usize {
            match n {
                SvoNode::Leaf(_) => 1,
                SvoNode::Internal(ch) => 1 + ch.iter().map(walk).sum::<usize>(),
            }
        }
        walk(&self.root)
    }

    /// Structural fingerprint: geometry plus surface state (audit F-17).
    pub fn fingerprint(&self) -> u64 {
        const PRIME: u64 = 0x0000_0100_0000_01B3;
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |bytes: &[u8]| {
            for b in bytes {
                h ^= *b as u64;
                h = h.wrapping_mul(PRIME);
            }
        };
        for v in [self.bounds.min.0, self.bounds.min.1, self.bounds.min.2, self.bounds.max.0, self.bounds.max.1, self.bounds.max.2, self.surface_z_nm, self.leaf_nm, self.radius_nm] {
            feed(&v.to_bits().to_le_bytes());
        }
        fn walk(n: &SvoNode, feed: &mut impl FnMut(&[u8])) {
            match n {
                SvoNode::Leaf(m) => {
                    feed(&[0]);
                    feed(&[*m as u8]);
                }
                SvoNode::Internal(ch) => {
                    feed(&[1]);
                    for c in ch.iter() {
                        walk(c, feed);
                    }
                }
            }
        }
        walk(&self.root, &mut feed);
        h
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tile() -> Wafer { Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap() }

    #[test]
    fn a_fresh_tile_is_nine_nodes() {
        let w = tile();
        assert_eq!(w.node_count(), 9, "root plus eight octants");
        assert_eq!(w.leaf_pitch_nm(), 16.0);
    }

    #[test]
    fn a_local_edit_subdivides_and_collapses_back() {
        let mut w = tile();
        let before = w.node_count();
        w.set_region(Bounds::new((-32.0, -32.0, 499_000.0), (32.0, 32.0, 500_000.0)), Material::SiliconDioxide);
        let after = w.node_count();
        assert!(after > before, "the edit subdivides: {before} -> {after}");
        assert_eq!(w.material_at(0.0, 0.0, 499_900.0), Some(Material::SiliconDioxide));
        assert_eq!(w.material_at(200.0, 200.0, 499_900.0), Some(Material::Silicon), "the edit is local");
        w.set_region(w.bounds, Material::Silicon);
        assert_eq!(w.node_count(), 1, "a homogeneous tree collapses to one leaf");
    }

    #[test]
    fn refinement_only_touches_where_the_test_says_yes() {
        let mut w = tile();
        let before = w.node_count();
        let split = w.refine_where(2, |b| b.min.0 < -100.0 && b.min.1 < -100.0);
        assert!(split > 0, "something must have split");
        assert!(w.node_count() > before, "the tree grew: {before} -> {}", w.node_count());
        let n = w.node_count();
        assert_eq!(w.refine_where(2, |_| false), 0);
        assert_eq!(w.node_count(), n, "a predicate that never fires must not split");
    }

    #[test]
    fn refinement_does_not_change_the_material_anywhere() {
        // refining is a representation change: the geometry must be identical
        // before and after, or it is an edit rather than a refinement
        let mut w = tile();
        w.set_region(Bounds::new((-50.0, -50.0, 0.0), (50.0, 50.0, 200.0)), Material::Silicon);
        let probe: Vec<f64> = vec![-80.0, -20.0, 0.0, 20.0, 80.0];
        let before: Vec<Option<Material>> = probe.iter().map(|x| w.material_at(*x, 0.0, 100.0)).collect();
        w.refine_where(2, |b| b.max.2 > 0.0);
        let after: Vec<Option<Material>> = probe.iter().map(|x| w.material_at(*x, 0.0, 100.0)).collect();
        assert_eq!(before, after, "refinement changed the answer");
    }

    #[test]
    fn column_runs_are_contiguous_bottom_to_top() {
        let mut w = tile();
        w.set_region(Bounds::new((-8.0, -8.0, 500_000.0), (8.0, 8.0, 500_048.0)), Material::SiliconDioxide);
        let col = w.column_runs(0.0, 0.0);
        assert!(col.len() >= 2, "{col:?}");
        for pair in col.windows(2) {
            assert!((pair[0].z1 - pair[1].z0).abs() < 1e-9, "runs must be contiguous: {pair:?}");
        }
        assert_eq!(col.first().unwrap().z0, w.bounds.min.2);
        assert_eq!(col.last().unwrap().z1, w.bounds.max.2);
    }

    #[test]
    fn outside_the_disk_is_not_material() {
        let w = tile();
        assert!(w.contains_xy(0.0, 0.0));
        assert!(!w.contains_xy(149e6, 149e6));
        assert_eq!(w.material_at(149e6, 149e6, 499_000.0), None);
    }

    #[test]
    fn fingerprint_tracks_geometry() {
        let a = tile();
        let mut b = tile();
        assert_eq!(a.fingerprint(), b.fingerprint());
        b.set_region(Bounds::new((-8.0, -8.0, 500_000.0), (8.0, 8.0, 500_016.0)), Material::Metal);
        assert_ne!(a.fingerprint(), b.fingerprint());
    }
}
