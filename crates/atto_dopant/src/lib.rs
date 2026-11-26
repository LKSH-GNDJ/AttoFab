// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The dopant field (build PRD D-01 to D-04).
//!
//! Geometry answers *where the material boundaries are*; a field answers *what
//! the value is inside them*. This one holds species, conductivity type and
//! concentration on the octree's leaf lattice, independently of whether the
//! octree has subdivided there — so an etch that later cuts a trench finds the
//! doping already at the right resolution.
//!
//! Net doping is derived, never stored: `N_net = N_D - N_A`. That is what makes
//! compensation work, and what makes a junction a sign change rather than a
//! material tag.

pub mod activation;
pub mod segregation;

use atto_units::{Nm};
use atto_geom::Run;
use atto_geom::svo::Material;
use atto_params::Species;
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Donor or acceptor.
pub enum Type {
    Donor,
    Acceptor,
}

impl Type {
    /// The name used in a report.
    pub fn label(self) -> &'static str {
        match self {
            Type::Donor => "n",
            Type::Acceptor => "p",
        }
    }
    /// +1 for a donor, -1 for an acceptor, for summing net doping.
    pub fn sign(self) -> f64 {
        match self {
            Type::Donor => 1.0,
            Type::Acceptor => -1.0,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// One implanted region: a species, its type and a concentration.
pub struct DopantRecord {
    pub species: Species,
    pub kind: Type,
    pub concentration_cm3: f64,
}

impl DopantRecord {
    /// Build a record.
    pub fn new(species: Species, kind: Type, concentration_cm3: f64) -> Self {
        DopantRecord { species, kind, concentration_cm3 }
    }
}

/// Integer lattice coordinate of a leaf. Stable under cloning and cheap to key.
pub type LeafKey = (i64, i64, i64);

#[derive(Debug, Clone, PartialEq)]
/// The dopant carried by the wafer: sparse cells plus a background.
pub struct DopantField {
    leaf_nm: f64,
    origin: (f64, f64, f64),
    /// Uniform doping present everywhere (the substrate background).
    background: Option<DopantRecord>,
    cells: BTreeMap<LeafKey, Vec<DopantRecord>>,
}

impl DopantField {
    /// Check the field's invariants, returning every violation (PRD A-11).
    ///
    /// An invariant that is never checked is a comment. These are the ones a
    /// corrupted or half-applied implant would break: a concentration that is
    /// negative or not a number, and a leaf size that cannot index anything.
    pub fn check(&self) -> Vec<String> {
        let mut out = Vec::new();
        if !(self.leaf_nm > 0.0) || !self.leaf_nm.is_finite() {
            out.push(format!("leaf size is not a positive number: {}", self.leaf_nm));
        }
        for (i, v) in [self.origin.0, self.origin.1, self.origin.2].iter().enumerate() {
            if !v.is_finite() {
                out.push(format!("origin axis {i} is not finite"));
            }
        }
        for (key, records) in &self.cells {
            for r in records {
                if !r.concentration_cm3.is_finite() || r.concentration_cm3 < 0.0 {
                    out.push(format!("cell {key:?} carries a concentration of {}", r.concentration_cm3));
                }
            }
        }
        if let Some(b) = &self.background {
            if !b.concentration_cm3.is_finite() || b.concentration_cm3 < 0.0 {
                out.push(format!("background concentration is {}", b.concentration_cm3));
            }
        }
        out
    }

    /// An empty field with a background and a leaf size.
    pub fn new(leaf: Nm, origin: (f64, f64, f64)) -> Self {
        let leaf_nm = leaf.value();
        DopantField { leaf_nm, origin, background: None, cells: BTreeMap::new() }
    }

    /// The cell size, nm.
    pub fn leaf_nm(&self) -> f64 { self.leaf_nm }
    /// How many cells carry an implant.
    pub fn occupied(&self) -> usize { self.cells.len() }
    pub fn background(&self) -> Option<DopantRecord> { self.background }

    /// The substrate background. A wafer is not undoped; it starts counter-doped.
    pub fn set_background(&mut self, rec: DopantRecord) { self.background = Some(rec); }

    /// The cell a point falls in.
    pub fn key_for(&self, x: f64, y: f64, z: f64) -> LeafKey {
        let f = |v: f64, o: f64| ((v - o) / self.leaf_nm).floor() as i64;
        (f(x, self.origin.0), f(y, self.origin.1), f(z, self.origin.2))
    }

    /// The centre of a cell.
    pub fn center_of(&self, key: LeafKey) -> (f64, f64, f64) {
        (
            self.origin.0 + (key.0 as f64 + 0.5) * self.leaf_nm,
            self.origin.1 + (key.1 as f64 + 0.5) * self.leaf_nm,
            self.origin.2 + (key.2 as f64 + 0.5) * self.leaf_nm,
        )
    }

    /// Add dopant to a leaf. Records of the same species and type accumulate,
    /// so two half-dose implants equal one full-dose implant (D-07).
    pub fn add(&mut self, key: LeafKey, rec: DopantRecord) {
        if !(rec.concentration_cm3.is_finite() && rec.concentration_cm3 > 0.0) {
            return;
        }
        let v = self.cells.entry(key).or_default();
        if let Some(existing) = v.iter_mut().find(|r| r.species == rec.species && r.kind == rec.kind) {
            existing.concentration_cm3 += rec.concentration_cm3;
        } else {
            v.push(rec);
        }
    }

    /// Replace the record for this species and type, or add it. Unlike `add`
    /// this does not accumulate, so a diffusion step can write its result back.
    /// A non-positive concentration removes the record.
    pub fn set(&mut self, key: LeafKey, rec: DopantRecord) {
        let v = self.cells.entry(key).or_default();
        match v.iter_mut().find(|r| r.species == rec.species && r.kind == rec.kind) {
            Some(existing) => existing.concentration_cm3 = rec.concentration_cm3,
            None => v.push(rec),
        }
        v.retain(|r| r.concentration_cm3 > 0.0 && r.concentration_cm3.is_finite());
        if v.is_empty() {
            self.cells.remove(&key);
        }
    }

    /// Every implant in a cell.
    pub fn records(&self, key: LeafKey) -> &[DopantRecord] {
        self.cells.get(&key).map(|v| v.as_slice()).unwrap_or(&[])
    }

    /// Signed net doping: donors positive, acceptors negative, cm^-3.
    pub fn net_cm3(&self, key: LeafKey) -> f64 {
        let mut net = 0.0;
        if let Some(b) = self.background {
            net += b.kind.sign() * b.concentration_cm3;
        }
        for r in self.records(key) {
            net += r.kind.sign() * r.concentration_cm3;
        }
        net
    }

    /// Net doping at a point, signed: positive is n-type.
    pub fn net_at(&self, x: f64, y: f64, z: f64) -> f64 {
        self.net_cm3(self.key_for(x, y, z))
    }

    /// Magnitude of the net doping, for reporting.
    pub fn concentration_at(&self, x: f64, y: f64, z: f64) -> f64 { self.net_at(x, y, z).abs() }

    /// The conductivity type at a point, if it is doped at all.
    pub fn type_at(&self, x: f64, y: f64, z: f64) -> Option<Type> {
        let n = self.net_at(x, y, z);
        if n > 0.0 {
            Some(Type::Donor)
        } else if n < 0.0 {
            Some(Type::Acceptor)
        } else {
            None
        }
    }

    /// The metallurgical junction in one column: the depth below `surface_z`
    /// where the net doping changes sign. Returns the shallowest such crossing.
    ///
    /// The column is walked at leaf resolution, so the result is quantised to
    /// the lattice — the same limit the geometry has.
    pub fn junction_nm(&self, column: &[Run], x: f64, y: f64, surface_z: f64) -> Option<f64> {
        let mut prev: Option<(f64, f64)> = None;
        for run in column.iter() {
            if run.material != Material::Silicon {
                continue;
            }
            let mut z = run.z1 - self.leaf_nm * 0.5;
            while z > run.z0 {
                let net = self.net_at(x, y, z);
                let depth = surface_z - z;
                if let Some((pdepth, pnet)) = prev {
                    if (pnet > 0.0) != (net > 0.0) && pnet != 0.0 && net != 0.0 {
                        // linear interpolation between the two sample depths
                        let f = pnet.abs() / (pnet.abs() + net.abs());
                        return Some(pdepth + f * (depth - pdepth));
                    }
                }
                prev = Some((depth, net));
                z -= self.leaf_nm;
            }
        }
        None
    }

    /// Every junction depth over a column grid, for a map.
    pub fn junction_map(&self, cols: &[Vec<Run>], active: &[bool], x0: f64, y0: f64, pitch: f64, nx: usize, surface_z: f64) -> Vec<Option<f64>> {
        let mut out = Vec::with_capacity(cols.len());
        for (k, col) in cols.iter().enumerate() {
            if !active[k] {
                out.push(None);
                continue;
            }
            let ix = k % nx;
            let iy = k / nx;
            let x = x0 + (ix as f64 + 0.5) * pitch;
            let y = y0 + (iy as f64 + 0.5) * pitch;
            out.push(self.junction_nm(col, x, y, surface_z));
        }
        out
    }

    /// Structural fingerprint of the field, for combining into the wafer's.
    pub fn fingerprint(&self) -> u64 {
        const PRIME: u64 = 0x0000_0100_0000_01B3;
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |bytes: &[u8]| {
            for b in bytes {
                h ^= *b as u64;
                h = h.wrapping_mul(PRIME);
            }
        };
        feed(&self.leaf_nm.to_bits().to_le_bytes());
        if let Some(b) = self.background {
            feed(&[1]);
            feed(&[b.species.label().len() as u8]);
            feed(&[b.kind as u8]);
            feed(&b.concentration_cm3.to_bits().to_le_bytes());
        }
        for (k, v) in self.cells.iter() {
            feed(&k.0.to_le_bytes());
            feed(&k.1.to_le_bytes());
            feed(&k.2.to_le_bytes());
            for r in v.iter() {
                feed(&[r.species.label().len() as u8, r.kind as u8]);
                feed(&r.concentration_cm3.to_bits().to_le_bytes());
            }
        }
        h
    }
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Um};
    use super::*;
    use atto_geom::column::ProcessWindow;
    use atto_geom::svo::Wafer;

    fn field() -> DopantField {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let mut f = DopantField::new(Nm::new(16.0), w.bounds.min);
        f.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e16));
        f
    }

    fn column(f: &DopantField) -> Vec<Run> {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let g = atto_geom::column::ColumnGrid::sample(&w, &ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] }).unwrap();
        g.cols[g.idx(8, 8)].clone()
    }

    #[test]
    fn keys_are_stable_and_local() {
        let f = field();
        let a = f.key_for(8.0, 8.0, 499_000.0);
        let b = f.key_for(8.5, 8.5, 499_000.0);
        let c = f.key_for(8.0, 8.0, 499_020.0);
        assert_eq!(a, b, "points in one leaf share a key");
        assert_ne!(a, c, "the next leaf up is a different key");
        let centre = f.center_of(a);
        assert!((centre.2 - 499_000.0).abs() <= f.leaf_nm());
    }

    #[test]
    fn net_doping_is_signed_and_compensates() {
        let mut f = field();
        let key = f.key_for(0.0, 0.0, 499_900.0);
        f.add(key, DopantRecord::new(Species::Phosphorus, Type::Donor, 5.0e16));
        // 5e16 donors against a 1e16 acceptor background -> +4e16 net donor
        assert!((f.net_cm3(key) - 4.0e16).abs() / 4.0e16 < 1e-12);
        assert_eq!(f.type_at(0.0, 0.0, 499_900.0), Some(Type::Donor));
        // deeper, undoped: still p from the background
        assert_eq!(f.type_at(0.0, 0.0, 400_000.0), Some(Type::Acceptor));
    }

    #[test]
    fn two_implants_accumulate() {
        let mut f = field();
        let key = f.key_for(0.0, 0.0, 499_900.0);
        f.add(key, DopantRecord::new(Species::Boron, Type::Acceptor, 2.0e17));
        f.add(key, DopantRecord::new(Species::Boron, Type::Acceptor, 2.0e17));
        assert_eq!(f.records(key).len(), 1, "same species and type merge");
        assert!((f.records(key)[0].concentration_cm3 - 4.0e17).abs() < 1.0);
    }

    #[test]
    fn opposite_types_stay_separate() {
        let mut f = field();
        let key = f.key_for(0.0, 0.0, 499_900.0);
        f.add(key, DopantRecord::new(Species::Boron, Type::Acceptor, 3.0e17));
        f.add(key, DopantRecord::new(Species::Phosphorus, Type::Donor, 1.0e17));
        assert_eq!(f.records(key).len(), 2);
        // net = 1e17 - 3e17 - 1e16 background = -2.1e17
        assert!((f.net_cm3(key) + 2.1e17).abs() / 2.1e17 < 1e-12);
    }

    #[test]
    fn a_junction_appears_where_the_sign_changes() {
        let mut f = field();
        // n+ region from the surface down to about 200 nm
        let surface = 500_000.0;
        for i in 0..12 {
            let z = surface - (i as f64 + 0.5) * 16.0;
            f.add(f.key_for(0.0, 0.0, z), DopantRecord::new(Species::Phosphorus, Type::Donor, 1.0e18));
        }
        let col = column(&f);
        let xj = f.junction_nm(&col, 0.0, 0.0, surface).expect("a junction");
        assert!((xj - 192.0).abs() <= 16.0, "junction at {xj} nm, expected about 192");
    }

    #[test]
    fn no_junction_when_the_doping_never_beats_the_background() {
        let mut f = field();
        let surface = 500_000.0;
        // 1e15 donors against a 1e16 acceptor background never flips the sign
        for i in 0..6 {
            let z = surface - (i as f64 + 0.5) * 16.0;
            f.add(f.key_for(0.0, 0.0, z), DopantRecord::new(Species::Phosphorus, Type::Donor, 1.0e15));
        }
        let col = column(&f);
        assert_eq!(f.junction_nm(&col, 0.0, 0.0, surface), None);
    }

    #[test]
    fn the_fingerprint_tracks_doping() {
        let a = field();
        let mut b = field();
        assert_eq!(a.fingerprint(), b.fingerprint());
        b.add(b.key_for(0.0, 0.0, 499_900.0), DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e17));
        assert_ne!(a.fingerprint(), b.fingerprint());
        let mut c = field();
        c.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, 2.0e16));
        assert_ne!(a.fingerprint(), c.fingerprint(), "the background is part of the state");
    }

    #[test]
    fn set_replaces_rather_than_accumulates() {
        let mut f = field();
        let key = f.key_for(0.0, 0.0, 499_900.0);
        f.add(key, DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e17));
        f.set(key, DopantRecord::new(Species::Boron, Type::Acceptor, 3.0e17));
        assert!((f.records(key)[0].concentration_cm3 - 3.0e17).abs() < 1.0);
        f.set(key, DopantRecord::new(Species::Boron, Type::Acceptor, 0.0));
        assert!(f.records(key).is_empty(), "a zero concentration removes the record");
        assert_eq!(f.occupied(), 0);
    }

    #[test]
    fn a_healthy_field_passes_its_own_check() {
        let mut f = DopantField::new(Nm::new(10.0), (0.0, 0.0, 0.0));
        let k = f.key_for(0.0, 0.0, 0.0);
        f.add(k, DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e18));
        assert!(f.check().is_empty(), "{:?}", f.check());
    }

    #[test]
    fn a_bad_concentration_cannot_be_stored_in_the_first_place() {
        // The check would catch a negative concentration, but `add` refuses to
        // store one, so the field never reaches that state. That is the stronger
        // property and it is the one asserted: prevention, not detection.
        let mut f = DopantField::new(Nm::new(10.0), (0.0, 0.0, 0.0));
        let k = f.key_for(0.0, 0.0, 0.0);
        f.add(k, DopantRecord::new(Species::Boron, Type::Acceptor, -1.0e18));
        assert!(f.check().is_empty(), "a negative was admitted: {:?}", f.check());
        assert!(f.net_at(0.0, 0.0, 0.0) >= 0.0, "and it did not become n-type");
    }

    #[test]
    fn a_field_with_a_broken_leaf_size_is_reported() {
        // this one the API cannot prevent: the leaf size is set at construction
        let mut g = DopantField::new(atto_units::Nm::new(f64::NAN), (0.0, 0.0, 0.0));
        let k = g.key_for(0.0, 0.0, 0.0);
        g.add(k, DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e18));
        let v = g.check();
        assert!(!v.is_empty(), "a NaN leaf size must be reported");
        assert!(v[0].contains("leaf size"), "{v:?}");
    }

    #[test]
    fn the_field_is_a_value() {
        let mut a = field();
        a.add(a.key_for(0.0, 0.0, 499_900.0), DopantRecord::new(Species::Arsenic, Type::Donor, 1.0e19));
        let b = a.clone();
        assert_eq!(a, b, "cloning the field gives an identical value");
        assert_eq!(a.fingerprint(), b.fingerprint());
    }
}
