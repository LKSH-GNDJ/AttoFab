// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The device cross-section: one row of the wafer, meshed for an electrical solve.

use atto_dopant::DopantField;
use atto_geom::column::{material_at_z, solid_top, ColumnGrid, ProcessWindow};
use atto_geom::svo::{Material, Wafer};
use atto_process::oxidation::InterfaceCharge;

#[derive(Debug, Clone, PartialEq)]
/// Why a cross-section could not be built.
pub enum CrossSectionError {
    Geometry(String),
    Empty(String),
    NoOxide(usize),
}

impl std::fmt::Display for CrossSectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CrossSectionError::Geometry(e) => write!(f, "cross-section geometry: {e}"),
            CrossSectionError::Empty(e) => write!(f, "cross-section has no silicon: {e}"),
            CrossSectionError::NoOxide(i) => write!(f, "column {i} carries no oxide, so it has no MOS interface"),
        }
    }
}
impl std::error::Error for CrossSectionError {}

/// One row of the wafer, meshed at the leaf pitch.
#[derive(Debug, Clone)]
pub struct CrossSection {
    pub x0: f64,
    pub y: f64,
    pub pitch: f64,
    pub nx: usize,
    pub nz: usize,
    pub z_bottom: f64,
    pub z_top: f64,
    /// Row-major `iz * nx + ix`, `iz = 0` at the bottom.
    pub material: Vec<Material>,
    /// Net doping at each cell centre, cm^-3 (donors positive).
    pub net_doping: Vec<f64>,
    /// Oxide thickness above the silicon, nm, per column.
    pub oxide_nm: Vec<f64>,
    /// Fixed charge at the silicon/oxide interface, cm^-2, per column.
    pub qf_cm2: Vec<f64>,
    /// Interface state density, cm^-2 eV^-1, per column.
    pub dit_cm2ev: Vec<f64>,
    /// Index of the first silicon cell from the top, per column.
    pub si_top: Vec<Option<usize>>,
}

impl CrossSection {
    /// Build a cross-section from the wafer, the doping and a window.
    pub fn build(w: &Wafer, field: &DopantField, win: &ProcessWindow, charge: InterfaceCharge) -> Result<Self, CrossSectionError> {
        let g = ColumnGrid::sample(w, win).map_err(CrossSectionError::Geometry)?;
        let row = g.ny / 2;
        let nz = ((g.z_top - g.z_bottom) / g.pitch).round() as usize;
        if nz == 0 {
            return Err(CrossSectionError::Empty("zero thickness".into()));
        }
        let nx = g.nx;
        let mut cs = CrossSection {
            x0: g.x0,
            y: g.y0 + (row as f64 + 0.5) * g.pitch,
            pitch: g.pitch,
            nx,
            nz,
            z_bottom: g.z_bottom,
            z_top: g.z_top,
            material: vec![Material::Air; nx * nz],
            net_doping: vec![0.0; nx * nz],
            oxide_nm: vec![0.0; nx],
            qf_cm2: vec![0.0; nx],
            dit_cm2ev: vec![0.0; nx],
            si_top: vec![None; nx],
        };
        let mut any_si = false;
        for ix in 0..nx {
            let k = g.idx(ix, row);
            if !g.active[k] {
                continue;
            }
            let col = &g.cols[k];
            let x = g.x0 + (ix as f64 + 0.5) * g.pitch;
            for iz in 0..nz {
                let z = g.z_bottom + (iz as f64 + 0.5) * g.pitch;
                let m = material_at_z(col, z);
                let idx = iz * nx + ix;
                cs.material[idx] = m;
                if m == Material::Silicon {
                    cs.net_doping[idx] = field.net_at(x, cs.y, z);
                    // the loop runs bottom-up, so overwriting keeps the topmost
                    cs.si_top[ix] = Some(iz);
                    any_si = true;
                }
            }
            // oxide thickness and interface charge, where the stack has an oxide
            let si_z = solid_top(col).map(|(t, m)| if m == Material::Silicon { Some(t) } else { None }).unwrap_or(None);
            let si_z = si_z.or_else(|| col.iter().rev().find(|r| r.material == Material::Silicon).map(|r| r.z1));
            if let Some(si) = si_z {
                let oxide: f64 = col
                    .iter()
                    .filter(|r| r.material == Material::SiliconDioxide && r.z0 >= si - 1e-9)
                    .map(|r| r.thickness())
                    .sum();
                cs.oxide_nm[ix] = oxide;
                if oxide > 0.0 {
                    cs.qf_cm2[ix] = charge.qf_cm2;
                    cs.dit_cm2ev[ix] = charge.dit_cm2ev;
                }
            }
        }
        if !any_si {
            return Err(CrossSectionError::Empty("no silicon in the window".into()));
        }
        Ok(cs)
    }

    /// Check the cross-section's invariants, returning every violation (PRD A-11).
    pub fn check(&self) -> Vec<String> {
        let mut out = Vec::new();
        if !(self.pitch > 0.0) || !self.pitch.is_finite() {
            out.push(format!("pitch is not a positive number: {}", self.pitch));
        }
        if self.nx == 0 || self.nz == 0 {
            out.push(format!("a {}x{} cross-section has no cells", self.nx, self.nz));
        }
        if !(self.z_top > self.z_bottom) {
            out.push(format!("the top {} is not above the bottom {}", self.z_top, self.z_bottom));
        }
        for (name, v) in [("x0", self.x0), ("y", self.y), ("z_bottom", self.z_bottom), ("z_top", self.z_top)] {
            if !v.is_finite() {
                out.push(format!("{name} is not finite"));
            }
        }
        out
    }

    /// The flat index of a cell.
    pub fn at(&self, ix: usize, iz: usize) -> usize { iz * self.nx + ix }

    /// The substrate doping a long-channel model wants: the net doping at the
    /// bottom of the silicon, well below any implant. Sampling a few cells under
    /// the surface instead reads the implant and gives a threshold voltage in
    /// the tens of volts.
    pub fn bulk_doping(&self, ix: usize) -> Option<f64> {
        let top = self.si_top[ix]?;
        let mut found = None;
        for k in 0..=top {
            if self.material[self.at(ix, k)] == Material::Silicon {
                found = Some(self.net_doping[self.at(ix, k)]);
                break; // the first silicon above the bottom is the deepest one
            }
        }
        found.filter(|v| v.abs() > 0.0)
    }
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Nm, Um};
    use super::*;
    use atto_dopant::{DopantRecord, Type};
    use atto_geom::column::grow_film;
    use atto_params::Species;

    fn setup() -> (Wafer, DopantField, ColumnGrid) {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let g = ColumnGrid::sample(&w, &win).unwrap();
        let mut f = DopantField::new(atto_units::Nm::new(16.0), w.bounds.min);
        f.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e16));
        (w, f, g)
    }

    #[test]
    fn the_cross_section_mirrors_the_column_stack() {
        let (mut w, f, g) = setup();
        let (ox, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(32.0), 0.0).unwrap();
        ox.write_back(&g, &mut w);
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let cs = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        assert!(cs.nx > 0 && cs.nz > 0);
        // the top cells are oxide, the cells below are silicon
        let ix = cs.nx / 2;
        let top = (0..cs.nz).rev().find(|iz| cs.material[cs.at(ix, *iz)] != Material::Air).expect("some solid");
        assert_eq!(cs.material[cs.at(ix, top)], Material::SiliconDioxide, "the topmost solid is the oxide");
        let si = cs.si_top[ix].expect("silicon present");
        assert_eq!(cs.material[cs.at(ix, si)], Material::Silicon);
        assert!((cs.oxide_nm[ix] - 32.0).abs() < 1e-9, "oxide {}", cs.oxide_nm[ix]);
    }

    #[test]
    fn a_healthy_cross_section_passes_its_own_check() {
        let (w, f, _g) = setup();
        let win = ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] };
        let cs = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        assert!(cs.check().is_empty(), "{:?}", cs.check());
    }

    #[test]
    fn a_cross_section_with_an_impossible_extent_is_reported() {
        let (w, f, _g) = setup();
        let win = ProcessWindow { center: [0.0, 0.0], size: [256.0, 256.0] };
        let mut cs = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        std::mem::swap(&mut cs.z_top, &mut cs.z_bottom);
        assert!(!cs.check().is_empty(), "a top below the bottom must be reported");
        let mut bad = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        bad.pitch = 0.0;
        assert!(bad.check().iter().any(|m| m.contains("pitch")), "{:?}", bad.check());
    }

    #[test]
    fn net_doping_reaches_the_mesh() {
        let (w, mut f, _g) = setup();
        // a donor layer just under the surface
        let key = f.key_for(0.0, 0.0, 499_990.0);
        f.add(key, DopantRecord::new(Species::Phosphorus, Type::Donor, 1.0e18));
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let cs = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        let ix = cs.nx / 2;
        let si = cs.si_top[ix].unwrap();
        let near = cs.net_doping[cs.at(ix, si)];
        assert!(near > 0.0, "the implanted layer is n-type on the mesh: {near}");
        let deep = cs.net_doping[cs.at(ix, si.saturating_sub(20))];
        assert!(deep < 0.0, "the substrate stays p-type: {deep}");
    }

    #[test]
    fn interface_charge_appears_only_where_there_is_oxide() {
        let (mut w, f, g) = setup();
        // oxide on the left half only
        let mut half = g.clone();
        for iy in 0..half.ny {
            for ix in 0..half.nx / 2 {
                let k = half.idx(ix, iy);
                let top = solid_top(&half.cols[k]).unwrap().0;
                atto_geom::column::set_interval(&mut half.cols[k], top, top + 32.0, Material::SiliconDioxide);
            }
        }
        half.write_back(&g, &mut w);
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let charge = InterfaceCharge { qf_cm2: 2.0e11, dit_cm2ev: 1.0e10 };
        let cs = CrossSection::build(&w, &f, &win, charge).unwrap();
        assert!((cs.qf_cm2[0] - 2.0e11).abs() < 1.0, "oxide columns carry the charge");
        let bare = cs.nx - 1;
        assert_eq!(cs.oxide_nm[bare], 0.0);
        assert_eq!(cs.qf_cm2[bare], 0.0, "no oxide, no interface charge");
    }

    #[test]
    fn bulk_doping_reads_the_substrate_not_the_implant() {
        // An n+ layer at the surface must not be mistaken for the substrate.
        let (w, mut f, _g) = setup();
        let key = f.key_for(0.0, 0.0, 499_990.0);
        f.add(key, DopantRecord::new(Species::Phosphorus, Type::Donor, 1.0e19));
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let cs = CrossSection::build(&w, &f, &win, InterfaceCharge::default()).unwrap();
        let b = cs.bulk_doping(cs.nx / 2).unwrap();
        assert!(b < 0.0, "the bulk is p-type even under an n+ layer: {b}");
    }
}
