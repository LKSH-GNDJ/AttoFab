//! Wafer2d: the 2D wafer cross-section grid (Phase 2).
//!
//! A dense, row-major flattened grid: `idx(x, y) = y * nx + x`. `x` runs
//! laterally across the wafer (columns), `y` runs with depth (0 = surface,
//! increasing = deeper into bulk silicon) - same depth convention as
//! `Wafer1d`, which this module treats as "one column" of a larger grid.
//!
//! Deliberately still a plain dense array (matching the "start dense,
//! prove correctness" plan): a 2D cross-section at reasonable resolution
//! (hundreds x hundreds of voxels) is small enough that RLE/chunked
//! sparsity isn't needed yet. That optimization is scoped for Phase 3
//! (whole-wafer 3D), once this dense 2D solver is validated - see
//! docs/architecture.md.

use crate::materials::{diffusion_params, Dopant, Material, SI_CONSUMPTION_RATIO};
use crate::oxidation::{grow_oxide, GrowOxideParams};
use crate::materials::{arrhenius, Ambient};
use serde::Serialize;
use std::collections::HashMap;

#[derive(Serialize)]

pub struct Wafer2d {
    pub nx: usize,
    pub ny: usize,
    pub dx_um: f64,
    pub dy_um: f64,
    pub material: Vec<Material>,       // len nx*ny, row-major (y*nx+x)
    pub species: HashMap<Dopant, Vec<f64>>, // len nx*ny each
    /// Oxide thickness is tracked per-column (x), since Deal-Grove is a 1D
    /// vertical process at each lateral position - this mirrors how a real
    /// LOCOS (local oxidation) process produces varying oxide thickness
    /// across x depending on which columns are exposed vs. masked.
    pub oxide_thickness_um: Vec<f64>, // len nx
    pub process_log: Vec<String>,
}

#[inline]
pub fn idx(x: usize, y: usize, nx: usize) -> usize {
    y * nx + x
}

impl Wafer2d {
    pub fn new(nx: usize, ny: usize, dx_um: f64, dy_um: f64, substrate_dopant: Dopant, substrate_conc: f64) -> Self {
        let n = nx * ny;
        let material = vec![Material::Silicon; n];
        let mut species = HashMap::new();
        species.insert(substrate_dopant, vec![substrate_conc; n]);

        Wafer2d {
            nx,
            ny,
            dx_um,
            dy_um,
            material,
            species,
            oxide_thickness_um: vec![0.0; nx],
            process_log: Vec::new(),
        }
    }

    #[inline]
    pub fn at(&self, x: usize, y: usize) -> usize {
        idx(x, y, self.nx)
    }

    fn ensure_species_channel(&mut self, dopant: Dopant) {
        let n = self.nx * self.ny;
        self.species.entry(dopant).or_insert_with(|| vec![0.0; n]);
    }

    /// Deposit a blanket layer of Photoresist across the top of every
    /// column (the first Void/topmost slot). Simplified: assumes a fixed
    /// resist thickness in voxel rows.
    pub fn spin_photoresist(&mut self, thickness_voxels: usize) {
        for x in 0..self.nx {
            for row in 0..thickness_voxels.min(self.ny) {
                let i = self.at(x, row);
                // Only coat where there's currently Silicon/Oxide surface
                // (don't overwrite existing resist or bury deeper layers
                // incorrectly) - Phase 2 keeps this intentionally simple:
                // resist always sits in the top N rows of a column.
                self.material[i] = Material::Photoresist;
            }
        }
        self.process_log.push(format!("spin_photoresist thickness={thickness_voxels}vox"));
    }

    /// Step 1: Thermal Oxidation (Deal-Grove), applied per-column. An
    /// optional `mask` (length nx, true = oxidize this column) supports
    /// LOCOS-style local oxidation where some columns are protected
    /// (e.g. by a nitride/resist layer) and don't grow oxide.
    pub fn oxidize(&mut self, temperature_c: f64, time_hours: f64, ambient: Ambient, mask: Option<&[bool]>) {
        for x in 0..self.nx {
            let allowed = mask.map(|m| m[x]).unwrap_or(true);
            if !allowed {
                continue;
            }
            let result = grow_oxide(&GrowOxideParams {
                temperature_c,
                time_hours,
                ambient,
                initial_oxide_um: self.oxide_thickness_um[x],
                orientation: crate::materials::CrystalOrientation::Si100,
            });
            let grown = result.final_oxide_um - self.oxide_thickness_um[x];
            self.oxide_thickness_um[x] = result.final_oxide_um;

            // Consume silicon rows proportional to grown oxide (converted
            // to a voxel count via dy_um), same 0.44x ratio as Wafer1d.
            let si_consumed_um = grown * SI_CONSUMPTION_RATIO;
            let rows_to_convert = (si_consumed_um / self.dy_um).round() as usize;
            let mut converted = 0;
            for y in 0..self.ny {
                if converted >= rows_to_convert {
                    break;
                }
                let i = self.at(x, y);
                if self.material[i] == Material::Silicon {
                    self.material[i] = Material::Oxide;
                    converted += 1;
                }
            }
        }
        self.process_log.push(format!(
            "oxidize T={temperature_c}C t={time_hours}hr ambient={ambient:?} masked_columns={}",
            mask.map(|m| m.iter().filter(|&&b| !b).count()).unwrap_or(0)
        ));
    }

    /// Step 2a: Ion Implantation (vertical Gaussian per column). An
    /// optional mask (length nx, true = implant this column) supports
    /// masked/selective doping. Phase 2 keeps the implant purely vertical
    /// (no lateral straggle) - lateral spread is instead captured
    /// naturally during the 2D anneal/diffusion step below.
    pub fn implant(&mut self, dopant: Dopant, dose_cm2: f64, energy_kev: f64, mask: Option<&[bool]>) {
        use crate::implant::{estimate_range, implant_profile, ImplantParams};
        self.ensure_species_channel(dopant);
        let range = estimate_range(dopant, energy_kev);

        let depth_grid_um: Vec<f64> = (0..self.ny).map(|y| y as f64 * self.dy_um).collect();
        let profile = implant_profile(&depth_grid_um, &ImplantParams { dose_cm2, range });

        let channel = self.species.get_mut(&dopant).unwrap();
        for x in 0..self.nx {
            let allowed = mask.map(|m| m[x]).unwrap_or(true);
            if !allowed {
                continue;
            }
            for y in 0..self.ny {
                let i = idx(x, y, self.nx);
                channel[i] += profile[y];
            }
        }
        self.process_log.push(format!("implant {dopant:?} dose={dose_cm2:.3e}/cm^2 energy={energy_kev}keV"));
    }

    /// Step 2b: Thermal Diffusion / Anneal - Fick's Second Law generalized
    /// to 2D:
    ///
    /// ```text
    /// dC/dt = D * (d^2C/dx^2 + d^2C/dy^2)
    /// ```
    ///
    /// Explicit FTCS with auto-selected stable dt. Stability condition
    /// (2D generalization of the 1D CFL-like limit):
    ///
    /// ```text
    /// D * dt * (1/dx^2 + 1/dy^2) <= 0.5
    /// ```
    ///
    /// Zero-flux (Neumann) boundaries on all four edges, same rationale
    /// as the 1D solver.
    pub fn anneal(&mut self, temperature_c: f64, time_minutes: f64, safety_factor: f64) {
        assert!(safety_factor < 0.5, "safety_factor must be < 0.5 for explicit FDM stability");

        let dx_cm = self.dx_um * 1e-4;
        let dy_cm = self.dy_um * 1e-4;
        let total_time_s = time_minutes * 60.0;
        let nx = self.nx;
        let ny = self.ny;

        let dopants: Vec<Dopant> = self.species.keys().copied().collect();
        for dopant in dopants {
            let params = diffusion_params(dopant);
            let t_kelvin = temperature_c + 273.15;
            let d_cm2_s = arrhenius(params.d0, params.ea, t_kelvin);

            let denom = 1.0 / (dx_cm * dx_cm) + 1.0 / (dy_cm * dy_cm);
            let dt_s = safety_factor / (d_cm2_s * denom);
            let n_steps = ((total_time_s / dt_s).ceil() as u64).max(1);
            let actual_dt = total_time_s / n_steps as f64;
            let lambda_x = d_cm2_s * actual_dt / (dx_cm * dx_cm);
            let lambda_y = d_cm2_s * actual_dt / (dy_cm * dy_cm);

            let mut c = self.species.get(&dopant).unwrap().clone();
            for _ in 0..n_steps {
                let mut c_next = c.clone();
                for y in 0..ny {
                    for x in 0..nx {
                        let i = idx(x, y, nx);
                        let left = if x == 0 { c[i] } else { c[idx(x - 1, y, nx)] };
                        let right = if x == nx - 1 { c[i] } else { c[idx(x + 1, y, nx)] };
                        let up = if y == 0 { c[i] } else { c[idx(x, y - 1, nx)] };
                        let down = if y == ny - 1 { c[i] } else { c[idx(x, y + 1, nx)] };
                        c_next[i] = c[i]
                            + lambda_x * (right - 2.0 * c[i] + left)
                            + lambda_y * (down - 2.0 * c[i] + up);
                    }
                }
                c = c_next;
            }
            self.species.insert(dopant, c);
        }

        self.process_log.push(format!("anneal T={temperature_c}C t={time_minutes}min"));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::Dopant;

    #[test]
    fn new_wafer_is_uniform_silicon_at_substrate_concentration() {
        let w = Wafer2d::new(10, 20, 0.01, 0.01, Dopant::Boron, 1e15);
        assert!(w.material.iter().all(|&m| m == Material::Silicon));
        let boron = &w.species[&Dopant::Boron];
        assert!(boron.iter().all(|&c| (c - 1e15).abs() < 1.0));
    }

    #[test]
    fn masked_oxidation_only_grows_on_open_columns() {
        let mut w = Wafer2d::new(5, 50, 0.005, 0.005, Dopant::Boron, 1e15);
        // Only column 2 is exposed (LOCOS-style local oxidation).
        let mask = [false, false, true, false, false];
        w.oxidize(1000.0, 0.75, Ambient::Dry, Some(&mask));

        assert!(w.oxide_thickness_um[2] > 0.0, "masked-open column should grow oxide");
        for &x in &[0usize, 1, 3, 4] {
            assert_eq!(w.oxide_thickness_um[x], 0.0, "masked-off column {x} must not grow oxide");
        }
    }

    #[test]
    fn diffusion_2d_spreads_symmetrically_from_a_point_source() {
        // A localized dopant spike at the grid center should diffuse
        // outward symmetrically in x given a symmetric grid and BCs.
        let nx = 21;
        let ny = 21;
        let mut w = Wafer2d::new(nx, ny, 0.01, 0.01, Dopant::Boron, 0.0);
        w.species.insert(Dopant::Phosphorus, vec![0.0; nx * ny]);
        let center = idx(nx / 2, ny / 2, nx);
        w.species.get_mut(&Dopant::Phosphorus).unwrap()[center] = 1e20;

        w.anneal(1000.0, 5.0, 0.2);

        let p = &w.species[&Dopant::Phosphorus];
        let left = p[idx(nx / 2 - 3, ny / 2, nx)];
        let right = p[idx(nx / 2 + 3, ny / 2, nx)];
        let up = p[idx(nx / 2, ny / 2 - 3, nx)];
        let down = p[idx(nx / 2, ny / 2 + 3, nx)];

        assert!(left > 0.0 && right > 0.0, "dopant should have spread laterally");
        assert!(
            (left - right).abs() / left < 1e-6,
            "spread should be symmetric left/right: left={left} right={right}"
        );
        assert!(
            (up - down).abs() / up < 1e-6,
            "spread should be symmetric up/down: up={up} down={down}"
        );
    }

    #[test]
    fn diffusion_2d_conserves_total_dose() {
        let nx = 15;
        let ny = 15;
        let mut w = Wafer2d::new(nx, ny, 0.01, 0.01, Dopant::Boron, 0.0);
        w.species.insert(Dopant::Phosphorus, vec![0.0; nx * ny]);
        let center = idx(nx / 2, ny / 2, nx);
        w.species.get_mut(&Dopant::Phosphorus).unwrap()[center] = 1e20;

        let before: f64 = w.species[&Dopant::Phosphorus].iter().sum();
        w.anneal(1000.0, 2.0, 0.2);
        let after: f64 = w.species[&Dopant::Phosphorus].iter().sum();

        // Zero-flux boundaries => total mass conserved (to numerical
        // precision), as long as the dopant hasn't hit the grid edges.
        let relative_error = (after - before).abs() / before;
        assert!(relative_error < 1e-6, "mass should be conserved, error={relative_error}");
    }
}
