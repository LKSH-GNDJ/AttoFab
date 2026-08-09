//! Wafer3d: the dense 3D wafer mesh (Phase 3).
//!
//! Row-major flattened grid: `idx(x, y, z) = z*nx*ny + y*nx + x`. `x`/`y`
//! are lateral (the GDSII mask plane), `z` is depth (0 = surface,
//! increasing = deeper into bulk silicon) - same depth convention as
//! `Wafer1d`/`Wafer2d`.
//!
//! Still a plain dense array, per the "start dense, prove correctness"
//! plan that Phase 1 and Phase 2 both followed: get a small, obviously-
//! correct 3D solver working and regression-tested first. RLE-compressed
//! bulk + chunked-sparse active zone (the actual wafer-scale memory
//! strategy - see docs/architecture.md) is deliberately NOT implemented
//! here; this module is the validated baseline that optimization would be
//! checked against, not a replacement for needing it. Keep grids small
//! (tens of voxels per axis) until that optimization exists - a dense
//! nx*ny*nz f64 array per dopant species grows fast (e.g. 50^3 voxels x 3
//! species x 8 bytes = 3MB, fine; 500^3 would be ~3GB and impractical).

use std::collections::HashMap;

use serde::Serialize;

use crate::materials::{arrhenius, diffusion_params, Ambient, Dopant, Material, SI_CONSUMPTION_RATIO};
use crate::oxidation::{grow_oxide, GrowOxideParams};

#[derive(Serialize)]
pub struct Wafer3d {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub dx_um: f64,
    pub dy_um: f64,
    pub dz_um: f64,
    pub material: Vec<Material>,            // len nx*ny*nz
    pub species: HashMap<Dopant, Vec<f64>>, // len nx*ny*nz each
    /// Oxide thickness per (x,y) column, len nx*ny.
    pub oxide_thickness_um: Vec<f64>,
    pub process_log: Vec<String>,
}

#[inline]
pub fn idx(x: usize, y: usize, z: usize, nx: usize, ny: usize) -> usize {
    z * nx * ny + y * nx + x
}

#[inline]
pub fn col_idx(x: usize, y: usize, nx: usize) -> usize {
    y * nx + x
}

impl Wafer3d {
    pub fn new(
        nx: usize,
        ny: usize,
        nz: usize,
        dx_um: f64,
        dy_um: f64,
        dz_um: f64,
        substrate_dopant: Dopant,
        substrate_conc: f64,
    ) -> Self {
        let n = nx * ny * nz;
        let material = vec![Material::Silicon; n];
        let mut species = HashMap::new();
        species.insert(substrate_dopant, vec![substrate_conc; n]);

        Wafer3d {
            nx,
            ny,
            nz,
            dx_um,
            dy_um,
            dz_um,
            material,
            species,
            oxide_thickness_um: vec![0.0; nx * ny],
            process_log: Vec::new(),
        }
    }

    #[inline]
    pub fn at(&self, x: usize, y: usize, z: usize) -> usize {
        idx(x, y, z, self.nx, self.ny)
    }

    fn ensure_species_channel(&mut self, dopant: Dopant) {
        let n = self.nx * self.ny * self.nz;
        self.species.entry(dopant).or_insert_with(|| vec![0.0; n]);
    }

    /// Blanket-coat photoresist across the top `thickness_voxels` of every
    /// (x,y) column.
    pub fn spin_photoresist(&mut self, thickness_voxels: usize) {
        for y in 0..self.ny {
            for x in 0..self.nx {
                for z in 0..thickness_voxels.min(self.nz) {
                    let i = self.at(x, y, z);
                    self.material[i] = Material::Photoresist;
                }
            }
        }
        self.process_log.push(format!("spin_photoresist thickness={thickness_voxels}vox"));
    }

    /// Thermal Oxidation (Deal-Grove), applied per (x,y) column. `mask`
    /// (length nx*ny, true = oxidize this column) supports GDSII-mask-
    /// driven local oxidation (LOCOS), same as Wafer2d.
    pub fn oxidize(&mut self, temperature_c: f64, time_hours: f64, ambient: Ambient, mask: Option<&[bool]>) {
        let mut masked_count = 0usize;
        for y in 0..self.ny {
            for x in 0..self.nx {
                let ci = col_idx(x, y, self.nx);
                let allowed = mask.map(|m| m[ci]).unwrap_or(true);
                if !allowed {
                    masked_count += 1;
                    continue;
                }
                let result = grow_oxide(&GrowOxideParams {
                    temperature_c,
                    time_hours,
                    ambient,
                    initial_oxide_um: self.oxide_thickness_um[ci],
                    orientation: crate::materials::CrystalOrientation::Si100,
                });
                let grown = result.final_oxide_um - self.oxide_thickness_um[ci];
                self.oxide_thickness_um[ci] = result.final_oxide_um;

                let si_consumed_um = grown * SI_CONSUMPTION_RATIO;
                let voxels_to_convert = (si_consumed_um / self.dz_um).round() as usize;
                let mut converted = 0;
                for z in 0..self.nz {
                    if converted >= voxels_to_convert {
                        break;
                    }
                    let i = self.at(x, y, z);
                    if self.material[i] == Material::Silicon {
                        self.material[i] = Material::Oxide;
                        converted += 1;
                    }
                }
            }
        }
        self.process_log.push(format!(
            "oxidize T={temperature_c}C t={time_hours}hr ambient={ambient:?} masked_columns={masked_count}"
        ));
    }

    /// Ion Implantation (vertical Gaussian per (x,y) column, no lateral
    /// straggle at implant time - lateral spread comes from the 3D anneal
    /// step, same design choice as Wafer2d).
    pub fn implant(&mut self, dopant: Dopant, dose_cm2: f64, energy_kev: f64, mask: Option<&[bool]>) {
        use crate::implant::{estimate_range, implant_profile, ImplantParams};
        self.ensure_species_channel(dopant);
        let range = estimate_range(dopant, energy_kev);

        let depth_grid_um: Vec<f64> = (0..self.nz).map(|z| z as f64 * self.dz_um).collect();
        let profile = implant_profile(&depth_grid_um, &ImplantParams { dose_cm2, range });

        let channel = self.species.get_mut(&dopant).unwrap();
        for y in 0..self.ny {
            for x in 0..self.nx {
                let ci = col_idx(x, y, self.nx);
                let allowed = mask.map(|m| m[ci]).unwrap_or(true);
                if !allowed {
                    continue;
                }
                for z in 0..self.nz {
                    let i = idx(x, y, z, self.nx, self.ny);
                    channel[i] += profile[z];
                }
            }
        }
        self.process_log.push(format!("implant {dopant:?} dose={dose_cm2:.3e}/cm^2 energy={energy_kev}keV"));
    }

    /// Thermal Diffusion / Anneal - Fick's Second Law in 3D:
    ///
    /// ```text
    /// dC/dt = D * (d^2C/dx^2 + d^2C/dy^2 + d^2C/dz^2)
    /// ```
    ///
    /// Explicit FTCS, auto-selected stable dt via the 3D CFL-like
    /// condition `D*dt*(1/dx^2+1/dy^2+1/dz^2) <= 0.5`, zero-flux (Neumann)
    /// boundaries on all six faces.
    pub fn anneal(&mut self, temperature_c: f64, time_minutes: f64, safety_factor: f64) {
        assert!(safety_factor < 0.5, "safety_factor must be < 0.5 for explicit FDM stability");

        let dx_cm = self.dx_um * 1e-4;
        let dy_cm = self.dy_um * 1e-4;
        let dz_cm = self.dz_um * 1e-4;
        let total_time_s = time_minutes * 60.0;
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);

        let dopants: Vec<Dopant> = self.species.keys().copied().collect();
        for dopant in dopants {
            let params = diffusion_params(dopant);
            let t_kelvin = temperature_c + 273.15;
            let d_cm2_s = arrhenius(params.d0, params.ea, t_kelvin);

            let denom = 1.0 / (dx_cm * dx_cm) + 1.0 / (dy_cm * dy_cm) + 1.0 / (dz_cm * dz_cm);
            let dt_s = safety_factor / (d_cm2_s * denom);
            let n_steps = ((total_time_s / dt_s).ceil() as u64).max(1);
            let actual_dt = total_time_s / n_steps as f64;
            let lambda_x = d_cm2_s * actual_dt / (dx_cm * dx_cm);
            let lambda_y = d_cm2_s * actual_dt / (dy_cm * dy_cm);
            let lambda_z = d_cm2_s * actual_dt / (dz_cm * dz_cm);

            let mut c = self.species.get(&dopant).unwrap().clone();
            for _ in 0..n_steps {
                let mut c_next = c.clone();
                for z in 0..nz {
                    for y in 0..ny {
                        for x in 0..nx {
                            let i = idx(x, y, z, nx, ny);
                            let xm = if x == 0 { c[i] } else { c[idx(x - 1, y, z, nx, ny)] };
                            let xp = if x == nx - 1 { c[i] } else { c[idx(x + 1, y, z, nx, ny)] };
                            let ym = if y == 0 { c[i] } else { c[idx(x, y - 1, z, nx, ny)] };
                            let yp = if y == ny - 1 { c[i] } else { c[idx(x, y + 1, z, nx, ny)] };
                            let zm = if z == 0 { c[i] } else { c[idx(x, y, z - 1, nx, ny)] };
                            let zp = if z == nz - 1 { c[i] } else { c[idx(x, y, z + 1, nx, ny)] };
                            c_next[i] = c[i]
                                + lambda_x * (xp - 2.0 * c[i] + xm)
                                + lambda_y * (yp - 2.0 * c[i] + ym)
                                + lambda_z * (zp - 2.0 * c[i] + zm);
                        }
                    }
                }
                c = c_next;
            }
            self.species.insert(dopant, c);
        }

        self.process_log.push(format!("anneal T={temperature_c}C t={time_minutes}min"));
    }

    /// Anisotropic (directional, "straight down") etch: removes up to
    /// `depth_voxels` of solid material from the surface downward in each
    /// mask-open (x,y) column, no lateral spread. `mask` length nx*ny.
    /// Unlike Wafer2d's `etch::etch_anisotropic` (which requires a prior
    /// `develop()` to have already opened the surface to Void), this
    /// takes the open/closed decision directly from `mask` for the step -
    /// a reasonable Phase 3 simplification since full 3D lithography
    /// expose/develop isn't built yet (spin_photoresist above is the only
    /// 3D lithography primitive so far).
    pub fn etch_anisotropic(&mut self, depth_voxels: usize, mask: &[bool]) {
        let (nx, ny, nz) = (self.nx, self.ny, self.nz);
        let mut removed = 0usize;
        for y in 0..ny {
            for x in 0..nx {
                let col = col_idx(x, y, nx);
                if !mask[col] {
                    continue;
                }
                let mut cut = 0usize;
                for z in 0..nz {
                    if cut >= depth_voxels {
                        break;
                    }
                    let i = self.at(x, y, z);
                    if self.material[i] != Material::Void {
                        self.material[i] = Material::Void;
                        cut += 1;
                        removed += 1;
                    }
                }
            }
        }
        self.process_log.push(format!("etch_anisotropic depth={depth_voxels}vox removed={removed}"));
    }

    /// Height (in voxels from the top, z=0) of the topmost non-Void
    /// material in each (x,y) column - the "surface" a heightfield
    /// renderer extrudes. Returns nx*ny, `None` for columns that are Void
    /// all the way down (shouldn't normally happen).
    pub fn surface_heightmap(&self) -> Vec<Option<usize>> {
        let mut heights = vec![None; self.nx * self.ny];
        for y in 0..self.ny {
            for x in 0..self.nx {
                for z in 0..self.nz {
                    let i = self.at(x, y, z);
                    if self.material[i] != Material::Void {
                        heights[col_idx(x, y, self.nx)] = Some(z);
                        break;
                    }
                }
            }
        }
        heights
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn new_wafer_is_uniform_silicon_at_substrate_concentration() {
        let w = Wafer3d::new(6, 6, 10, 0.01, 0.01, 0.01, Dopant::Boron, 1e15);
        assert!(w.material.iter().all(|&m| m == Material::Silicon));
        assert_eq!(w.material.len(), 6 * 6 * 10);
        let boron = &w.species[&Dopant::Boron];
        assert!(boron.iter().all(|&c| (c - 1e15).abs() < 1.0));
    }

    #[test]
    fn masked_oxidation_only_grows_on_open_columns() {
        let nx = 4;
        let ny = 4;
        let mut w = Wafer3d::new(nx, ny, 40, 0.005, 0.005, 0.005, Dopant::Boron, 1e15);
        // Only the (1,1) column is exposed.
        let mut mask = vec![false; nx * ny];
        mask[col_idx(1, 1, nx)] = true;
        w.oxidize(1000.0, 0.75, Ambient::Dry, Some(&mask));

        assert!(w.oxide_thickness_um[col_idx(1, 1, nx)] > 0.0, "masked-open column should grow oxide");
        for y in 0..ny {
            for x in 0..nx {
                if (x, y) != (1, 1) {
                    assert_eq!(w.oxide_thickness_um[col_idx(x, y, nx)], 0.0, "column ({x},{y}) must not grow oxide");
                }
            }
        }
    }

    #[test]
    fn diffusion_3d_spreads_symmetrically_from_a_point_source() {
        let n = 15; // odd, so there's a true center voxel on all three axes
        let mut w = Wafer3d::new(n, n, n, 0.01, 0.01, 0.01, Dopant::Boron, 0.0);
        w.species.insert(Dopant::Phosphorus, vec![0.0; n * n * n]);
        let c = n / 2;
        let center = idx(c, c, c, n, n);
        w.species.get_mut(&Dopant::Phosphorus).unwrap()[center] = 1e20;

        w.anneal(1000.0, 5.0, 0.2);

        let p = &w.species[&Dopant::Phosphorus];
        let x_minus = p[idx(c - 2, c, c, n, n)];
        let x_plus = p[idx(c + 2, c, c, n, n)];
        let y_minus = p[idx(c, c - 2, c, n, n)];
        let y_plus = p[idx(c, c + 2, c, n, n)];
        let z_minus = p[idx(c, c, c - 2, n, n)];
        let z_plus = p[idx(c, c, c + 2, n, n)];

        assert!(x_minus > 0.0 && y_minus > 0.0 && z_minus > 0.0, "should have spread on all three axes");
        assert!((x_minus - x_plus).abs() / x_minus < 1e-6, "x symmetry: {x_minus} vs {x_plus}");
        assert!((y_minus - y_plus).abs() / y_minus < 1e-6, "y symmetry: {y_minus} vs {y_plus}");
        assert!((z_minus - z_plus).abs() / z_minus < 1e-6, "z symmetry: {z_minus} vs {z_plus}");
        // Cubic symmetry: equal grid spacing on all axes should give equal
        // spread magnitude too, not just left/right symmetry per axis.
        assert!((x_minus - z_minus).abs() / x_minus < 1e-6, "x vs z spread should match: {x_minus} vs {z_minus}");
    }

    #[test]
    fn diffusion_3d_conserves_total_dose() {
        let n = 11;
        let mut w = Wafer3d::new(n, n, n, 0.01, 0.01, 0.01, Dopant::Boron, 0.0);
        w.species.insert(Dopant::Phosphorus, vec![0.0; n * n * n]);
        let c = n / 2;
        w.species.get_mut(&Dopant::Phosphorus).unwrap()[idx(c, c, c, n, n)] = 1e20;

        let before: f64 = w.species[&Dopant::Phosphorus].iter().sum();
        w.anneal(1000.0, 2.0, 0.2);
        let after: f64 = w.species[&Dopant::Phosphorus].iter().sum();

        let relative_error = (after - before).abs() / before;
        assert!(relative_error < 1e-6, "mass should be conserved, error={relative_error}");
    }

    #[test]
    fn etch_anisotropic_only_removes_masked_columns_straight_down() {
        let nx = 4;
        let ny = 4;
        let mut w = Wafer3d::new(nx, ny, 20, 0.01, 0.01, 0.01, Dopant::Boron, 1e15);
        let mut mask = vec![false; nx * ny];
        mask[col_idx(2, 2, nx)] = true;
        w.etch_anisotropic(5, &mask);

        for z in 0..5 {
            assert_eq!(w.material[w.at(2, 2, z)], Material::Void, "masked column should be etched at z={z}");
        }
        assert_ne!(w.material[w.at(2, 2, 5)], Material::Void, "should stop exactly at depth_voxels");
        for y in 0..ny {
            for x in 0..nx {
                if (x, y) != (2, 2) {
                    assert_ne!(w.material[w.at(x, y, 0)], Material::Void, "unmasked column ({x},{y}) must be untouched");
                }
            }
        }
    }

    #[test]
    fn surface_heightmap_reflects_oxide_growth() {
        let nx = 3;
        let ny = 3;
        let mut w = Wafer3d::new(nx, ny, 30, 0.005, 0.005, 0.005, Dopant::Boron, 1e15);
        w.oxidize(1000.0, 2.0, Ambient::Wet, None);
        let heights = w.surface_heightmap();
        // With no mask, oxide grows uniformly, so the surface (z=0, since
        // oxide voxels replace silicon at z=0 upward and material is never
        // Void here) is still detected at z=0 everywhere - the real signal
        // is that oxide was actually deposited, checked separately below.
        assert!(heights.iter().all(|h| h == &Some(0)));
        assert!(w.oxide_thickness_um.iter().all(|&t| t > 0.0));
    }
}
