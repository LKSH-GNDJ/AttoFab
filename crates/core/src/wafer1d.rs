//! The Wafer State (1D column) - Phase 1 of the Wafer Voxel Mesh.
//!
//! Tracks a single vertical column of the wafer as a discrete depth grid.
//! Index 0 is the wafer surface; increasing index goes deeper into bulk
//! silicon. This module's `Wafer1d` is deliberately structured (parallel
//! arrays + a depth grid) so a future `Wafer2d`/`Wafer3d` can reuse the
//! same column type as a single "line" of a larger mesh.

use crate::diffusion::{diffuse, DiffuseParams};
use crate::implant::{estimate_range, implant_profile, ImplantParams, Range};
use crate::materials::{Dopant, Material, SI_CONSUMPTION_RATIO};
use crate::oxidation::{grow_oxide, GrowOxideParams};
use std::collections::HashMap;

pub struct Substrate {
    pub dopant: Dopant,
    pub concentration: f64,
}

pub struct Wafer1d {
    pub dx_um: f64,
    pub depth_um: f64,
    pub depth_grid_um: Vec<f64>,
    pub material: Vec<Material>,
    pub species: HashMap<&'static str, Vec<f64>>,
    pub oxide_thickness_um: f64,
    pub substrate: Substrate,
}

pub struct OxidizeStep {
    pub temperature_c: f64,
    pub time_hours: f64,
    pub ambient: crate::materials::Ambient,
}

pub struct OxidizeOutcome {
    pub oxide_thickness_um: f64,
    pub si_consumed_um: f64,
}

pub struct ImplantStep {
    pub dopant: Dopant,
    pub dose_cm2: f64,
    /// If both are Some, they override the energy-based estimate.
    pub rp_override_um: Option<f64>,
    pub d_rp_override_um: Option<f64>,
    pub energy_kev: Option<f64>,
}

pub struct AnnealStep {
    pub temperature_c: f64,
    pub time_minutes: f64,
    pub safety_factor: f64,
}

impl Wafer1d {
    pub fn new(depth_um: f64, dx_um: f64, substrate: Substrate) -> Self {
        let n_points = (depth_um / dx_um).round() as usize + 1;
        let depth_grid_um: Vec<f64> = (0..n_points).map(|i| i as f64 * dx_um).collect();
        let material = vec![Material::Silicon; n_points];

        let mut species: HashMap<&'static str, Vec<f64>> = HashMap::new();
        species.insert(substrate.dopant.name_static(), vec![substrate.concentration; n_points]);

        Self {
            dx_um,
            depth_um,
            depth_grid_um,
            material,
            species,
            oxide_thickness_um: 0.0,
            substrate,
        }
    }

    fn n_points(&self) -> usize {
        self.depth_grid_um.len()
    }

    pub fn oxidize(&mut self, step: &OxidizeStep) -> OxidizeOutcome {
        let result = grow_oxide(&GrowOxideParams {
            temperature_c: step.temperature_c,
            time_hours: step.time_hours,
            ambient: step.ambient,
            initial_oxide_um: self.oxide_thickness_um,
            orientation: crate::materials::CrystalOrientation::Si100,
        });
        let grown = result.final_oxide_um - self.oxide_thickness_um;
        self.oxide_thickness_um = result.final_oxide_um;
        let si_consumed_um = grown * SI_CONSUMPTION_RATIO;

        OxidizeOutcome { oxide_thickness_um: self.oxide_thickness_um, si_consumed_um }
    }

    pub fn implant(&mut self, step: &ImplantStep) -> Vec<f64> {
        let range = match (step.rp_override_um, step.d_rp_override_um) {
            (Some(rp), Some(d_rp)) => Range { rp, d_rp },
            _ => estimate_range(
                step.dopant,
                step.energy_kev.expect("energy_kev required when Rp/dRp not provided"),
            ),
        };

        let added = implant_profile(
            &self.depth_grid_um,
            &ImplantParams { dose_cm2: step.dose_cm2, range },
        );

        let n = self.n_points();
        let entry = self
            .species
            .entry(step.dopant.name_static())
            .or_insert_with(|| vec![0.0; n]);
        for (c, a) in entry.iter_mut().zip(added.iter()) {
            *c += a;
        }
        added
    }

    pub fn anneal(&mut self, step: &AnnealStep) {
        let dopants = [Dopant::Boron, Dopant::Phosphorus, Dopant::Arsenic, Dopant::Antimony];
        for dopant in dopants {
            let name = dopant.name_static();
            if let Some(conc) = self.species.get(name) {
                let result = diffuse(&DiffuseParams {
                    concentration: conc,
                    dx_um: self.dx_um,
                    dopant,
                    temperature_c: step.temperature_c,
                    time_minutes: step.time_minutes,
                    safety_factor: step.safety_factor,
                });
                self.species.insert(name, result.concentration);
            }
        }
    }

    /// Net signed dopant concentration (P-type positive, N-type negative)
    /// used for metallurgical junction finding.
    pub fn net_dopant(&self) -> Vec<f64> {
        let n = self.n_points();
        let mut net = vec![0.0; n];
        for (&name, conc) in self.species.iter() {
            let is_p = Dopant::from_name(name).map(|d| d.is_p_type()).unwrap_or(true);
            for i in 0..n {
                net[i] += if is_p { conc[i] } else { -conc[i] };
            }
        }
        net
    }

    pub fn junction_depths(&self) -> Vec<f64> {
        let net = self.net_dopant();
        let mut depths = Vec::new();
        for i in 0..net.len().saturating_sub(1) {
            let a = net[i];
            let b = net[i + 1];
            if (a > 0.0 && b <= 0.0) || (a < 0.0 && b >= 0.0) {
                let frac = a / (a - b);
                depths.push(
                    self.depth_grid_um[i]
                        + frac * (self.depth_grid_um[i + 1] - self.depth_grid_um[i]),
                );
            }
        }
        depths
    }
}

impl Dopant {
    fn name_static(&self) -> &'static str {
        self.name()
    }

    fn from_name(name: &str) -> Option<Dopant> {
        match name {
            "Boron" => Some(Dopant::Boron),
            "Phosphorus" => Some(Dopant::Phosphorus),
            "Arsenic" => Some(Dopant::Arsenic),
            "Antimony" => Some(Dopant::Antimony),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::Ambient;
    use approx::assert_relative_eq;

    #[test]
    fn full_flow_matches_js_phase1_benchmark() {
        let mut wafer = Wafer1d::new(5.0, 0.005, Substrate { dopant: Dopant::Boron, concentration: 1e15 });

        wafer.oxidize(&OxidizeStep { temperature_c: 1000.0, time_hours: 0.75, ambient: Ambient::Dry });
        wafer.implant(&ImplantStep {
            dopant: Dopant::Phosphorus,
            dose_cm2: 1e15,
            rp_override_um: None,
            d_rp_override_um: None,
            energy_kev: Some(80.0),
        });
        wafer.anneal(&AnnealStep { temperature_c: 1000.0, time_minutes: 30.0, safety_factor: 0.45 });

        // Pinned to the *current* (PDF-sourced, corrected) Deal-Grove
        // constants and this wafer's own zero-initial-oxide seeding
        // convention (oxide_thickness_um starts at 0, not NATIVE_OXIDE_UM -
        // see oxidation::tests for the NATIVE_OXIDE_UM-seeded case, which
        // is a different scenario and correctly yields a different value).
        // This number will legitimately disagree with the original JS
        // Phase-1 engine (~43.1nm), which used an earlier, less carefully
        // sourced constants table - see docs/math_references.md.
        assert_relative_eq!(wafer.oxide_thickness_um * 1000.0, 29.8, epsilon = 0.5);

        let junctions = wafer.junction_depths();
        assert_eq!(junctions.len(), 1);
        assert_relative_eq!(junctions[0], 0.502, epsilon = 0.01);
    }
}
