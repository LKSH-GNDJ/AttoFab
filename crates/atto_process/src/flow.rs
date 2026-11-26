// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! One layout, one recipe, one cross-section (build PRD L-02, D-01, X-01).
//!
//! This is the only place that runs the whole chain in order: coat, expose from
//! a GDS layer, develop to resist, etch, implant through what the resist left
//! open, anneal, oxidise, and read the junction back out. Every step is a
//! library call that has its own tests; what is tested here is that they compose.

use atto_units::{Kelvin, MTorr, Mm, Nm, Um, Watts};
use crate::etch::EtchModel;
use crate::implant::{implant_into_field_masked, mask_from_top, Profile};
use crate::litho::{apply_resist, expose_layout, Optics};
use crate::oxidation::oxidation_nm;
use crate::anneal::{diffuse, AnnealConfig};
use atto_dopant::{DopantField, Type as DopantType};
use atto_geom::column::{etch, set_interval, solid_top, ColumnGrid, EtchSpec, ProcessWindow};
use atto_geom::gds::Layout;
use atto_geom::svo::{Material, Wafer};
use atto_params::{Ambient, Species as DopantSpecies};

/// Everything a layout needs to become a device.
#[derive(Debug, Clone)]
pub struct FlowSpec {
    pub layer: i16,
    pub optics: Optics,
    pub field_nm: f64,
    pub threshold: f64,
    pub oxide_nm: f64,
    pub resist_nm: f64,
    pub etch_power_w: f64,
    pub etch_pressure_mtorr: f64,
    pub etch_min: f64,
    /// Fraction of etched material that lands back on the walls (PRD E-07).
    pub redeposition_sticking: f64,
    /// Wall angle above which a deposit counts as a facet, degrees.
    pub redeposition_facet_deg: f64,
    pub implant_energy_kev: f64,
    pub implant_dose_cm2: f64,
    pub anneal_temp_c: f64,
    pub anneal_min: f64,
    pub oxidation_temp_c: f64,
    pub oxidation_min: f64,
    /// Bulk doping the junction is measured against, cm^-3.
    pub background_cm3: f64,
}

impl Default for FlowSpec {
    fn default() -> Self {
        FlowSpec {
            layer: 1,
            optics: Optics { lambda_nm: 193.0, na: 0.93, sigma: 0.6, defocus_nm: 0.0 },
            field_nm: 3000.0,
            threshold: 0.35,
            oxide_nm: 120.0,
            resist_nm: 80.0,
            etch_power_w: 200.0,
            etch_pressure_mtorr: 20.0,
            redeposition_sticking: 0.15,
            redeposition_facet_deg: 8.0,
            etch_min: 1.0,
            implant_energy_kev: 60.0,
            implant_dose_cm2: 5.0e14,
            anneal_temp_c: 1000.0,
            anneal_min: 10.0,
            oxidation_temp_c: 1000.0,
            oxidation_min: 30.0,
            background_cm3: 1.0e15,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
/// What one layout through the whole recipe produced.
pub struct FlowResult {
    /// Fraction of the field left covered by resist.
    pub resist_fraction: f64,
    /// Deepest etch, nm.
    pub etch_depth_nm: f64,
    /// What the etch put back on the walls (PRD E-07). `None` when nothing was etched.
    pub redeposition: Option<crate::redeposition::Redeposition>,
    /// Columns the resist protected from the implant.
    pub masked_columns: usize,
    /// Columns the beam reached.
    pub open_columns: usize,
    /// Junction depth below the silicon surface in a representative open column.
    pub junction_nm: Option<f64>,
    /// Why the junction is absent, when it is. An absent number that does not say
    /// why reads as "there is no junction" rather than "the window was too short".
    pub junction_reason: String,
    /// Oxide grown at the end, nm.
    pub tox_nm: f64,
}

/// What the etch left on the walls, from the depths it actually removed (PRD E-07).
///
/// The depth profile is the deepest removal per open column, so the model sees the etch
/// the grid reports rather than a nominal rate. `None` when nothing was removed, which is
/// the honest answer: there is no redeposition without an etch.
fn redeposition_of(
    stats: &atto_geom::column::EtchStats,
    spec: &FlowSpec,
) -> Option<crate::redeposition::Redeposition> {
    if !(stats.feature_depth_max_nm > 0.0) {
        return None;
    }
    let params = crate::redeposition::RedepositionParams::new(
        atto_units::Fraction::new(spec.redeposition_sticking),
        spec.redeposition_facet_deg,
    )?;
    let depth = stats.feature_depth_max_nm;
    crate::redeposition::redistribute(&[depth], atto_units::Nm::new(spec.field_nm), &params)
}

/// Run the whole recipe for one layout layer.
pub fn run_flow(layout: &Layout, spec: &FlowSpec) -> Result<FlowResult, String> {
    if !(spec.field_nm > 0.0) || !(spec.oxide_nm > 0.0) || !(spec.resist_nm > 0.0) {
        return Err("field, oxide and resist must be positive".into());
    }
    let mut w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).map_err(|e| format!("{e:?}"))?;
    let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };

    // 1. oxide, then resist from the layout
    let g0 = ColumnGrid::sample(&w, &win).map_err(|e| format!("{e:?}"))?;
    let mut coated = g0.clone();
    for k in 0..coated.cols.len() {
        let top = solid_top(&coated.cols[k]).map(|(t, _)| t).unwrap_or(coated.z_bottom);
        set_interval(&mut coated.cols[k], top, top + spec.oxide_nm, Material::SiliconDioxide);
    }
    let pattern = expose_layout(layout, spec.layer, &spec.optics, spec.field_nm, 500.0, 256, spec.threshold, true)
        .map_err(|e| e.to_string())?;
    apply_resist(&mut coated, &pattern, spec.field_nm, spec.resist_nm);
    let resist_fraction = pattern.iter().filter(|b| **b).count() as f64 / pattern.len().max(1) as f64;

    // 2. etch
    let model = EtchModel::bundled_oxide();
    let rates = model.rates(Watts::new(spec.etch_power_w), MTorr::new(spec.etch_pressure_mtorr));
    let es = EtchSpec { rates, time_min: spec.etch_min };
    let (etched, stats) = etch(&coated, &es, 1.0).map_err(|e| format!("{e:?}"))?;
    etched.write_back(&coated, &mut w);

    // 3. implant through what the resist left open
    let g1 = ColumnGrid::sample(&w, &win).map_err(|e| format!("{e:?}"))?;
    // The resist masks the implant as well as the etch, and the etch consumes
    // it, so the pattern has to be read before the etch runs.
    let open: Vec<bool> = mask_from_top(&coated).iter().map(|m| !m).collect();
    let open_columns = open.iter().filter(|o| **o).count();
    let masked_columns = open.len() - open_columns;
    let mut field = DopantField::new(atto_units::Nm::new(16.0), (0.0, 0.0, 0.0));
    field.set_background(atto_dopant::DopantRecord::new(DopantSpecies::Boron, DopantType::Acceptor, spec.background_cm3));
    let profile = Profile::from_species(DopantSpecies::Phosphorus, spec.implant_energy_kev)
        .map_err(|e| format!("{e:?}"))?;
    implant_into_field_masked(
        &mut field,
        &g1,
        &profile,
        spec.implant_dose_cm2,
        DopantSpecies::Phosphorus,
        DopantType::Donor,
        &open,
        40.0,
    );

    // 4. anneal a representative open column and read the junction
    let mut best: Option<f64> = None;
    let mut why = String::from("no open column to measure");
    // the first column the beam actually reached, not the middle of the field:
    // the middle may sit under resist, and then there is nothing to measure
    let picked = (0..g1.cols.len()).find(|k| open.get(*k).copied().unwrap_or(false));
    let (ix, iy) = match picked {
        Some(k) => (k % g1.nx, k / g1.nx),
        None => (0, 0),
    };
    // The junction is measured from the *silicon* surface, not from the top of
    // the stack: measuring from the oxide top displaces the whole profile by the
    // oxide thickness and puts the peak 120 nm too deep.
    let surface_z = g1.cols[g1.idx(ix, iy)]
        .iter()
        .filter(|r| matches!(r.material, Material::Silicon))
        .map(|r| r.z1)
        .fold(f64::MIN, f64::max);
    let surface_z = if surface_z.is_finite() { surface_z } else { g1.z_bottom };
    let x = g1.x0 + (ix as f64 + 0.5) * g1.pitch;
    let y = g1.y0 + (iy as f64 + 0.5) * g1.pitch;
    if picked.is_some() {
        let leaf = field.leaf_nm();
        // the window has to be deep enough for a long anneal: a 30 min drive
        // pushes the junction past a micron, and a short window reports None
        let depth = 2000.0;
        let n = (depth / leaf) as usize;
        // Keep the sign: the junction is where the net doping changes from
        // donor to acceptor, and taking the magnitude destroys exactly the
        // information that defines it.
        let mut c: Vec<f64> = (0..n)
            .map(|i| field.net_at(x, y, surface_z - (i as f64 + 0.5) * leaf))
            .collect();
        let cfg = AnnealConfig::new(DopantSpecies::Phosphorus, DopantType::Donor, Kelvin::from_celsius(spec.anneal_temp_c), spec.anneal_min);
        diffuse(&mut c, atto_units::Nm::new(leaf), &cfg).map_err(|e| format!("{e:?}"))?;
        why = format!("no donor-to-acceptor crossing within the {depth:.0} nm window");
        for i in 1..n {
            if c[i - 1] > 0.0 && c[i] <= 0.0 {
                let f = c[i - 1] / (c[i - 1] - c[i]);
                best = Some((i as f64 - 0.5 + f) * leaf);
                why = String::from("junction found");
                break;
            }
        }
    }

    // 5. oxidation
    //
    // The thin regime is included: below about 30 nm Deal-Grove alone under-predicts dry oxide
    // growth, and a flow that starts from bare silicon spends most of its time in exactly that
    // regime. The constants behind it are registered as assumed (ox.massoud.c / .l).
    let tox = crate::massoud::oxidation_nm_thin(
        Ambient::Dry,
        atto_units::Kelvin::from_celsius(spec.oxidation_temp_c),
        spec.oxidation_min,
        atto_units::Nm(0.0),
        &crate::massoud::Massoud::bundled(),
    )
    .map_err(|e| format!("{e}"))?
    .0;

    Ok(FlowResult {
        resist_fraction,
        etch_depth_nm: stats.feature_depth_max_nm,
        redeposition: redeposition_of(&stats, spec),
        masked_columns,
        open_columns,
        junction_nm: best,
        junction_reason: why,
        tox_nm: tox,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use atto_geom::gds::{Polygon, Structure};

    fn lines() -> Layout {
        let mut polygons = Vec::new();
        for i in 0..6 {
            let x = i as f64 * 0.5;
            polygons.push(Polygon {
                layer: 1,
                datatype: 0,
                points: vec![(x, 0.0), (x + 0.25, 0.0), (x + 0.25, 1.0), (x, 1.0), (x, 0.0)],
            });
        }
        Layout { units_per_micron: 1000.0, structures: vec![Structure { name: "TOP".into(), polygons }] }
    }

    #[test]
    fn the_flow_grows_more_oxide_than_deal_grove_alone() {
        // the thin regime was wired into step 5 of the flow; without this test the change is
        // a number nobody pins. The same conditions through pure Deal-Grove must come out lower.
        let spec = FlowSpec::default();
        let r = run_flow(&lines(), &spec).expect("the flow runs");
        let dg = oxidation_nm(
            Ambient::Dry,
            atto_units::Kelvin::from_celsius(spec.oxidation_temp_c),
            spec.oxidation_min,
            0.0,
        )
        .expect("deal-grove");
        assert!(r.tox_nm > dg, "the flow includes the thin regime: {} vs {}", r.tox_nm, dg);
        assert!((r.tox_nm - dg) / dg > 1e-3, "and it is a measurable difference: {} vs {}", r.tox_nm, dg);
    }

    #[test]
    fn a_layout_runs_the_whole_recipe() {
        let r = run_flow(&lines(), &FlowSpec::default()).expect("the flow runs");
        assert!(r.resist_fraction > 0.35 && r.resist_fraction < 0.65, "resist {:.2}", r.resist_fraction);
        assert!(r.etch_depth_nm > 0.0, "etch depth {:.1}", r.etch_depth_nm);
        assert!(r.open_columns > 0 && r.masked_columns > 0, "the mask must split the field");
        assert!(r.tox_nm > 35.0 && r.tox_nm < 130.0, "tox {:.1}", r.tox_nm);
        assert!(r.junction_nm.is_some(), "an open column must have a junction: {}", r.junction_reason);
        let xj = r.junction_nm.unwrap();
        // A band is not a check: an 800 nm window swallowed a 1.5x error, so the
        // window is narrowed to what this recipe can produce.
        assert!(xj > 290.0 && xj < 355.0, "junction {xj:.0} nm");
    }

    #[test]
    fn the_flow_reports_what_the_etch_put_back() {
        // E-07 wired, not merely implemented: the recipe's own etch depths drive the model
        let r = run_flow(&lines(), &FlowSpec::default()).expect("the flow runs");
        let rd = r.redeposition.expect("an etch deep enough to matter must report one");
        assert_eq!(rd.sidewall_nm.len(), 1, "one segment, one wall thickness");
        assert!(rd.sidewall_nm[0] > 0.0, "material came back: {}", rd.sidewall_nm[0]);
        assert!(rd.returned_fraction <= 0.15 + 1e-12, "never more than the sticking allows");
    }

    #[test]
    fn no_etch_means_no_redeposition() {
        // the honest absence: nothing removed, nothing to put back
        let mut spec = FlowSpec::default();
        spec.etch_min = 0.0;
        let r = run_flow(&lines(), &spec).expect("the flow still runs");
        // a zero-length etch removes nothing, so the model has nothing to put back and
        // must say so rather than reporting an empty deposit as a result
        assert_eq!(r.etch_depth_nm, 0.0, "no time, no depth");
        assert!(r.redeposition.is_none(), "no etch means no redeposition");
    }

    #[test]
    fn the_flow_is_repeatable() {
        let a = run_flow(&lines(), &FlowSpec::default()).unwrap();
        let b = run_flow(&lines(), &FlowSpec::default()).unwrap();
        assert_eq!(a, b, "the same layout and recipe must give the same result");
    }

    #[test]
    fn a_longer_anneal_deepens_the_junction() {
        let mut a = FlowSpec::default();
        a.anneal_min = 2.0;
        let mut b = FlowSpec::default();
        b.anneal_min = 10.0;
        let ra = run_flow(&lines(), &a).unwrap();
        let rb = run_flow(&lines(), &b).unwrap();
        assert_eq!(ra.etch_depth_nm, rb.etch_depth_nm, "the anneal does not change the etch");
        let (xa, xb) = (ra.junction_nm.expect("2 min"), rb.junction_nm.expect("30 min"));
        assert!(xb > xa, "10 min must beat 2 min: {xa:.0} -> {xb:.0} nm");
    }

    #[test]
    fn a_layout_that_masks_everything_implants_nothing() {
        // An empty layer is an opaque field: nothing prints clear, so resist
        // covers everything and the beam never reaches the silicon. A polygon
        // would be the opposite — clear — and would remove the resist.
        let dark = Layout { units_per_micron: 1000.0, structures: vec![] };
        let r = run_flow(&dark, &FlowSpec::default()).unwrap();
        assert!(r.resist_fraction > 0.9, "resist {:.2}", r.resist_fraction);
        assert_eq!(r.open_columns, 0, "nothing is open");
    }

    #[test]
    fn the_default_recipe_finishes_inside_its_budget() {
        // X-04. This is a tripwire, not a benchmark: the bound is loose enough
        // that machine speed does not decide it, and tight enough that an
        // accidental O(n^2) inside a loop does.
        use std::time::Instant;
        let spec = FlowSpec::default();
        let layout = lines();
        let t0 = Instant::now();
        let r = run_flow(&layout, &spec).expect("the default recipe runs");
        let elapsed = t0.elapsed().as_secs_f64();
        assert!(r.open_columns > 0 || r.masked_columns > 0, "the recipe did something");
        assert!(
            elapsed < 30.0,
            "the default recipe took {elapsed:.1} s, beyond its 30 s budget"
        );
    }

    #[test]
    fn bad_input_is_rejected() {
        let mut spec = FlowSpec::default();
        spec.field_nm = 0.0;
        assert!(run_flow(&lines(), &spec).is_err());
    }
}
