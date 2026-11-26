// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! End-to-end runner: the seven-stage recipe on the default tile, with a
//! report of what each stage did. This is the same evaluator the tests drive.

mod mk3;
use atto_units::{Kelvin, MTorr, Mm, Nm, PerCm2, Um, Watts};
use atto_geom::column::{coat_planarizing, etch, grow_film, strip_resist, ColumnGrid, EtchSpec, ProcessWindow};
use atto_geom::mask::{aerial_image, develop_map, rayleigh_nm, MaskPattern};
use atto_geom::svo::{Material, Wafer};
use atto_measure::{etch_depth_nm, film_thickness_nm, step_height_nm, top_run_widths_nm};
use atto_device::{
    cross_section::CrossSection,
    mos::{threshold_voltage, MosParams},
    poisson::{cv_curve, PoissonColumn},
};
use atto_dopant::{DopantField, DopantRecord, Type};
use atto_params::{Ambient, Species};
use atto_process::{
    anneal::{anneal_field, AnnealConfig, ThermalHistory, ThermalStep},
    etch::EtchModel,
    implant::{implant_into_field_masked, mask_from_top, Profile},
    oxidation::oxidise_geometry,
};

/// Mean junction depth over the window, for before/after comparison.
fn junction_summary(field: &DopantField, g: &ColumnGrid, surface_z: f64) -> Option<f64> {
    let m = field.junction_map(&g.cols, &g.active, g.x0, g.y0, g.pitch, g.nx, surface_z);
    let found: Vec<f64> = m.iter().filter_map(|v| *v).collect();
    if found.is_empty() {
        None
    } else {
        Some(found.iter().sum::<f64>() / found.len() as f64)
    }
}

fn main() {
    let t0 = std::time::Instant::now();
    let mut w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).expect("wafer");
    let win = ProcessWindow { center: [0.0, 0.0], size: [1200.0, 1200.0] };
    // The substrate is p-type; the field carries that background everywhere.
    let mut field = DopantField::new(atto_units::Nm::new(w.leaf_pitch_nm()), w.bounds.min);
    field.set_background(DopantRecord::new(Species::Boron, Type::Acceptor, 1.0e15));
    println!("AttoFab core — end-to-end recipe");
    println!("wafer   300 mm Si, tile 2048 nm at 16 nm, {} nodes", w.node_count());
    println!("window  {:.0} x {:.0} nm\n", win.size[0], win.size[1]);
    println!("{:<26} {:<38} {}", "stage", "result", "nodes");

    // 1. substrate prep: a 1.5 nm native oxide is recorded, not drawn
    println!("{:<26} {:<38} {}", "1 substrate prep", "native oxide 1.5 nm, Ra 0.2 nm", w.node_count());

    // 2. deposition: 100 nm of SiO2, directional
    let g0 = ColumnGrid::sample(&w, &win).expect("grid");
    let (g1, f) = grow_film(&g0, Material::SiliconDioxide, Nm::new(100.0), 0.0).expect("deposition");
    g1.write_back(&g0, &mut w);
    println!("{:<26} {:<38} {}", "2 deposition", format!("SiO2 {:.1} nm mean, nonunif {:.1}%", f.mean_nm, f.nonuniformity_pct), w.node_count());
    // atto_bridge exists to carry the pipeline's plain-data seam. The material code is the
    // stable identifier an external tool or a plugin sees, so reporting it is the one place the
    // CLI has a reason to call the bridge.
    println!("{:<26} {:<38} {}", "   material code", format!("{} (atto_bridge)", atto_bridge::material_code(Material::SiliconDioxide)), "");
    // S-07: what the same thickness would cost as an ALD. The point of ALD is that it is linear
    // in cycles and independent of dose above saturation, not that it is faster.
    if let Ok(n) = atto_process::ald::cycles_for(atto_units::Nm(f.mean_nm), &atto_process::ald::AldProcess::oxide()) {
        println!("{:<26} {:<38} {}", "   ALD equivalent", format!("{} cycles at {:.2} nm/cycle (assumed)", n, atto_process::ald::AldProcess::oxide().nm_per_cycle), "");
    }

    // 3. lithography: coat, expose, develop
    let g2 = ColumnGrid::sample(&w, &win).expect("grid");
    let (coated, c) = coat_planarizing(&g2, 144.0).expect("coat");
    let r = rayleigh_nm(0.4, atto_units::Nm::new(193.0), 0.85).expect("rayleigh");
    let mask = MaskPattern::LineSpace { line_nm: 150.0, pitch_nm: 320.0, count: 3, vertical: true, clear_field: true };
    let img = aerial_image(&g2, &mask, [0.0, 0.0], r / 2.355).expect("aerial image");
    let removed = develop_map(&img, 2.0, true);
    let (dev, _) = atto_geom::column::develop(&coated, &removed, 1.0).expect("develop");
    dev.write_back(&g2, &mut w);
    let cd = top_run_widths_nm(&dev, Material::Photoresist);
    println!("{:<26} {:<38} {}", "3 lithography", format!("R {:.1} nm, resist {:.0} nm, CD {:?}", r, c.mean_nm, cd.iter().map(|v| v.round()).collect::<Vec<_>>()), w.node_count());
    // L-10: shot-noise roughness at a nominal EUV dose. The square-root law is the point:
    // doubling the dose improves the roughness by 1.41, not by 2.
    match atto_process::stochastic::lwr_nm(atto_units::MjPerCm2(20.0), atto_units::Nm(16.0), atto_units::Nm(atto_process::stochastic::LWR_RESPONSE_NM)) {
        Ok(v) => println!("{:<26} {:<38} {}", "   shot-noise floor", format!("{:.1} nm at 20 mJ/cm2, pixel {:.0} nm (photon count only)", v, 16.0), ""),
        Err(e) => println!("{:<26} {:<38} {}", "   shot-noise floor", format!("not computed: {e}"), ""),
    }

    // 4. etch: plasma-driven rates, anisotropic
    let g3 = ColumnGrid::sample(&w, &win).expect("grid");
    let model = EtchModel::bundled_oxide();
    let rates = model.rates(Watts::new(200.0), MTorr::new(20.0));
    let spec = EtchSpec { rates, time_min: 1.8 };
    let (etched, st) = etch(&g3, &spec, 1.0).expect("etch");
    etched.write_back(&g3, &mut w);
    let depth = etch_depth_nm(&g3, &etched).unwrap_or(0.0);
    println!("{:<26} {:<38} {}", "4 etching", format!("200 W / 20 mT, depth {:.1} nm, mask lost {} cols", st.feature_depth_max_nm, st.mask_failed_columns), w.node_count());
    // E-10: the profile defects of the deepest feature. The width is taken as the grid pitch,
    // which is the narrowest feature this run can represent, so the aspect ratio reported is the
    // worst case the grid allows.
    if st.feature_depth_max_nm > 0.0 {
        let width = g3.pitch;
        match atto_process::profile::defects(atto_units::Nm(width), atto_units::Nm(st.feature_depth_max_nm), &atto_process::profile::DefectCoefficients::silicon()) {
            Ok(d) => println!("{:<26} {:<38} {}", "   profile", d.label(), ""),
            Err(e) => println!("{:<26} {:<38} {}", "   profile", format!("not computed: {e}"), ""),
        }
    }
    // E-12: what the same etch does to the mask. Resist is the mask here, so its rate is the
    // mask rate and the oxide rate is the film rate.
    if let Ok(e) = atto_process::mask_erosion::mask_erosion(atto_units::NmPerMin(rates.resist), atto_units::NmPerMin(rates.sio2), atto_units::Minutes(spec.time_min), atto_units::Nm(1000.0)) {
        println!("{:<26} {:<38} {}", "   mask erosion", e.label(), if e.survived() { "mask survived" } else { "mask consumed" });
    }

    // 5. doping: n+ implant written into the dopant field at leaf resolution
    let prof = Profile::from_species(Species::Phosphorus, 60.0).expect("profile");
    let dose = 5.0e14;
    // Re-sample after the etch: the implant depth reference is each column's own
    // silicon surface, so it must follow the geometry the etch produced.
    let g_dope = ColumnGrid::sample(&w, &win).expect("grid");
    // the resist still on the lines masks the beam, and lateral straggle lets
    // dopant creep under the edge (D-06)
    let open: Vec<bool> = mask_from_top(&g_dope).iter().map(|m| !m).collect();
    let masked = open.iter().filter(|o| !**o).count();
    let stats = implant_into_field_masked(&mut field, &g_dope, &prof, dose, Species::Phosphorus, Type::Donor, &open, 40.0);
    let analytic_xj = prof.junction_depth_nm(PerCm2::new(dose), 1.0e15).unwrap_or(0.0);
    println!("{:<26} {:<38} {}", "5 doping", format!("P 60 keV {:.0e}, {:?}, xj {:.0} nm, {} masked", dose, prof.shape, analytic_xj, masked), w.node_count());

    // 5b. anneal: diffuse the field and track the thermal budget
    let xj_before = junction_summary(&field, &g_dope, w.surface_z_nm);
    let mut history = ThermalHistory::new();
    history.push(ThermalStep { temp_c: 1000.0, minutes: 10.0 }).expect("thermal step");
    let cfg = AnnealConfig::new(Species::Phosphorus, Type::Donor, Kelvin::from_celsius(1000.0), 10.0);
    let ast = anneal_field(&mut field, &g_dope, &cfg).expect("anneal");
    println!("{:<26} {:<38} {}", "5b anneal", format!("1000 C 10 min, {} columns, Dt {:.2e} cm2", ast.columns, ast.dt_cm2), "");
    let xj_after = junction_summary(&field, &g_dope, w.surface_z_nm);
    if let (Some(b), Some(a)) = (xj_before, xj_after) {
        println!("{:<26} {:<38} {}", "junction shift", format!("mean {:.0} -> {:.0} nm ({:+.0} nm)", b, a, a - b), "");
    }
    println!("{:<26} {:<38} {}", "thermal budget", format!("{} step(s), {:.0} min, peak {:.0} C", history.steps.len(), history.total_minutes(), history.peak_c()), "");

    // 6. oxidation: dry oxide, blocked by the resist still on the lines
    let g_ox = ColumnGrid::sample(&w, &win).expect("grid");
    let (oxidised, ost) = oxidise_geometry(&g_ox, Ambient::Dry, 1000.0, 40.0, 1.0).expect("oxidation");
    oxidised.write_back(&g_ox, &mut w);
    println!("{:<26} {:<38} {}", "6 oxidation", format!("dry 1000 C 40 min, +{:.0} nm rise, {}/{} masked", ost.mean_rise_nm, ost.masked_columns, ost.columns), w.node_count());
    // D-15: how much of the dopant the grown oxide pulls out of the silicon. The
    // segregation coefficient is assumed (m = 0.3, boron-like: the oxide is a sink).
    if let Ok(sg) = atto_dopant::segregation::segregate(atto_units::PerCm2(dose), atto_units::Nm(200.0), atto_units::Nm(ost.mean_rise_nm.max(1.0)), atto_dopant::segregation::BORON_SEGREGATION_M) {
        println!("{:<26} {:<38} {}", "   segregation", sg.label(), "");
    }

    // 7. metrology
    let g4 = ColumnGrid::sample(&w, &win).expect("grid");
    let (_stripped, _) = strip_resist(&g4, 1.0).expect("strip");
    let film = film_thickness_nm(&g4, Material::SiliconDioxide).unwrap_or(0.0);
    let step = step_height_nm(&g4).unwrap_or(0.0);
    println!("{:<26} {:<38} {}", "7 metrology", format!("SiO2 {:.0} nm, step {:.0} nm, depth {:.0} nm", film, step, depth), w.node_count());

    // 8. device: threshold voltage from the cross-section
    let cs = CrossSection::build(&w, &field, &win, ost.interface).expect("cross-section");
    let mut vts: Vec<f64> = Vec::new();
    for ix in 0..cs.nx {
        if let Ok(r) = threshold_voltage(&cs, ix, &MosParams::default()) {
            vts.push(r.vt);
        }
    }
    if !vts.is_empty() {
        let mean = vts.iter().sum::<f64>() / vts.len() as f64;
        let lo = vts.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = vts.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("{:<26} {:<38} {}", "8 device", format!("Vt {:.3} V mean over {} columns ({:.3}..{:.3})", mean, vts.len(), lo, hi), "");
        if let Ok(r) = threshold_voltage(&cs, cs.nx / 2, &MosParams::default()) {
            println!("{:<26} {:<38} {}", "   centre column", r.label(), "");
            // avalanche breakdown of the same junction, from its bulk doping. The coefficients
            // are registered as assumed (aval.a, aval.b), and the model is a one-sided abrupt
            // junction, so the number is an order-of-magnitude statement, not a design value.
            let n = r.n_bulk.abs();
            if n > 0.0 {
                match atto_device::avalanche::breakdown_voltage_v(
                    atto_units::PerCm3(n),
                    &atto_device::avalanche::Ionisation::silicon(),
                ) {
                    Ok(v) => println!(
                        "{:<26} {:<38} {}",
                        "   breakdown",
                        format!("V_bd {:.0} V at N {:.1e} cm^-3 (abrupt, assumed)", v, n),
                        ""
                    ),
                    Err(e) => println!("{:<26} {:<38} {}", "   breakdown", format!("not computed: {e}"), ""),
                }
            }
        }
    } else {
        println!("{:<26} {:<38} {}", "8 device", "no column carries both oxide and doping", "");
    }

    // 9. C-V: solve Poisson across the gate sweep and differentiate
    if let Ok(pcol) = PoissonColumn::from_cross_section(&cs, cs.nx / 2, &MosParams::default()) {
        match cv_curve(&pcol, -1.5, 1.5, 7) {
            Ok(curve) => {
                print!("{:<26} ", "C-V (nF/cm2)");
                for q in &curve {
                    print!(" {:.0}@{:+.1}", q.c * 1e9, q.vg);
                }
                println!();
                let n_surf = pcol.doping.first().copied().unwrap_or(0.0);
                println!("{:<26} {:<38} {}", "   grid", format!("{} nodes at {:.0} nm, C_ox {:.0} nF/cm2", pcol.doping.len(), pcol.h_nm, pcol.cox() * 1e9), "");
                println!("{:<26} {:<38} {}", "   why flat", format!("surface N {:+.1e} is the blanket n+ layer; a p-surface would deplete", n_surf), "");
                // V-14: where this C-V rolls off, from the oxide capacitance and an assumed series
                // resistance. The corner and the time constant are the same number two ways.
                let area_cm2 = (win.size[0] * 1.0e-7) * (win.size[1] * 1.0e-7);
                let c_farad = pcol.cox() * area_cm2;
                if let Ok(fc) = atto_device::small_signal::corner_hz(atto_device::small_signal::DEFAULT_SERIES_R_OHM, c_farad) {
                    // The area is the simulated window, not a C-V pad — about 1.4e-8 cm^2, so the
                    // corner lands in the terahertz. Reported with the area in the label, because a
                    // bare "3e12 Hz" reads as an error rather than as a small capacitor.
                    println!("{:<26} {:<38} {}", "   C-f corner", format!("{:.1e} Hz at R {:.0} ohm, area {:.1e} cm2 (window, not a pad)", fc, atto_device::small_signal::DEFAULT_SERIES_R_OHM, area_cm2), "");
                }
            }
            Err(e) => println!("{:<26} {}", "C-V", e),
        }
    }

    // junction map from the field
    let jmap = field.junction_map(&g4.cols, &g4.active, g4.x0, g4.y0, g4.pitch, g4.nx, w.surface_z_nm);
    let found: Vec<f64> = jmap.iter().filter_map(|v| *v).collect();
    if !found.is_empty() {
        let mean = found.iter().sum::<f64>() / found.len() as f64;
        let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
        for v in &found {
            lo = lo.min(*v);
            hi = hi.max(*v);
        }
        println!("{:<26} {:<38} {}", "junction map", format!("{}/{} columns, {:.0}..{:.0} nm, mean {:.0}", found.len(), jmap.len(), lo, hi, mean), "");
    }
    println!("{:<26} {:<38} {}", "dopant field", format!("{} leaves, {:.0}% of dose retained", field.occupied(), stats.retained_fraction * 100.0), "");
    // net doping down the centre column
    let x = g4.x0 + (g4.nx / 2) as f64 * g4.pitch;
    let y = g4.y0 + (g4.ny / 2) as f64 * g4.pitch;
    print!("{:<26} ", "net doping (cm-3)");
    for d in [40.0, 80.0, 120.0, 160.0, 200.0, 300.0, 600.0] {
        let n = field.net_at(x, y, w.surface_z_nm - d);
        print!(" {:+.1e}@{}", n, d as i32);
    }
    println!();

    // Sheet resistance from the surface layer (M-07). The dose is the surface concentration over
    // the junction depth, activated at the anneal temperature, and the mobility is the
    // registered assumed constant — all three of those are stated in the label.
    if analytic_xj > 0.0 {
        let n_surf = field.net_at(x, y, w.surface_z_nm - 40.0).abs();
        if n_surf > 0.0 {
            let act = atto_dopant::activation::activate(
                atto_units::PerCm3(n_surf),
                atto_params::Species::Phosphorus,
                atto_units::Kelvin::from_celsius(1000.0),
            )
            .ok();
            let fraction = act.as_ref().map(|a| a.fraction).unwrap_or(1.0);
            let dose = atto_units::PerCm2::new(n_surf * fraction * analytic_xj * 1.0e-7);
            match atto_device::sheet::sheet_resistance(dose, atto_device::sheet::DEFAULT_MOBILITY_CM2_VS) {
                Ok(rs) => println!(
                    "{:<26} {:<38} {}",
                    "sheet resistance",
                    format!("{:.0} ohm/sq (active {:.0}% of surface, mu {:.0} assumed)", rs, fraction * 100.0, atto_device::sheet::DEFAULT_MOBILITY_CM2_VS),
                    ""
                ),
                Err(e) => println!("{:<26} {:<38} {}", "sheet resistance", format!("not computed: {e}"), ""),
            }
        }
    }

    let state = w.fingerprint() ^ field.fingerprint().rotate_left(17);
    println!("\nstate fingerprint {:016x}  (geometry ^ dopant)   elapsed {:.1} ms", state, t0.elapsed().as_secs_f64() * 1e3);
    println!("notes: the dopant field is at leaf resolution and carries species and type;");
    println!("       oxidation now consumes silicon and raises the surface, with a beak at mask edges.");
    mk3::report();
}

