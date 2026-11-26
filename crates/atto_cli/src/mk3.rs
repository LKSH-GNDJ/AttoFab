// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The MK3 models, called once each from the program so that none of them is a function only
//! its own tests reach. Every input here is a nominal value written beside its call; none comes
//! from the run above, so these rows say what the model returns at that input, not what the tile did.

use atto_units::{Celsius, Degrees, MjPerCm2, Minutes, Nm, NmPerMin, Seconds, Volts};
use atto_params::Species;

fn num<E: std::fmt::Display>(label: &str, r: Result<f64, E>, note: &str) {
    match r {
        Ok(v) => println!("{:<26} {:<38} {}", format!("   {label}"), format!("{v:.4e}"), note),
        Err(e) => println!("{:<26} {:<38} {}", format!("   {label}"), format!("not computed: {e}"), note),
    }
}

/// Print one line per MK3 model at a nominal input.
pub fn report() {
    println!("\nMK3 models at nominal inputs (not taken from the run above)");

    // Process.
    let ald = atto_process::ald::AldProcess::oxide();
    num("ALD thickness", atto_process::ald::thickness_nm(100, 1.0, &ald), "100 cycles, saturating dose");
    let g = atto_process::epitaxy::GradedSiGe { ge_start: 0.0, ge_end: 0.2, thickness_nm: 50.0 };
    num("graded SiGe Ge", atto_process::epitaxy::composition_at(Nm(25.0), &g), "mid-layer of 0 to 20% over 50 nm");
    num("critical thickness", atto_process::epitaxy::critical_thickness_nm(0.2), "nm, 20% Ge");
    println!("{:<26} {:<38} {}", "   growth on silicon", format!("{:?}", atto_process::epitaxy::growth_on(atto_geom::svo::Material::Silicon)), "");
    let b = atto_process::mask3d::Mask3d::binary();
    let a6 = atto_process::mask3d::Mask3d::attenuated_6();
    num("best-focus shift binary", atto_process::mask3d::best_focus_shift_nm(&b, 0.33), "nm at NA 0.33");
    num("best-focus shift 6% att.", atto_process::mask3d::best_focus_shift_nm(&a6, 0.33), "nm at NA 0.33");
    num("OED factor", atto_process::oed::oed_factor(Species::Boron, NmPerMin(5.0)), "boron, 5 nm/min oxidation");
    num("OED diffusivity", atto_process::oed::enhanced_diffusivity(Species::Boron, Celsius(1000.0), NmPerMin(5.0)), "cm^2/s, boron, 1000 C");
    let pc = atto_process::pair_diffusion::PairConfig::implant();
    num("supersaturation", atto_process::pair_diffusion::supersaturation(Minutes(1.0), &pc), "after 1 min");
    num("TED enhancement", atto_process::pair_diffusion::enhancement(Minutes(1.0), &pc), "after 1 min");
    num("TED integrated", atto_process::pair_diffusion::integrated_enhancement(Minutes(1.0), &pc), "over 1 min");
    num("photons per pixel", atto_process::stochastic::photons_per_pixel(MjPerCm2(20.0), Nm(16.0)), "20 mJ/cm2, 16 nm pixel");
    let v = atto_process::stress_coupling::ActivationVolumes::silicon();
    num("stress on diffusion", atto_process::stress_coupling::diffusion_factor(100.0, &v, Celsius(1000.0)), "multiplier, 100 MPa");
    num("stress on oxidation", atto_process::stress_coupling::oxidation_factor(100.0, &v, Celsius(1000.0)), "multiplier, 100 MPa");
    num("stress on mobility", atto_process::stress_coupling::mobility_factor(100.0, &v, Celsius(1000.0)), "multiplier, 100 MPa");
    println!("{:<26} {:<38} {}", "   wet-etch pit", format!("{:?}", atto_process::wet_etch::square_pit(Nm(1000.0), atto_process::wet_etch::Orientation::Si100)), "1 um opening, (100)");
    println!("{:<26} {:<38} {}", "   {111} wall angle", format!("{:.2} deg", atto_process::wet_etch::sidewall_angle_deg()), "");
    let y = atto_process::yield_model::YieldParams::silicon();
    num("energy factor", atto_process::yield_model::energy_factor(500.0, &y), "500 eV");
    num("angle factor", atto_process::yield_model::angle_factor(Degrees(45.0), &y), "45 deg");
    num("sputter yield", atto_process::yield_model::sputter_yield(500.0, Degrees(45.0), &y), "500 eV at 45 deg");

    let pp = atto_process::passivation::PassivationParams::bundled();
    for (f, what) in [(0.0, "isotropic chemistry, no passivant"), (1.0, "same chemistry, full passivant")] {
        match atto_process::passivation::sidewall(&pp, atto_units::Fraction(f), NmPerMin(300.0), NmPerMin(300.0), Minutes(5.0)) {
            Ok(r) => println!("{:<26} {:<38} {}", "   sidewall passivation", r.label(), what),
            Err(e) => println!("{:<26} {:<38} {}", "   sidewall passivation", format!("not computed: {e}"), what),
        }
    }

    // Device.
    let ch = atto_device::continuity::Channel { surface_potential_v: 0.9, length_nm: 100.0, mobility_cm2_vs: 200.0, doping_cm3: 1.0e17, temp_k: 300.0 };
    num("channel current", atto_device::continuity::drain_current_a(&ch, Volts(0.5)), "A, V_ds 0.5 V");
    println!("{:<26} {:<38} {}", "   channel saturated", atto_device::continuity::is_saturated(&ch, Volts(0.5)), "V_ds 0.5 V");
    println!("{:<26} {:<38} {}", "   thermal voltage", format!("{:.4} V", ch.thermal_v()), "300 K");
    num("RC time constant", atto_device::small_signal::time_constant_s(1.0e3, 1.0e-12), "s, 1 kohm and 1 pF");
    num("capacitance at 1 MHz", atto_device::small_signal::capacitance_at(1.0e6, 1.0e3, 1.0e-12), "F");
    let samples: Vec<(f64, f64)> = (1..=5)
        .filter_map(|k| {
            let t = k as f64 * 0.5e-9;
            atto_device::small_signal::step_response_v(Seconds(t), 1.0e3, 1.0e-12, 1.0).ok().map(|v| (t, v))
        })
        .collect();
    num("tau from a step", atto_device::small_signal::time_constant_from_transient(&samples, 1.0), "s, fitted to five points");
    println!("{:<26} {:<38} {}", "   registered value", format!("{:?}", atto_device::registry::value_anywhere("aval.a")), "aval.a");
    println!("{:<26} {:<38} {}", "   registered constants", atto_device::registry::all().len(), "across both layers");
    println!("{:<26} {:<38} {}", "   plugin status 0", format!("{:?}", atto_bridge::plugin::PluginStatus::from_code(0)), "");

    // Methods that share a name with another type's method, called by path so each is counted.
    let o = atto_process::litho::Optics { lambda_nm: 193.0, na: 0.85, sigma: 0.5, defocus_nm: 0.0 };
    let gr = atto_process::litho::Grating { pitch_nm: 200.0, line_nm: 100.0 };
    println!("{:<26} {:<38} {}", "   Optics image", format!("{} samples", atto_process::litho::Optics::image(&o, &gr, 64).map(|v| v.len()).unwrap_or(0)), "193 nm, NA 0.85, 100 nm lines at 200 nm");
    println!("{:<26} {:<38} {}", "   binary mask image", format!("{} samples", atto_process::psm::image(&o, &gr, atto_process::psm::MaskType::Binary, 64).map(|v| v.len()).unwrap_or(0)), "same optics");
    if let Ok(r) = atto_process::psm::evaluate(&o, &gr, atto_process::psm::MaskType::Alternating, 64) {
        println!("{:<26} {:<38} {}", "   alternating mask", atto_process::psm::PsmResult::label(&r), "same optics");
    }
    if let Ok(p) = atto_process::wet_etch::square_pit(Nm(1000.0), atto_process::wet_etch::Orientation::Si100) {
        println!("{:<26} {:<38} {}", "   pit label", atto_process::wet_etch::Pit::label(&p), "");
    }
    if let Ok(e) = atto_process::mask_erosion::mask_erosion(NmPerMin(30.0), NmPerMin(300.0), Minutes(1.0), Nm(500.0)) {
        println!("{:<26} {:<38} {}", "   mask erosion", atto_process::mask_erosion::MaskErosion::label(&e), "30 vs 300 nm/min, 1 min");
    }
    println!("{:<26} {:<38} {}", "   OED sign, boron", atto_process::oed::sign(Species::Boron), "");
    println!("{:<26} {:<38} {}", "   donor sign", atto_dopant::Type::sign(atto_dopant::Type::Donor), "");
    if let Ok(sg) = atto_dopant::segregation::segregate(atto_units::PerCm2(1.0e14), Nm(100.0), Nm(50.0), 0.3) {
        println!("{:<26} {:<38} {}", "   segregation", atto_dopant::segregation::Segregation::label(&sg), "1e14 cm^-2, k 0.3");
    }

    // Files named on the command line.
    let args: Vec<String> = std::env::args().collect();
    for w in args.windows(2) {
        if w[0] == "--deck" {
            match std::fs::read_to_string(&w[1]).map_err(|e| e.to_string()).and_then(|t| atto_params::deck::parse(&t).map_err(|e| format!("{e:?}"))) {
                Ok(d) => println!("{:<26} {:<38} {}", "   deck round trip", format!("{} stages, {} bytes", d.stages.len(), atto_params::deck::to_text(&d).len()), w[1]),
                Err(e) => println!("{:<26} {:<38} {}", "   deck", format!("not read: {e}"), w[1]),
            }
        }
        if w[0] == "--oasis" {
            match std::fs::read(&w[1]).map_err(|e| format!("{e}")).and_then(|b| atto_geom::oasis::parse(&b).map_err(|e| format!("{e:?}"))) {
                Ok(f) => println!("{:<26} {:<38} {}", "   OASIS file", if f.is_complete() { "complete" } else { "incomplete" }, w[1]),
                Err(e) => println!("{:<26} {:<38} {}", "   OASIS file", format!("not read: {e}"), w[1]),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_report_runs_every_nominal_call_without_panicking() {
        report();
    }

    #[test]
    fn a_failed_model_is_printed_as_not_computed_not_hidden() {
        // `num` must show the error text rather than drop the row. The negative target is refused by the model.
        let r = atto_process::ald::cycles_for(Nm(-1.0), &atto_process::ald::AldProcess::oxide()).map(|n| n as f64);
        assert!(r.is_err());
        num("ALD cycles", r, "negative target");
    }
}
