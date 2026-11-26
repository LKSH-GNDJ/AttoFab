// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Design of experiments and Monte Carlo (build PRD W-01, W-02).
//!
//! Prints the normalised sensitivity of every output to every knob, then the
//! share of the threshold-voltage spread that each declared uncertainty owns.

use atto_study::{monte_carlo, run, sensitivity, Knob, Knobs, Perturbation};

fn main() {
    let base = Knobs::default();
    let Some(o) = run(&base) else {
        eprintln!("the base case does not run");
        std::process::exit(1);
    };
    println!("AttoFab study");
    println!("  base   xj {:.1} nm, tox {:.1} nm, Vt {:+.4} V\n", o.xj_nm, o.tox_nm, o.vt_v);

    println!("  relative sensitivity  (dy/y) / (dx/x)");
    println!("  {:<28} {:>8} {:>8} {:>8}", "knob", "xj", "tox", "Vt");
    for (k, s) in sensitivity(&base) {
        println!("  {:<28} {:>+8.3} {:>+8.3} {:>+8.3}", k.label(), s.xj_nm, s.tox_nm, s.vt_v);
    }

    let pert = [
        Perturbation { knob: Knob::AnnealTemp, rel_sigma: 0.02, label: "anneal temperature 2%" },
        Perturbation { knob: Knob::OxMinutes, rel_sigma: 0.05, label: "oxidation time 5%" },
        Perturbation { knob: Knob::Qf, rel_sigma: 0.20, label: "fixed charge 20% (assumed)" },
        Perturbation { knob: Knob::Bulk, rel_sigma: 0.01, label: "bulk doping 1%" },
    ];
    let Some(mc) = monte_carlo(&base, &pert, 2000, 12345) else {
        eprintln!("the Monte Carlo did not run");
        std::process::exit(1);
    };
    println!("\n  Monte Carlo, {} samples, seed 12345", mc.samples);
    // Vt sits near zero here, so a percentage of the mean would be meaningless
    println!("  Vt   mean {:+.4} V   sd {:.1} mV", mc.mean.vt_v, mc.sd.vt_v * 1e3);
    println!("  xj   mean {:.1} nm   sd {:.2} nm", mc.mean.xj_nm, mc.sd.xj_nm);
    println!("\n  share of the Vt variance");
    for c in &mc.contributions {
        println!("    {:<34} {:>6.1}%   corr {:+.2}", c.label, 100.0 * c.share, c.correlation);
    }
    println!("\n  note: the fixed charge carries a 20% uncertainty because it is an assumed");
    println!("        constant. That is the honest cost of not having measured it.");
}
