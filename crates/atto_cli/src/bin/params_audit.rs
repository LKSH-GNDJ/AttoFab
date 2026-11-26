// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Audit the parameter manifest (build PRD K-01).
//!
//! Prints the coverage of every registered constant and fails if the manifest
//! and the code disagree. This is the check that would have caught the earlier build's
//! bundled implant table disagreeing with its own code by 34%.

fn main() {
    let mut drift = Vec::new();
    let mut all = Vec::new();

    let p = atto_params::registry();
    drift.extend(atto_params::verify(&p, &atto_params::live));
    all.extend(p);

    let p = atto_process::registry();
    drift.extend(atto_params::verify(&p, &atto_process::live));
    all.extend(p);

    let p = atto_device::registry();
    drift.extend(atto_params::verify(&p, &atto_device::live));
    all.extend(p);

    let c = atto_params::coverage(&all);
    println!("AttoFab parameter manifest");
    println!("  registered  {}", c.total);
    println!("  measured    {}", c.measured);
    println!("  published   {}", c.published);
    println!("  fitted      {}", c.fitted);
    println!("  assumed     {}", c.assumed);
    println!("  drift       {}", drift.len());
    for d in &drift {
        println!("    {d}");
    }
    if c.measured == 0 {
        println!("\n  note: no constant is measured. Every model number here is published,");
        println!("        fitted or assumed, which is why the sign-off gate says VERIFIED ONLY.");
    }
    if !drift.is_empty() {
        eprintln!("\nthe manifest has drifted from the code");
        std::process::exit(1);
    }
}
