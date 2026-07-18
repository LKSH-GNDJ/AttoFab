//! LOCOS-style Phase 2 demonstration.
//!
//! Runs a small but complete 2D process flow exercising every Phase 2
//! solver, then writes a JSON snapshot of the final wafer state to
//! `output/locos_2d.json` for the canvas cross-section renderer
//! (`web/index.html`).
//!
//! Flow:
//!   1. Blanket photoresist coat.
//!   2. Mask + expose + develop, opening "field" regions while a central
//!      "active area" stays resist-protected (classic LOCOS masking).
//!   3. Masked thermal oxidation - field oxide grows only in the open
//!      (unprotected) regions ("oxide mountains"); the active area
//!      underneath the resist stays bare silicon.
//!   4. Strip resist (develop already removed it from field regions;
//!      here we clear what's left in the active region to reveal it).
//!   5. Masked Phosphorus implant into the now-exposed active area only.
//!   6. 2D anneal - dopant diffuses both down *and* sideways, partially
//!      spreading under the adjacent field oxide edge.
//!   7. A second lithography + anisotropic etch cycle cuts a narrow
//!      contact trench straight down into the field oxide (no undercut).
//!
//! Run with: `cargo run --example locos_2d_demo -p attofab-core`

use attofab_core::etch::{deposit, etch_anisotropic};
use attofab_core::grid2d::Wafer2d;
use attofab_core::lithography::{develop, expose};
use attofab_core::materials::{Ambient, Dopant, Material};
use std::fs;
use std::io::Write;

fn main() {
    let nx = 80;
    let ny = 100;
    let dx_um = 0.01;
    let dy_um = 0.01;

    let mut wafer = Wafer2d::new(nx, ny, dx_um, dy_um, Dopant::Boron, 1e15);

    // --- Step 1: coat resist ---
    wafer.spin_photoresist(10);

    // --- Step 2: mask - protect the central third (the "active area") ---
    let active_lo = nx / 3;
    let active_hi = 2 * nx / 3;
    let field_open: Vec<bool> = (0..nx).map(|x| x < active_lo || x >= active_hi).collect();

    expose(&mut wafer, &field_open);
    develop(&mut wafer);

    // --- Step 3: LOCOS field oxidation (masked - only field regions grow oxide) ---
    wafer.oxidize(1000.0, 2.0, Ambient::Wet, Some(&field_open));

    // --- Step 4: strip remaining resist from the active area ---
    for m in wafer.material.iter_mut() {
        if *m == Material::Photoresist {
            *m = Material::Void;
        }
    }

    // --- Step 5: masked implant into the active area only ---
    let active_mask: Vec<bool> = (0..nx).map(|x| x >= active_lo && x < active_hi).collect();
    wafer.implant(Dopant::Phosphorus, 1e15, 60.0, Some(&active_mask));

    // --- Step 6: anneal (2D diffusion - spreads down AND sideways) ---
    wafer.anneal(1000.0, 20.0, 0.2);

    // --- Step 7: contact trench - a second litho + anisotropic etch cycle ---
    wafer.spin_photoresist(6);
    let trench_open: Vec<bool> = (0..nx).map(|x| x >= nx / 2 - 2 && x < nx / 2 + 2).collect();
    expose(&mut wafer, &trench_open);
    develop(&mut wafer);
    etch_anisotropic(&mut wafer, 12);

    // --- Step 8: fill the trench with metal (contact plug) ---
    deposit(&mut wafer, Material::Metal, 12);

    // --- Step 9: strip resist (ash/clean) - reveal the real topology for
    // the renderer instead of showing a blanket resist coat on top ---
    for m in wafer.material.iter_mut() {
        if *m == Material::Photoresist || *m == Material::PhotoresistExposed {
            *m = Material::Void;
        }
    }

    // --- Export ---
    fs::create_dir_all("output").expect("create output dir");
    let json = serde_json::to_string_pretty(&wafer).expect("serialize wafer");
    let mut f = fs::File::create("output/locos_2d.json").expect("create output file");
    f.write_all(json.as_bytes()).expect("write output file");

    println!("=== AttoFab Phase 2: LOCOS 2D Demo ===");
    println!("Grid: {nx} x {ny} voxels ({:.2} x {:.2} um)", nx as f64 * dx_um, ny as f64 * dy_um);
    println!("Field oxide thickness (edge column): {:.1} nm", wafer.oxide_thickness_um[0] * 1000.0);
    println!("Active area oxide thickness (should be ~0): {:.4} nm", wafer.oxide_thickness_um[nx / 2] * 1000.0);
    println!("\nProcess log:");
    for entry in &wafer.process_log {
        println!("  - {entry}");
    }
    println!("\nJSON snapshot written to output/locos_2d.json");
}
