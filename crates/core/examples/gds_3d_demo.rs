//! Phase 3 demonstration: parses a GDSII layout (two rectangles on
//! different layers - an "active area" opening and a "contact" opening),
//! rasterizes each into a lithography mask, and runs a 3D process flow
//! (masked oxidation, masked implant, 3D anneal) on a `Wafer3d`. Exports
//! JSON (material + species + a surface heightmap) for the WebGL
//! topography viewer (`web/topography_3d.html`).
//!
//! This demo builds its own small synthetic GDSII file in memory (two
//! BOUNDARY rectangles on layers 1 and 2) rather than requiring an
//! external .gds file, so it's runnable with zero external assets - the
//! same rationale as `gds::tests::synthetic_gds_with_rectangle`, just at
//! demo scale instead of unit-test scale. To drive this from a real GDS
//! file (e.g. an OpenLane-synthesized standard cell), swap
//! `synthetic_layout_gds()` for `std::fs::read("your_layout.gds")`.
//!
//! Run with: `cargo run --example gds_3d_demo -p attofab-core`

use attofab_core::gds::{encode_real8, parse_gds, rasterize_polygon_um};
use attofab_core::grid3d::Wafer3d;
use attofab_core::materials::{Ambient, Dopant};
use std::fs;

fn record(rtype: u8, dtype: u8, data: &[u8]) -> Vec<u8> {
    let len = 4 + data.len();
    let mut out = Vec::with_capacity(len);
    out.extend_from_slice(&(len as u16).to_be_bytes());
    out.push(rtype);
    out.push(dtype);
    out.extend_from_slice(data);
    out
}

fn xy_record(points_dbu: &[(i32, i32)]) -> Vec<u8> {
    let mut data = Vec::new();
    for &(x, y) in points_dbu {
        data.extend_from_slice(&x.to_be_bytes());
        data.extend_from_slice(&y.to_be_bytes());
    }
    record(0x10, 0x03, &data)
}

/// Builds a small synthetic GDSII library: layer 1 = a large "active area"
/// rectangle, layer 2 = a small "contact" rectangle inside it. Coordinates
/// in nanometers (1 dbu = 1nm).
fn synthetic_layout_gds() -> Vec<u8> {
    let mut buf = Vec::new();
    buf.extend(record(0x00, 0x02, &600i16.to_be_bytes())); // HEADER
    buf.extend(record(0x01, 0x02, &[0u8; 24])); // BGNLIB
    buf.extend(record(0x02, 0x06, b"DEMOLIB\0")); // LIBNAME

    let mut units_data = Vec::new();
    units_data.extend_from_slice(&encode_real8(0.001));
    units_data.extend_from_slice(&encode_real8(1e-9));
    buf.extend(record(0x03, 0x05, &units_data)); // UNITS

    buf.extend(record(0x05, 0x02, &[0u8; 24])); // BGNSTR
    buf.extend(record(0x06, 0x06, b"TOP\0")); // STRNAME

    // Layer 1: active area, 600nm x 600nm, centered in a 800nm x 800nm die.
    buf.extend(record(0x08, 0x00, &[])); // BOUNDARY
    buf.extend(record(0x0d, 0x02, &1i16.to_be_bytes())); // LAYER
    buf.extend(xy_record(&[(100, 100), (700, 100), (700, 700), (100, 700), (100, 100)]));
    buf.extend(record(0x11, 0x00, &[])); // ENDEL

    // Layer 2: contact opening, 150nm x 150nm, centered.
    buf.extend(record(0x08, 0x00, &[]));
    buf.extend(record(0x0d, 0x02, &2i16.to_be_bytes()));
    buf.extend(xy_record(&[(325, 325), (475, 325), (475, 475), (325, 475), (325, 325)]));
    buf.extend(record(0x11, 0x00, &[]));

    buf.extend(record(0x07, 0x00, &[])); // ENDSTR
    buf.extend(record(0x04, 0x00, &[])); // ENDLIB
    buf
}

fn main() {
    println!("=== AttoFab Phase 3: GDSII-driven 3D demo ===\n");

    let gds_bytes = synthetic_layout_gds();
    let lib = parse_gds(&gds_bytes).expect("synthetic GDS should parse");
    println!("Parsed GDSII: {} structure(s), {:.4} um/dbu", lib.structures.len(), lib.um_per_dbu);

    let top = &lib.structures[0];
    for p in &top.polygons {
        println!("  layer {} polygon: {} vertices", p.layer, p.points_dbu.len());
    }

    // Grid: 800nm x 800nm die at 10nm/voxel laterally, 5nm/voxel in depth.
    let nx = 80;
    let ny = 80;
    let nz = 60;
    let dx_um = 0.01;
    let dy_um = 0.01;
    let dz_um = 0.005;

    // Rasterize each GDS layer into a lithography mask at the wafer grid resolution.
    let mut active_mask = vec![false; nx * ny];
    let mut contact_mask = vec![false; nx * ny];
    for p in &top.polygons {
        let points_um: Vec<(f64, f64)> =
            p.points_dbu.iter().map(|&(x, y)| (x as f64 * lib.um_per_dbu, y as f64 * lib.um_per_dbu)).collect();
        match p.layer {
            1 => rasterize_polygon_um(&points_um, nx, ny, dx_um, dy_um, &mut active_mask),
            2 => rasterize_polygon_um(&points_um, nx, ny, dx_um, dy_um, &mut contact_mask),
            _ => {}
        }
    }
    let active_count = active_mask.iter().filter(|&&b| b).count();
    let contact_count = contact_mask.iter().filter(|&&b| b).count();
    println!("\nRasterized masks: active area = {active_count} voxel-columns, contact = {contact_count} voxel-columns");

    // Field oxide grows everywhere EXCEPT the active area (LOCOS pattern,
    // same concept as the 2D demo, now driven by a real GDS layer instead
    // of a hand-computed range).
    let field_mask: Vec<bool> = active_mask.iter().map(|&b| !b).collect();

    let mut wafer = Wafer3d::new(nx, ny, nz, dx_um, dy_um, dz_um, Dopant::Boron, 1e15);
    wafer.oxidize(1000.0, 1.5, Ambient::Wet, Some(&field_mask));
    wafer.implant(Dopant::Phosphorus, 1e15, 60.0, Some(&active_mask));
    wafer.anneal(1000.0, 15.0, 0.2);

    // Cut a contact trench through the field oxide at the GDS-layer-2
    // opening, straight down 15 voxels (75nm) - this is what actually
    // gives the topology real height variation for the WebGL viewer:
    // oxidation replaces material in-place (silicon -> oxide) rather than
    // modeling the ~1.9x volume expansion that would push the surface
    // upward, so without an etch step the surface height is flat
    // everywhere despite the material varying. A real trench (Void
    // region) is what a "topography" viewer actually has something to
    // show.
    wafer.etch_anisotropic(15, &contact_mask);

    println!("\nProcess log:");
    for entry in &wafer.process_log {
        println!("  {entry}");
    }

    let field_ox_col = 5 * nx + 5; // a corner column, should be field oxide
    let active_ox_col = (ny / 2) * nx + nx / 2; // center column, should be active (no field oxide)
    println!(
        "\nField oxide (corner column): {:.1} nm",
        wafer.oxide_thickness_um[field_ox_col] * 1000.0
    );
    println!(
        "Field oxide (active center column): {:.1} nm",
        wafer.oxide_thickness_um[active_ox_col] * 1000.0
    );

    let json = serde_json::to_string(&wafer).expect("serialization should not fail");
    fs::create_dir_all("output").ok();
    fs::write("output/gds_3d_demo.json", &json).expect("failed to write output JSON");
    println!("\nExported {} bytes to output/gds_3d_demo.json", json.len());
}
