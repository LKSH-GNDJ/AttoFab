// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Runs the loader and the process pipeline on a real, hierarchical, multi-scale layout
//! (`tests/data/attofab_master_wafer.gds`: 1 nm database unit, 3 cells, 23 AREF + 9 SREF, features from 50 nm to 5 mm).

use sim_core::geom::{MaskPattern, ProcessWindow};
use sim_core::io::gds::{parse, summarize, FlatPoly};
use sim_core::params::StageParams;
use sim_core::physics::Stage;
use sim_core::pipeline::evaluate_stage;
use sim_core::wafer::WaferState;

fn bytes() -> Vec<u8> { std::fs::read(concat!(env!("CARGO_MANIFEST_DIR"), "/tests/data/attofab_master_wafer.gds")).expect("fixture") }
fn bb(p: &FlatPoly) -> [f64; 4] { p.pts.iter().fold([f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY], |b, v| [b[0].min(v.0), b[1].min(v.1), b[2].max(v.0), b[3].max(v.1)]) }

#[test]
fn the_real_layout_loads_with_its_hierarchy_units_and_layers() {
    let lay = parse(&bytes()).unwrap();
    assert!((lay.db_unit_m - 1e-9).abs() < 1e-18);
    assert_eq!(lay.structures.len(), 3);
    assert_eq!(lay.default_top().as_deref(), Some("ATTOFAB_MASTER_WAFER"));
    let flat = lay.flatten(None).unwrap();
    let layers = summarize(&flat);
    eprintln!("flattened polygons: {}, layers: {}", flat.len(), layers.len());
    assert!(flat.len() > 1000, "arrays and references must expand the 643 stored boundaries: {}", flat.len());
    let l = |a: i32, b: i32| layers.iter().find(|x| x.layer == a && x.datatype == b).unwrap_or_else(|| panic!("layer {a}/{b} missing: {:?}", layers.iter().map(|x| (x.layer, x.datatype)).collect::<Vec<_>>()));
    assert!(l(2, 0).feature_p10_nm <= 60.0, "layer 2/0 carries ~50 nm features: {}", l(2, 0).feature_p10_nm);
    assert!(l(1, 4).feature_p10_nm > 3000.0, "layer 1/4 carries µm-scale features: {}", l(1, 4).feature_p10_nm);
    let (w, h) = l(1, 0).size_nm();
    assert!(w > 9.0e7 && h > 8.0e7, "the die outline spans a ~100 mm wafer: {w} x {h} nm");
}

#[test]
fn region_flatten_matches_the_full_flatten_on_the_real_file() {
    let lay = parse(&bytes()).unwrap();
    let full = lay.flatten(None).unwrap();
    let a = summarize(&full).into_iter().find(|x| x.layer == 2 && x.datatype == 0).unwrap().anchor;
    let region = [a.0 - 3000.0, a.1 - 3000.0, a.0 + 3000.0, a.1 + 3000.0];
    let part = lay.flatten_region(None, Some(region)).unwrap();
    let expect = full.iter().filter(|p| { let b = bb(p); !(b[2] < region[0] || b[0] > region[2] || b[3] < region[1] || b[1] > region[3]) }).count();
    assert_eq!(part.len(), expect);
    assert!(!part.is_empty() && part.len() < full.len() / 10, "{} of {}", part.len(), full.len());
}

fn run_litho(mask: MaskPattern, leaf: f64, levels: u32, window: f64) -> sim_core::pipeline::StageEval {
    let wafer = WaferState::new_wafer_tile(300.0, 500.0, leaf, levels, (0.0, 0.0)).unwrap();
    let mut lp = if let StageParams::Lithography(l) = StageParams::default_for(Stage::Lithography) { l } else { unreachable!() };
    lp.mask = mask;
    lp.resist_thickness_nm = 100.0;
    let win = ProcessWindow { center: [0.0, 0.0], size: [window, window] };
    evaluate_stage(&wafer, &win, &StageParams::Lithography(lp), 1, 1.0).expect("lithography evaluates")
}

#[test]
fn nm_scale_optical_features_print_at_a_fine_cell() {
    let lay = parse(&bytes()).unwrap();
    let full = lay.flatten(None).unwrap();
    let small = full.iter().filter(|p| p.layer == 2 && p.datatype == 0).min_by(|a, b| { let d = |p: &FlatPoly| { let b = bb(p); (b[2] - b[0]).min(b[3] - b[1]) }; d(a).partial_cmp(&d(b)).unwrap() }).unwrap();
    let b = bb(small);
    let c = ((b[0] + b[2]) / 2.0, (b[1] + b[3]) / 2.0);
    let region = [c.0 - 1500.0, c.1 - 1500.0, c.0 + 1500.0, c.1 + 1500.0];
    let (mask, kept) = MaskPattern::from_layout(&lay.flatten_region(None, Some(region)).unwrap(), 2, Some(0), c, 1500.0, false).unwrap();
    assert!(kept >= 1);
    let ev = run_litho(mask, 16.0, 7, 1200.0);
    let (opened, retained) = (ev.outcome.metric("Opened columns").unwrap(), ev.outcome.metric("Retained columns").unwrap());
    assert!(opened > 0.0 && retained > 0.0, "a real feature must produce a pattern: opened {opened}, retained {retained}");
}

#[test]
fn micrometre_scale_features_print_at_an_advised_coarser_cell() {
    let lay = parse(&bytes()).unwrap();
    let full = lay.flatten(None).unwrap();
    let a = summarize(&full).into_iter().find(|x| x.layer == 1 && x.datatype == 4).unwrap();
    let advice = sim_core::scale::advise(a.feature_p10_nm, 100.0);
    assert!(advice.cell_nm <= 50.0 && advice.cell_nm >= 8.0, "{advice:?}");
    let half = advice.window_nm * 0.5 + 500.0;
    let reg = lay.flatten_region(None, Some([a.anchor.0 - half, a.anchor.1 - half, a.anchor.0 + half, a.anchor.1 + half])).unwrap();
    let (mask, kept) = MaskPattern::from_layout(&reg, 1, Some(4), a.anchor, half, false).unwrap();
    assert!(kept >= 1);
    let ev = run_litho(mask, advice.cell_nm, advice.tile_levels, advice.window_nm);
    let (opened, retained) = (ev.outcome.metric("Opened columns").unwrap(), ev.outcome.metric("Retained columns").unwrap());
    assert!(opened > 0.0 && retained > 0.0, "advised scale {advice:?}: opened {opened}, retained {retained}");
}
