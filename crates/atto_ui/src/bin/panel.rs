// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Drive the presentation layer with the real process models.
//!
//! This is the headless equivalent of clicking through the recipe: commit a
//! deposition, read the measurement card with its provenance, change a
//! parameter and watch the status derive itself, then undo.

use atto_units::{Mm, Nm, Um};
use atto_geom::column::{grow_film, ColumnGrid, ProcessWindow};
use atto_geom::svo::{Material, Wafer};
use atto_measure::{film_thickness_nm, step_height_nm, Measurement};
use atto_ui::{compare, measure_rows, rows_are_honest, Session, StageId, StageResult, StageStatus};

fn main() {
    let mut session = Session::new();
    session.set_param(1, "thickness_nm", 100.0).expect("param");
    session.set_param(1, "coverage", 0.0).expect("param");

    // the evaluator is the real geometry, not a stub
    let eval = |spec: &atto_ui::StageSpec, fp: u64| -> Result<StageResult, String> {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).map_err(|e| format!("{e:?}"))?;
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        let g = ColumnGrid::sample(&w, &win).map_err(|e| format!("{e:?}"))?;
        let t = spec.get("thickness_nm").ok_or("thickness_nm is missing")?;
        let cov = spec.get("coverage").unwrap_or(0.0);
        let (out, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(t), cov).map_err(|e| format!("{e:?}"))?;
        let film = match film_thickness_nm(&out, Material::SiliconDioxide) {
            Some(v) => Measurement::against("Film thickness", "nm", v, t, 0.1 * t, 1.0),
            None => Measurement::missing("Film thickness", "nm", "the window carries no film"),
        };
        let step = match step_height_nm(&out) {
            Some(v) => Measurement::against("Step height", "nm", v, 0.0, 1.0, 0.5),
            None => Measurement::missing("Step height", "nm", "no active columns"),
        };
        Ok(StageResult {
            measurements: vec![film, step],
            model: "atto_geom::grow_film".into(),
            parameters: spec.params.iter().map(|(k, _)| k.clone()).collect(),
            inputs_fingerprint: fp,
        })
    };

    session.commit(1, &eval).expect("commit");

    println!("RECIPE");
    for r in session.recipe_rows() {
        let mark = if r.selected { ">" } else { " " };
        println!("  {mark} {}. {:<18} {}", r.index + 1, r.label, r.status.word());
    }

    let rows = measure_rows(session.result(1).expect("a result"));
    println!("\nMEASUREMENT CARD (stage {})", StageId::Deposition.label());
    for r in &rows {
        println!("  {:<16} {:<26} [{}]", r.label, r.display(), r.status.label());
    }
    println!("  honest rows: {}", rows_are_honest(&rows));
    if let Some(p) = rows.first() {
        println!("\nPROVENANCE  {}", p.provenance.summary());
        if let Some(s) = p.provenance.uncertainty {
            println!("            sigma {s:.2} nm");
        }
    }

    let before = rows.clone();
    println!("\nchange thickness 100 -> 140 nm");
    session.set_param(1, "thickness_nm", 140.0).expect("param");
    println!("  stage 1 status: {}", session.status(1).word());
    println!("  stage 2 status: {}", session.status(2).word());

    session.undo().expect("undo");
    println!("undo -> stage 1 status: {}", session.status(1).word());
    assert_eq!(session.status(1), StageStatus::Committed);

    session.set_param(1, "thickness_nm", 140.0).expect("param");
    session.commit(1, &eval).expect("commit");
    let after = measure_rows(session.result(1).expect("a result"));
    println!("\nCOMPARE (before -> after)");
    for d in compare(&before, &after) {
        println!("  {:<16} {}", d.label, d.display());
    }
}
