// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! An interactive shell over the panel layer.
//!
//! This is the interface this environment can build and test: no window, no
//! toolkit, but the same session, the same measurements, the same provenance and
//! the same project file a windowed shell would drive. Commands read from stdin:
//!
//!     show                 the recipe and the current stage's card
//!     set <stage> <k> <v>  change a parameter
//!     commit <stage>       run it (real geometry, not a stub)
//!     undo / redo
//!     save <file>          write the project file (the command log)
//!     load <file>          replay it
//!     help / quit

use atto_units::{Mm, Nm, Um};
use atto_geom::column::{grow_film, ColumnGrid, ProcessWindow};
use atto_geom::svo::{Material, Wafer};
use atto_measure::{film_thickness_nm, step_height_nm, Measurement};
use atto_ui::journal::{Command, Journal};
use atto_ui::{measure_rows, rows_are_honest, Session, StageResult, StageStatus};
use std::io::{BufRead, Write};

fn evaluate(spec: &atto_ui::StageSpec, fp: u64) -> Result<StageResult, String> {
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
}

fn show(s: &Session) {
    println!("\nRECIPE");
    for r in s.recipe_rows() {
        println!("  {} {}. {:<18} {}", if r.selected { ">" } else { " " }, r.index + 1, r.label, r.status.word());
    }
    let i = s.selected();
    match s.result(i) {
        Some(res) => {
            let rows = measure_rows(res);
            println!("\nCARD  stage {}", i + 1);
            for r in &rows {
                println!("  {:<16} {:<34} [{}]", r.label, r.display(), r.status.label());
            }
            if let Some(p) = rows.first() {
                println!("  provenance: {}", p.provenance.summary());
            }
            println!("  honest rows: {}", rows_are_honest(&rows));
        }
        None => {
            let why = match s.status(i) {
                StageStatus::NotRun => "this stage has not been committed",
                StageStatus::Stale => "the parameters changed since the last commit",
                StageStatus::Failed => "the last commit failed",
                StageStatus::Blocked => "an upstream stage has not been committed",
                StageStatus::Committed => "committed, but produced no measurements",
            };
            println!("\nCARD  stage {}: nothing to show — {why}", i + 1);
        }
    }
    println!("\nstate {:016x}", s.fingerprint());
}

fn main() {
    let mut s = Session::new();
    let mut j = Journal::new();
    s.set_param(1, "thickness_nm", 100.0).ok();
    s.set_param(1, "coverage", 0.0).ok();
    j.record(Command::SetParam { stage: 1, key: "thickness_nm".into(), value: 100.0 });
    j.record(Command::SetParam { stage: 1, key: "coverage".into(), value: 0.0 });
    println!("AttoFab shell — type `help`, then `show`.");
    let stdin = std::io::stdin();
    loop {
        print!("> ");
        let _ = std::io::stdout().flush();
        let mut line = String::new();
        if stdin.lock().read_line(&mut line).unwrap_or(0) == 0 {
            break;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }
        match parts[0] {
            "quit" | "q" => break,
            "help" | "h" => {
                println!("show | set <stage> <key> <value> | commit <stage> | undo | redo | save <file> | load <file> | quit");
                println!("stages: 1..8, e.g. `set 2 thickness_nm 140` then `commit 2`");
            }
            "show" | "s" => show(&s),
            "set" => {
                if parts.len() < 4 {
                    println!("usage: set <stage> <key> <value>");
                    continue;
                }
                let stage = match parts[1].parse::<usize>() {
                    Ok(v) if v >= 1 && v <= 8 => v - 1,
                    _ => { println!("stage must be 1..8"); continue; }
                };
                let value = match parts[3].parse::<f64>() {
                    Ok(v) => v,
                    Err(_) => { println!("value must be a number"); continue; }
                };
                match s.set_param(stage, parts[2], value) {
                    Ok(()) => {
                        j.record(Command::SetParam { stage, key: parts[2].to_string(), value });
                        println!("stage {}: {} = {}", stage + 1, parts[2], value);
                    }
                    Err(e) => println!("{e:?}"),
                }
            }
            "commit" | "c" => {
                let stage = parts.get(1).and_then(|p| p.parse::<usize>().ok()).map(|v| v.saturating_sub(1)).unwrap_or(s.selected());
                match s.commit(stage, &evaluate) {
                    Ok(Some(_)) => {
                        j.record(Command::Commit { stage });
                        println!("stage {} committed", stage + 1);
                        show(&s);
                    }
                    Ok(None) => println!("stage {} ran and failed", stage + 1),
                    Err(e) => println!("{e:?}"),
                }
            }
            "undo" => match s.undo() {
                Ok(()) => { j.record(Command::Undo); println!("undone"); }
                Err(e) => println!("{e:?}"),
            },
            "redo" => match s.redo() {
                Ok(()) => { j.record(Command::Redo); println!("redone"); }
                Err(e) => println!("{e:?}"),
            },
            "save" => {
                let path = parts.get(1).copied().unwrap_or("project.afp.json");
                match std::fs::write(path, j.to_json()) {
                    Ok(()) => println!("wrote {path} ({} commands)", j.len()),
                    Err(e) => println!("cannot write {path}: {e}"),
                }
            }
            "load" => {
                let path = parts.get(1).copied().unwrap_or("project.afp.json");
                match std::fs::read_to_string(path) {
                    Ok(text) => match Journal::from_json(&text) {
                        Ok(journal) => match journal.replay(&evaluate) {
                            Ok(replayed) => {
                                s = replayed;
                                j = journal;
                                println!("loaded {path} ({} commands), state {:016x}", j.len(), s.fingerprint());
                                show(&s);
                            }
                            Err(e) => println!("replay failed: {e:?}"),
                        },
                        Err(e) => println!("{path} is not a project file: {e}"),
                    },
                    Err(e) => println!("cannot read {path}: {e}"),
                }
            }
            other => println!("unknown command '{other}' — type help"),
        }
    }
}
