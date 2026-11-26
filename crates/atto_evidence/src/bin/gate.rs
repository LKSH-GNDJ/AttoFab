// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Run the sign-off gate over a repository's evidence.
//!
//! Usage: gate <repo-root>
//!
//! Reads `<root>/data/calibration/manifest.json` and `<root>/data/blind_test/`,
//! plus the external signature at `<root>/data/calibration.sig.json` if one
//! exists. Every failure to read is reported, never skipped.

use atto_evidence::{blind, calibration, gate, hex_digest, lint};
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        eprintln!("usage: gate <repo-root>");
        std::process::exit(2);
    }
    let root = PathBuf::from(&args[1]);
    let criterion = 5.0;

    // ---- check 1: no silent error swallowing in the domain crates -----------
    let mut violations = 0usize;
    let mut scanned = 0usize;
    // walk the whole tree so the check works for both the `crates/` layout
    // and the `sim_core/`+`sim_ui/` one
    let mut stack = vec![root.clone()];
    while let Some(dir) = stack.pop() {
        if let Ok(entries) = std::fs::read_dir(&dir) {
            for e in entries.flatten() {
                let p = e.path();
                let name = p.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default();
                // same exclusions as the lint test: shells and the lint module
                // itself are not the domain
                if name == "atto_cli" || name == "bin" || name == "lint.rs" {
                    continue;
                }
                if p.is_dir() {
                    if name == "target" || name == ".git" || name == "tests" || name == "benches" {
                        continue;
                    }
                    stack.push(p);
                } else if p.extension().map(|x| x == "rs").unwrap_or(false) {
                    if let Ok(text) = std::fs::read_to_string(&p) {
                        scanned += 1;
                        let v = lint::scan(&text);
                        if !v.is_empty() {
                            println!("  scan        {}:{} {}", p.strip_prefix(&root).unwrap_or(&p).display(), v[0].line, v[0].pattern);
                        }
                        violations += v.len();
                    }
                }
            }
        }
    }
    let checks = vec![gate::Finding {
        id: "R-06 source scan".into(),
        // a scan that found no source has not verified anything
        status: if scanned == 0 {
            gate::CheckStatus::Blocked
        } else if violations == 0 {
            gate::CheckStatus::Pass
        } else {
            gate::CheckStatus::Fail
        },
        message: format!("{scanned} files scanned, {violations} silent fallbacks"),
    }];

    // ---- calibration -------------------------------------------------------
    let cal_dir = root.join("data/calibration");
    let manifest_path = cal_dir.join("manifest.json");
    let override_path = root.join("data/calibration.sig.json");
    let read = |p: &str| std::fs::read(root.join(p)).ok();
    let sig = std::fs::read_to_string(&override_path).ok();
    let calib = match std::fs::read_to_string(&manifest_path) {
        Ok(text) => match calibration::audit(&text, &read, &hex_digest, sig.as_deref()) {
            Ok(a) => a,
            Err(e) => {
                eprintln!("manifest is unusable: {e}");
                std::process::exit(1);
            }
        },
        Err(e) => {
            eprintln!("cannot read {}: {e}", manifest_path.display());
            std::process::exit(1);
        }
    };

    // ---- blind datasets ----------------------------------------------------
    let mut reports = Vec::new();
    let mut errors = Vec::new();
    let blind_dir = root.join("data/blind_test");
    if let Ok(rd) = std::fs::read_dir(&blind_dir) {
        let mut paths: Vec<PathBuf> = rd.flatten().map(|e| e.path()).filter(|p| p.extension().map(|x| x == "json").unwrap_or(false)).collect();
        paths.sort();
        for p in paths {
            let name = p.file_name().map(|f| f.to_string_lossy().into_owned()).unwrap_or_default();
            match std::fs::read_to_string(&p) {
                Err(e) => errors.push(format!("{name}: {e}")),
                Ok(text) => match blind::load(&text) {
                    Err(e) => errors.push(format!("{name}: {e}")),
                    Ok(ds) => match blind::evaluate(&ds, criterion) {
                        Err(e) => errors.push(format!("{name}: {e}")),
                        Ok(r) => reports.push(r),
                    },
                },
            }
        }
    } else {
        errors.push(format!("no blind directory at {}", blind_dir.display()));
    }

    let rep = gate::evaluate(gate::GateInputs {
        checks: &checks,
        calibration: &calib,
        blind: &reports,
        blind_errors: &errors,
        criterion_pct: criterion,
    });

    println!("AttoFab sign-off gate");
    println!("  root        {}", root.display());
    println!("  readiness   {}", rep.readiness.label());
    println!("  calibration {:?} ({} files)", rep.evidence, calib.files.len());
    for c in &checks {
        println!("  check       [{}] {} — {}", if c.status == gate::CheckStatus::Pass { "pass" } else { "FAIL" }, c.id, c.message);
    }
    for f in &calib.files {
        println!("  file        {:<40} {:?}", f.path, f.hash);
    }
    for b in &reports {
        println!("  dataset     {:<42} {}{}", b.name, if b.synthetic { "SYNTHETIC " } else { "" },
            format!("max err {:.2}% tol {:.1}%", b.max_error_pct, b.tolerance_used));
    }
    for e in &errors {
        println!("  ERROR       {e}");
    }
    println!("\n  what stands between this and a pass:");
    if rep.blockers.is_empty() {
        println!("    nothing");
    }
    for b in &rep.blockers {
        println!("    - {b}");
    }
    let code = match rep.readiness {
        gate::Readiness::Passed | gate::Readiness::VerifiedOnly => 0,
        _ => 1,
    };
    std::process::exit(code);
}
