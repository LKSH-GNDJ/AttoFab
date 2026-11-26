// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Command-line runner for AttoFab's blind-validation report framework.
use sim_core::validation::blind::evaluate_file;
use std::path::PathBuf;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: blind_validator <dataset.json> <output-dir>");
        std::process::exit(2);
    }
    match evaluate_file(&PathBuf::from(&args[1]), &PathBuf::from(&args[2])) {
        Ok(report) => {
            println!(
                "dataset={} passed={}/{} mean_abs_relative_error={:.4}% max={:.4}% overall_pass={}",
                report.dataset,
                report.passed,
                report.case_count,
                report.mean_absolute_relative_error_percent,
                report.max_absolute_relative_error_percent,
                report.overall_pass
            );
            if !report.overall_pass {
                std::process::exit(1);
            }
        }
        Err(error) => {
            eprintln!("blind validation failed: {error}");
            std::process::exit(2);
        }
    }
}
