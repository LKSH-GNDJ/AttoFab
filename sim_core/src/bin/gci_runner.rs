//! CLI for three-grid Grid Convergence Index calculations.
//!
//! Input JSON shape: {"coarse":{"spacing":4.0,"value":1.16},
//! "medium":{"spacing":2.0,"value":1.04}, "fine":{"spacing":1.0,"value":1.01},
//! "safety_factor":1.25}.
use serde::{Deserialize, Serialize};
use sim_core::validation::gci::{calculate_gci, GciReport, GridResult};
use std::{env, fs, path::PathBuf};

#[derive(Debug, Deserialize)]
struct Input {
    coarse: GridResult,
    medium: GridResult,
    fine: GridResult,
    #[serde(default = "default_safety_factor")]
    safety_factor: f64,
}
fn default_safety_factor() -> f64 {
    1.25
}
#[derive(Serialize)]
struct Output {
    report: GciReport,
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: gci_runner <input.json> <output.json>");
        std::process::exit(2);
    }
    let input_path = PathBuf::from(&args[1]);
    let output_path = PathBuf::from(&args[2]);
    let result = (|| -> Result<(), Box<dyn std::error::Error>> {
        let input: Input = serde_json::from_slice(&fs::read(input_path)?)?;
        let report = calculate_gci(input.coarse, input.medium, input.fine, input.safety_factor)?;
        if let Some(parent) = output_path.parent() {
            if !parent.as_os_str().is_empty() {
                fs::create_dir_all(parent)?;
            }
        }
        fs::write(output_path, serde_json::to_vec_pretty(&Output { report })?)?;
        Ok(())
    })();
    if let Err(error) = result {
        eprintln!("GCI run failed: {error}");
        std::process::exit(1);
    }
}
