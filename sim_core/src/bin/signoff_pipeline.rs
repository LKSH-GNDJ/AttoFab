//! `signoff_pipeline [config.json] [--require-passed]`
//!
//! Runs the verification suite, audits the calibration pack, evaluates blind datasets and
//! writes `signoff.json` + `signoff.html`. Exit code: 1 on `Failed`/`Blocked`, 0 otherwise;
//! with `--require-passed` anything but `Passed` exits 1 (use for release gates).

use sim_core::signoff::{evaluate, status_counts, to_html, Readiness, SignoffConfig};
use std::{env, fs, path::PathBuf, process::ExitCode};

fn main() -> ExitCode {
    let args: Vec<String> = env::args().skip(1).collect();
    let require = args.iter().any(|a| a == "--require-passed");
    let cfg_path = PathBuf::from(args.iter().find(|a| !a.starts_with("--")).cloned().unwrap_or_else(|| "config/signoff.json".into()));
    let cfg = match SignoffConfig::load(&cfg_path) { Ok(c) => c, Err(e) => { eprintln!("cannot read config: {e}"); return ExitCode::from(2); } };
    let root = cfg_path.parent().and_then(|p| p.parent()).map(|p| if p.as_os_str().is_empty() { PathBuf::from(".") } else { p.to_path_buf() }).unwrap_or_else(|| PathBuf::from("."));
    let report = evaluate(&cfg, &root);
    let out = root.join(&cfg.output_dir);
    if let Err(e) = fs::create_dir_all(&out) { eprintln!("cannot create {}: {e}", out.display()); return ExitCode::from(2); }
    let _ = fs::write(out.join("signoff.json"), serde_json::to_string_pretty(&report).unwrap());
    let _ = fs::write(out.join("signoff.html"), to_html(&report));
    let (p, f, b) = status_counts(&report.verification);
    println!("sign-off: {}   (verification: {p} pass, {f} fail, {b} blocked)", report.readiness.label());
    for x in &report.blockers { println!("  - {x}"); }
    println!("report: {}", out.join("signoff.html").display());
    match report.readiness {
        Readiness::Failed | Readiness::Blocked => ExitCode::from(1),
        Readiness::VerifiedOnly if require => ExitCode::from(1),
        _ => ExitCode::SUCCESS,
    }
}
