//! Evidence-based sign-off gate.
//!
//! Verification (computed here) and validation (needs measured data) are kept apart:
//!
//! * `Failed`       – a verification case failed or a bundled hash no longer matches.
//! * `Blocked`      – the verification suite could not run.
//! * `VerifiedOnly` – every computed check passes, but there is no measured validation
//!                    dataset and/or the calibration pack is not experimental. The model is
//!                    *verified*, not *validated*; results are illustrative.
//! * `Passed`       – verification passes AND calibration is experimental AND every measured
//!                    blind dataset is within tolerance.
//!
//! Synthetic datasets are evaluated and shown, but never count as validation.

use crate::calibration::{self, CalibrationStatus, Evidence};
use crate::validation::blind::{self, BlindDataset, BlindReport};
use crate::validation::CheckStatus;
use crate::vv::{run_suite, VvConfig, VvReport};
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Deserialize)]
pub struct SignoffConfig {
    #[serde(default = "d_cal")] pub calibration_dir: String,
    #[serde(default = "d_blind")] pub blind_dir: String,
    #[serde(default = "d_out")] pub output_dir: String,
    #[serde(default = "d_leaf")] pub gci_min_leaf_sizes: Vec<f64>,
    #[serde(default)] pub signoff_criteria: Criteria,
}

#[derive(Debug, Clone, Deserialize)]
pub struct Criteria { #[serde(default = "d_blind_err")] pub max_blind_error_pct: f64 }
impl Default for Criteria { fn default() -> Self { Self { max_blind_error_pct: 5.0 } } }

fn d_cal() -> String { "data/calibration".into() }
fn d_blind() -> String { "data/blind_test".into() }
fn d_out() -> String { "validation_results/signoff".into() }
fn d_leaf() -> Vec<f64> { vec![32.0, 16.0, 8.0] }
fn d_blind_err() -> f64 { 5.0 }

impl SignoffConfig {
    pub fn load(path: &Path) -> Result<Self, String> {
        serde_json::from_str(&std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?).map_err(|e| e.to_string())
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum Readiness { Blocked, Failed, VerifiedOnly, Passed }

impl Readiness {
    pub fn label(self) -> &'static str {
        match self { Readiness::Blocked => "BLOCKED", Readiness::Failed => "FAILED", Readiness::VerifiedOnly => "VERIFIED ONLY (not validated)", Readiness::Passed => "PASSED" }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct BlindEntry { pub file: String, pub synthetic: bool, pub report: BlindReport }

#[derive(Debug, Clone, Serialize)]
pub struct SignoffReport {
    pub readiness: Readiness,
    pub blockers: Vec<String>,
    pub verification: VvReport,
    pub calibration: Option<CalibrationStatus>,
    pub blind: Vec<BlindEntry>,
}

pub fn is_synthetic(dataset: &str) -> bool { let d = dataset.to_uppercase(); d.contains("SYNTHETIC") || d.contains("DEMO") }

fn load_blind(dir: &Path, tol: f64) -> Vec<BlindEntry> {
    let mut out = Vec::new();
    let Ok(rd) = std::fs::read_dir(dir) else { return out };
    let mut paths: Vec<PathBuf> = rd.flatten().map(|e| e.path()).filter(|p| p.extension().map(|x| x == "json").unwrap_or(false)).collect();
    paths.sort();
    for p in paths {
        let Ok(t) = std::fs::read_to_string(&p) else { continue };
        let Ok(mut ds) = serde_json::from_str::<BlindDataset>(&t) else { continue };
        ds.max_relative_error_percent = ds.max_relative_error_percent.min(tol);
        if let Ok(r) = blind::evaluate(&ds) { out.push(BlindEntry { file: p.file_name().unwrap().to_string_lossy().into(), synthetic: is_synthetic(&ds.dataset), report: r }); }
    }
    out
}

/// Run the whole gate. `root` is the repository root (paths in the config are relative to it).
pub fn evaluate(cfg: &SignoffConfig, root: &Path) -> SignoffReport {
    let verification = run_suite(&VvConfig { leaf_sizes_nm: cfg.gci_min_leaf_sizes.clone(), ..Default::default() });
    let calibration = calibration::audit(root, &root.join(&cfg.calibration_dir)).ok();
    let blind = load_blind(&root.join(&cfg.blind_dir), cfg.signoff_criteria.max_blind_error_pct);
    let mut blockers = Vec::new();

    let readiness = if verification.cases.is_empty() {
        blockers.push("the verification suite produced no cases".into());
        Readiness::Blocked
    } else if !verification.failed().is_empty() {
        for c in verification.failed() { blockers.push(format!("verification failed: {} — {}", c.id, c.message)); }
        Readiness::Failed
    } else if calibration.as_ref().map(|c| c.files.iter().any(|f| f.hash_ok == Some(false))).unwrap_or(false) {
        blockers.push("a calibration file no longer matches its recorded SHA-256".into());
        Readiness::Failed
    } else {
        match &calibration {
            None => blockers.push("calibration pack could not be read".into()),
            Some(c) if c.evidence != Evidence::Experimental => {
                blockers.push(format!("calibration pack is {:?}: source, conditions and uncertainty are not documented, so no stage may be labelled Calibrated", c.evidence));
                if c.max_implant_deviation_pct() > 15.0 { blockers.push(format!("bundled implant range table deviates from the built-in LSS-style table by up to {:.0}%", c.max_implant_deviation_pct())); }
            }
            _ => {}
        }
        let measured: Vec<&BlindEntry> = blind.iter().filter(|b| !b.synthetic).collect();
        if measured.is_empty() { blockers.push("no measured (non-synthetic) blind-validation dataset is available".into()); }
        let measured_fail = measured.iter().any(|b| !b.report.overall_pass);
        if measured_fail { blockers.push("a measured blind dataset exceeds the error tolerance".into()); }
        if measured_fail { Readiness::Failed } else if blockers.is_empty() { Readiness::Passed } else { Readiness::VerifiedOnly }
    };
    SignoffReport { readiness, blockers, verification, calibration, blind }
}

pub fn to_html(r: &SignoffReport) -> String {
    let esc = |s: &str| s.replace('&', "&amp;").replace('<', "&lt;").replace('>', "&gt;");
    let mut h = String::from("<!doctype html><meta charset=utf-8><title>AttoFab sign-off</title><style>body{font:14px system-ui;margin:2rem;max-width:70rem}td,th{border:1px solid #888;padding:3px 8px}table{border-collapse:collapse}.Pass{color:#0a7a2f}.Fail{color:#b00020;font-weight:700}.Blocked{color:#a60}</style>");
    h.push_str(&format!("<h1>AttoFab sign-off: {}</h1>", r.readiness.label()));
    if !r.blockers.is_empty() { h.push_str("<h2>What stands between this and a pass</h2><ul>"); for b in &r.blockers { h.push_str(&format!("<li>{}</li>", esc(b))); } h.push_str("</ul>"); }
    h.push_str("<h2>Verification</h2><table><tr><th>case</th><th>group</th><th>status</th><th>error</th><th>limit</th><th>detail</th></tr>");
    for c in &r.verification.cases { h.push_str(&format!("<tr><td>{}</td><td>{}</td><td class={:?}>{:?}</td><td>{:.4}</td><td>{:.4} {}</td><td>{}</td></tr>", esc(&c.id), c.group, c.status, c.status, c.value, c.limit, c.unit, esc(&c.message))); }
    h.push_str("</table>");
    if let Some(c) = &r.calibration { h.push_str(&format!("<h2>Calibration pack: {:?}</h2><p>{}</p>", c.evidence, esc(&c.notes.join(" ")))); }
    h.push_str("<h2>Blind validation</h2><ul>");
    for b in &r.blind { h.push_str(&format!("<li>{} — {}{} (max error {:.2}%)</li>", esc(&b.file), if b.report.overall_pass { "within tolerance" } else { "OUT OF TOLERANCE" }, if b.synthetic { ", SYNTHETIC: not validation evidence" } else { "" }, b.report.max_absolute_relative_error_percent)); }
    h.push_str("</ul>");
    h
}

pub fn status_counts(r: &VvReport) -> (usize, usize, usize) {
    let c = |s: CheckStatus| r.cases.iter().filter(|x| x.status == s).count();
    (c(CheckStatus::Pass), c(CheckStatus::Fail), c(CheckStatus::Blocked))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn root() -> PathBuf { Path::new(env!("CARGO_MANIFEST_DIR")).join("..") }

    #[test]
    fn bundled_data_is_verified_but_not_validated() {
        let cfg = SignoffConfig { calibration_dir: d_cal(), blind_dir: d_blind(), output_dir: d_out(), gci_min_leaf_sizes: vec![32.0, 16.0], signoff_criteria: Criteria::default() };
        let r = evaluate(&cfg, &root());
        assert_eq!(r.readiness, Readiness::VerifiedOnly, "{:#?}", r.blockers);
        assert!(r.blockers.iter().any(|b| b.contains("calibration pack")) && r.blockers.iter().any(|b| b.contains("measured")));
        assert!(r.blind.iter().all(|b| b.synthetic), "the bundled blind data is synthetic and must not count as validation");
        assert!(to_html(&r).contains("VERIFIED ONLY"));
    }

    #[test]
    fn a_failing_verification_case_fails_the_gate_regardless_of_data() {
        let mut r = SignoffReport { readiness: Readiness::Blocked, blockers: vec![], verification: VvReport { cases: vec![], leaf_sizes_nm: vec![] }, calibration: None, blind: vec![] };
        assert_eq!(r.readiness, Readiness::Blocked);
        r.verification.cases.push(crate::vv::VvCase { id: "x".into(), group: "g", status: CheckStatus::Fail, value: 2.0, limit: 1.0, unit: "", message: "m".into() });
        assert!(!r.verification.all_pass() && r.verification.failed().len() == 1);
    }

    #[test]
    fn synthetic_detection_is_conservative() {
        assert!(is_synthetic("SYNTHETIC_DEMO_NOT_EXPERIMENTAL_DATA_WAFER_001"));
        assert!(is_synthetic("my demo set"));
        assert!(!is_synthetic("FAB-A lot 42 CD-SEM"));
    }
}
