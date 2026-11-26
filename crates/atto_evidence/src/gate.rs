// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The sign-off gate (build PRD R-01 to R-06).
//!
//! Readiness is decided from checks that can fail. Every way of losing evidence
//! — a file that will not load, a hash that will not verify, a signature that is
//! absent — is a blocker, because the failure mode this gate exists to prevent
//! is reporting a better status than the evidence supports.

use crate::blind::BlindReport;
use crate::calibration::{CalibrationAudit, Evidence, HashStatus};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The outcome of one check.
pub enum CheckStatus {
    Pass,
    Fail,
    Blocked,
}

#[derive(Debug, Clone, PartialEq)]
/// One check's result, with the reason it gave.
pub struct Finding {
    pub id: String,
    pub status: CheckStatus,
    pub message: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The strongest claim the evidence supports.
pub enum Readiness {
    /// The gate could not run.
    Blocked,
    /// A check failed.
    Failed,
    /// Everything computable passes, but there is no measured evidence.
    VerifiedOnly,
    /// Verification passes and measured data is within tolerance.
    Passed,
}

impl Readiness {
    /// The name used in a report.
    pub fn label(self) -> &'static str {
        match self {
            Readiness::Blocked => "BLOCKED",
            Readiness::Failed => "FAILED",
            Readiness::VerifiedOnly => "VERIFIED ONLY (not validated)",
            Readiness::Passed => "PASSED",
        }
    }
}

#[derive(Debug, Clone)]
/// The gate's verdict and everything behind it.
pub struct GateReport {
    pub readiness: Readiness,
    pub blockers: Vec<String>,
    pub checks: Vec<Finding>,
    pub blind: Vec<BlindReport>,
    pub evidence: Evidence,
}

#[derive(Debug, Clone, Copy)]
/// What the gate is asked to judge.
pub struct GateInputs<'a> {
    pub checks: &'a [Finding],
    pub calibration: &'a CalibrationAudit,
    pub blind: &'a [BlindReport],
    /// Datasets that could not be loaded. Never dropped.
    pub blind_errors: &'a [String],
    pub criterion_pct: f64,
}

/// Run every check and report the strongest claim they support.
pub fn evaluate(inputs: GateInputs<'_>) -> GateReport {
    let mut blockers = Vec::new();

    // R-06: the gate must have run at all.
    if inputs.checks.is_empty() {
        blockers.push("the verification suite produced no checks".to_string());
        return GateReport {
            readiness: Readiness::Blocked,
            blockers,
            checks: inputs.checks.to_vec(),
            blind: inputs.blind.to_vec(),
            evidence: inputs.calibration.evidence,
        };
    }

    // R-01: evidence that cannot be read is a blocker.
    for e in inputs.blind_errors {
        blockers.push(format!("blind dataset could not be used: {e}"));
    }

    // R-02: a file that is missing, unreadable or altered is a blocker.
    for f in &inputs.calibration.files {
        match &f.hash {
            HashStatus::Verified => {}
            HashStatus::Mismatch => blockers.push(format!("{} no longer matches its recorded SHA-256", f.path)),
            HashStatus::Unreadable(why) => blockers.push(format!("calibration file is unusable: {why}")),
        }
    }

    for c in inputs.checks {
        match c.status {
            CheckStatus::Fail => blockers.push(format!("verification failed: {} — {}", c.id, c.message)),
            // a check that could not run has verified nothing
            CheckStatus::Blocked => blockers.push(format!("verification could not run: {} — {}", c.id, c.message)),
            CheckStatus::Pass => {}
        }
    }

    // R-03: a dataset looser than the criterion was rejected by the loader, so
    // it arrives as an error above. A dataset that failed is a failure.
    let measured: Vec<&BlindReport> = inputs.blind.iter().filter(|b| !b.synthetic).collect();
    for b in &measured {
        if !b.overall_pass {
            blockers.push(format!("measured dataset '{}' exceeds {:.1}% (worst {:.2}%)", b.name, b.tolerance_used, b.max_error_pct));
        }
    }

    if inputs.calibration.evidence != Evidence::Experimental {
        blockers.push(format!(
            "calibration pack is {:?}: source, conditions and uncertainty are not all verifiable, or no external signature",
            inputs.calibration.evidence
        ));
    }
    if measured.is_empty() {
        blockers.push("no measured (non-synthetic) blind-validation dataset is available".to_string());
    }
    for b in inputs.blind.iter().filter(|b| b.synthetic) {
        blockers.push(format!("dataset '{}' is synthetic and does not count as validation", b.name));
    }

    let hard_fail = inputs.checks.iter().any(|c| c.status == CheckStatus::Fail)
        || inputs.checks.iter().any(|c| c.status == CheckStatus::Blocked)
        || !inputs.blind_errors.is_empty()
        || inputs.calibration.files.iter().any(|f| !matches!(f.hash, HashStatus::Verified))
        || measured.iter().any(|b| !b.overall_pass);

    let readiness = if hard_fail {
        Readiness::Failed
    } else if measured.is_empty() || inputs.calibration.evidence != Evidence::Experimental {
        Readiness::VerifiedOnly
    } else {
        Readiness::Passed
    };

    GateReport {
        readiness,
        blockers,
        checks: inputs.checks.to_vec(),
        blind: inputs.blind.to_vec(),
        evidence: inputs.calibration.evidence,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::blind::{load, evaluate as eval_blind};
    use crate::calibration::audit;

    const MANIFEST: &str = r#"{"files":[{"path":"a.json","sha256":"aaa","evidence_status":"unverified",
        "source":null,"conditions":null,"uncertainty":null}]}"#;
    const GOOD_MANIFEST: &str = r#"{"files":[{"path":"a.json","sha256":"aaa","evidence_status":"experimental",
        "source":{"kind":"measured","id":"lot 42"},"conditions":{"tool":"A"},"uncertainty":2.5}]}"#;
    const SIG: &str = r#"{"signed_by":"L. Gundoju","pack":"data/calibration"}"#;
    const DS: &str = r#"{"dataset":"FAB-A lot 42","max_relative_error_percent":5.0,
        "cases":[{"id":"a","reference":100,"predicted":102,"unit":"nm"}]}"#;
    const SYNTH: &str = r#"{"dataset":"SYNTHETIC_DEMO","max_relative_error_percent":5.0,
        "cases":[{"id":"a","reference":100,"predicted":102,"unit":"nm"}]}"#;

    fn digest(b: &[u8]) -> String { String::from_utf8_lossy(b).to_string() }
    fn read(_p: &str) -> Option<Vec<u8>> { Some(b"aaa".to_vec()) }
    fn pass() -> Finding { Finding { id: "C1".into(), status: CheckStatus::Pass, message: "ok".into() } }
    fn fail() -> Finding { Finding { id: "C2".into(), status: CheckStatus::Fail, message: "bad".into() } }

    #[test]
    fn bundled_style_evidence_is_verified_only() {
        let cal = audit(MANIFEST, &read, &digest, None).unwrap();
        let b = eval_blind(&load(SYNTH).unwrap(), 5.0).unwrap();
        let r = evaluate(GateInputs { checks: &[pass()], calibration: &cal, blind: &[b], blind_errors: &[], criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::VerifiedOnly);
        assert!(r.blockers.iter().any(|x| x.contains("synthetic")));
        assert!(r.blockers.iter().any(|x| x.contains("measured")));
    }

    #[test]
    fn a_failed_check_fails_the_gate() {
        let cal = audit(GOOD_MANIFEST, &read, &digest, Some(SIG)).unwrap();
        let b = eval_blind(&load(DS).unwrap(), 5.0).unwrap();
        let r = evaluate(GateInputs { checks: &[fail()], calibration: &cal, blind: &[b], blind_errors: &[], criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::Failed);
    }

    #[test]
    fn an_unloadable_dataset_fails_the_gate() {
        let cal = audit(GOOD_MANIFEST, &read, &digest, Some(SIG)).unwrap();
        let b = eval_blind(&load(DS).unwrap(), 5.0).unwrap();
        let errs = vec!["wafer_006.json: expected ':'".to_string()];
        let r = evaluate(GateInputs { checks: &[pass()], calibration: &cal, blind: &[b], blind_errors: &errs, criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::Failed, "losing a dataset must not be silent");
        assert!(r.blockers.iter().any(|x| x.contains("wafer_006.json")));
    }

    #[test]
    fn an_unreadable_calibration_file_fails_the_gate() {
        let read_missing = |p: &str| if p == "a.json" { None } else { Some(b"aaa".to_vec()) };
        let cal = audit(GOOD_MANIFEST, &read_missing, &digest, Some(SIG)).unwrap();
        let b = eval_blind(&load(DS).unwrap(), 5.0).unwrap();
        let r = evaluate(GateInputs { checks: &[pass()], calibration: &cal, blind: &[b], blind_errors: &[], criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::Failed);
        assert!(r.blockers.iter().any(|x| x.contains("unusable")));
    }

    #[test]
    fn everything_clean_reaches_passed() {
        let cal = audit(GOOD_MANIFEST, &read, &digest, Some(SIG)).unwrap();
        let b = eval_blind(&load(DS).unwrap(), 5.0).unwrap();
        let r = evaluate(GateInputs { checks: &[pass()], calibration: &cal, blind: &[b], blind_errors: &[], criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::Passed);
        assert!(r.blockers.is_empty(), "{:?}", r.blockers);
    }

    #[test]
    fn an_empty_check_list_blocks_rather_than_passes() {
        let cal = audit(GOOD_MANIFEST, &read, &digest, Some(SIG)).unwrap();
        let r = evaluate(GateInputs { checks: &[], calibration: &cal, blind: &[], blind_errors: &[], criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::Blocked);
    }

    #[test]
    fn a_failing_measured_dataset_fails_the_gate() {
        let cal = audit(GOOD_MANIFEST, &read, &digest, Some(SIG)).unwrap();
        let bad = DS.replace("\"predicted\":102", "\"predicted\":130");
        let b = eval_blind(&load(&bad).unwrap(), 5.0).unwrap();
        let r = evaluate(GateInputs { checks: &[pass()], calibration: &cal, blind: &[b], blind_errors: &[], criterion_pct: 5.0 });
        assert_eq!(r.readiness, Readiness::Failed);
    }
}
