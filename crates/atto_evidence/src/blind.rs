// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Blind validation datasets (build PRD R-01, R-03, K-03).
//!
//! Two rules the audit found broken: a dataset that cannot be loaded is a
//! blocker rather than a silent omission, and a dataset's own tolerance is never
//! rewritten — if it is looser than the criterion the dataset is rejected.

use crate::json;

#[derive(Debug, Clone, PartialEq)]
/// One blind case: a prediction and the reference it is compared against.
pub struct BlindCase {
    pub id: String,
    pub reference: f64,
    pub predicted: f64,
    pub unit: String,
}

#[derive(Debug, Clone, PartialEq)]
/// A named set of blind cases with its declared tolerance.
pub struct BlindDataset {
    pub name: String,
    pub tolerance_pct: f64,
    pub cases: Vec<BlindCase>,
    /// True when the name marks it as fabricated. Synthetic data is shown but
    /// never counts as validation evidence.
    pub synthetic: bool,
    /// The commit or run that produced the predictions. A dataset without one
    /// cannot be checked against the code that made it.
    pub generated_from: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
/// How a dataset compared.
pub struct BlindReport {
    pub name: String,
    pub synthetic: bool,
    pub tolerance_used: f64,
    pub case_count: usize,
    pub passed: usize,
    pub max_error_pct: f64,
    pub overall_pass: bool,
    pub generated_from: Option<String>,
}

/// Whether a dataset name marks its contents as generated rather than measured.
pub fn is_synthetic(name: &str) -> bool {
    let d = name.to_uppercase();
    d.contains("SYNTHETIC") || d.contains("DEMO")
}

/// Load a dataset. Every failure is an error the caller must report, never a
/// dataset that quietly disappears.
pub fn load(text: &str) -> Result<BlindDataset, String> {
    let root = json::parse(text)?;
    let name = root.get("dataset").and_then(|d| d.as_str()).unwrap_or("").to_string();
    if name.trim().is_empty() {
        return Err("dataset has no name".into());
    }
    let tolerance_pct = root.get("max_relative_error_percent").and_then(|t| t.as_f64()).ok_or("dataset has no tolerance")?;
    if !(tolerance_pct.is_finite() && tolerance_pct > 0.0) {
        return Err("dataset tolerance must be finite and positive".into());
    }
    let cases_json = root.get("cases").and_then(|c| c.as_array()).ok_or("dataset has no cases")?;
    if cases_json.is_empty() {
        return Err("dataset has no cases".into());
    }
    let mut cases = Vec::with_capacity(cases_json.len());
    for c in cases_json {
        let id = c.get("id").and_then(|v| v.as_str()).unwrap_or("").to_string();
        let reference = c.get("reference").and_then(|v| v.as_f64()).ok_or_else(|| format!("case '{id}' has no reference"))?;
        let predicted = c.get("predicted").and_then(|v| v.as_f64()).ok_or_else(|| format!("case '{id}' has no prediction"))?;
        let unit = c.get("unit").and_then(|v| v.as_str()).unwrap_or("").to_string();
        if id.trim().is_empty() || unit.trim().is_empty() {
            return Err(format!("case '{id}' needs an id and a unit"));
        }
        if !(reference.is_finite() && predicted.is_finite()) {
            return Err(format!("case '{id}' has a non-finite value"));
        }
        if reference == 0.0 {
            return Err(format!("case '{id}' has a zero reference, so a relative error is undefined"));
        }
        cases.push(BlindCase { id, reference, predicted, unit });
    }
    let generated_from = root.get("generated_from").and_then(|g| g.as_str()).map(|s| s.to_string());
    Ok(BlindDataset { synthetic: is_synthetic(&name), name, tolerance_pct, cases, generated_from })
}

/// Judge a dataset against a sign-off criterion.
///
/// The dataset's own tolerance is used as written. If it is looser than the
/// criterion, that is an error: silently tightening it would change the
/// evidence, and silently accepting it would weaken the gate.
pub fn evaluate(ds: &BlindDataset, criterion_pct: f64) -> Result<BlindReport, String> {
    if !(criterion_pct.is_finite() && criterion_pct > 0.0) {
        return Err("criterion must be finite and positive".into());
    }
    if ds.tolerance_pct > criterion_pct + 1e-12 {
        return Err(format!(
            "dataset '{}' declares a {:.1}% tolerance, looser than the {:.1}% sign-off criterion",
            ds.name, ds.tolerance_pct, criterion_pct
        ));
    }
    let mut passed = 0;
    let mut max_error: f64 = 0.0;
    for c in &ds.cases {
        let e = ((c.predicted - c.reference) / c.reference).abs() * 100.0;
        if e <= ds.tolerance_pct {
            passed += 1;
        }
        max_error = max_error.max(e);
    }
    Ok(BlindReport {
        name: ds.name.clone(),
        synthetic: ds.synthetic,
        tolerance_used: ds.tolerance_pct,
        case_count: ds.cases.len(),
        passed,
        max_error_pct: max_error,
        overall_pass: passed == ds.cases.len(),
        generated_from: ds.generated_from.clone(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const OK: &str = r#"{"dataset":"FAB-A lot 42","max_relative_error_percent":5.0,
        "generated_from":"commit abc123",
        "cases":[{"id":"a","reference":100,"predicted":102,"unit":"nm"},
                 {"id":"b","reference":50,"predicted":48.5,"unit":"nm"}]}"#;

    #[test]
    fn a_valid_dataset_loads_and_passes() {
        let ds = load(OK).unwrap();
        assert_eq!(ds.cases.len(), 2);
        assert!(!ds.synthetic);
        assert_eq!(ds.generated_from.as_deref(), Some("commit abc123"));
        let r = evaluate(&ds, 5.0).unwrap();
        assert!(r.overall_pass);
        assert!((r.max_error_pct - 3.0).abs() < 1e-9);
    }

    #[test]
    fn malformed_input_is_an_error_not_a_missing_dataset() {
        assert!(load("{").is_err());
        assert!(load(r#"{"dataset":"x"}"#).is_err());
        assert!(load(r#"{"dataset":"x","max_relative_error_percent":5,"cases":[]}"#).is_err());
        assert!(load(r#"{"dataset":"x","max_relative_error_percent":5,"cases":[{"id":"a","reference":0,"predicted":1,"unit":"nm"}]}"#).is_err());
    }

    #[test]
    fn a_looser_tolerance_is_rejected_not_clamped() {
        let ds = load(OK).unwrap();
        let e = evaluate(&ds, 2.0).unwrap_err();
        assert!(e.contains("looser"), "{e}");
        // and the dataset is unchanged
        assert!((ds.tolerance_pct - 5.0).abs() < 1e-12);
    }

    #[test]
    fn a_stricter_tolerance_is_honoured() {
        let strict = OK.replace("\"max_relative_error_percent\":5.0", "\"max_relative_error_percent\":1.0");
        let ds = load(&strict).unwrap();
        let r = evaluate(&ds, 5.0).unwrap();
        assert!(!r.overall_pass, "2% error against a 1% tolerance fails");
        assert_eq!(r.passed, 0);
    }

    #[test]
    fn synthetic_names_are_detected_conservatively() {
        assert!(is_synthetic("SYNTHETIC_DEMO_NOT_EXPERIMENTAL_DATA"));
        assert!(is_synthetic("my demo set"));
        assert!(!is_synthetic("FAB-A lot 42 CD-SEM"));
    }
}
