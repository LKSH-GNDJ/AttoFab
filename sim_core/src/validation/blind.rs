// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Blind validation report generation from precomputed predictions and references.
//!
//! Input records must be held out from calibration/tuning. This module does not
//! create or claim to contain real experimental data.

use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;
use thiserror::Error;

/// One scalar comparison for an independently held-out wafer or measurement.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlindCase {
    /// Unique case identifier.
    pub id: String,
    /// Reference measurement and unit.
    pub reference: f64,
    /// Model prediction in the same unit.
    pub predicted: f64,
    /// Unit label, e.g. `nm` or `ions/cm^3`.
    pub unit: String,
}

/// Input file for a blind-validation run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlindDataset {
    /// Dataset name/version.
    pub dataset: String,
    /// Maximum allowed absolute relative error in percent.
    pub max_relative_error_percent: f64,
    /// Held-out cases.
    pub cases: Vec<BlindCase>,
}

/// Per-case result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlindCaseResult {
    /// Case identifier.
    pub id: String,
    /// Reference value.
    pub reference: f64,
    /// Predicted value.
    pub predicted: f64,
    /// Absolute relative error (%).
    pub absolute_relative_error_percent: f64,
    /// Acceptance result for this case.
    pub pass: bool,
    /// Measurement unit.
    pub unit: String,
}

/// Aggregate report.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlindReport {
    /// Dataset name/version.
    pub dataset: String,
    /// Cases evaluated.
    pub case_count: usize,
    /// Number within tolerance.
    pub passed: usize,
    /// Mean absolute relative error (%).
    pub mean_absolute_relative_error_percent: f64,
    /// Worst absolute relative error (%).
    pub max_absolute_relative_error_percent: f64,
    /// Acceptance tolerance (%).
    pub tolerance_percent: f64,
    /// All cases met tolerance.
    pub overall_pass: bool,
    /// Individual case results.
    pub results: Vec<BlindCaseResult>,
}

/// Errors from invalid datasets or report output.
#[derive(Debug, Error)]
pub enum BlindError {
    /// Invalid input dataset.
    #[error("invalid blind dataset: {0}")]
    Invalid(String),
    /// JSON or filesystem error.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// JSON encoding/decoding error.
    #[error(transparent)]
    Json(#[from] serde_json::Error),
}

/// Evaluate a held-out dataset against its declared tolerance.
pub fn evaluate(dataset: &BlindDataset) -> Result<BlindReport, BlindError> {
    if dataset.dataset.trim().is_empty() || dataset.cases.is_empty() {
        return Err(BlindError::Invalid(
            "dataset name and at least one case are required".into(),
        ));
    }
    if !dataset.max_relative_error_percent.is_finite() || dataset.max_relative_error_percent < 0.0 {
        return Err(BlindError::Invalid(
            "tolerance must be finite and non-negative".into(),
        ));
    }

    let mut results = Vec::with_capacity(dataset.cases.len());
    for case in &dataset.cases {
        if case.id.trim().is_empty()
            || case.unit.trim().is_empty()
            || !case.reference.is_finite()
            || !case.predicted.is_finite()
            || case.reference == 0.0
        {
            return Err(BlindError::Invalid(format!(
                "case '{}' has missing or invalid values; reference must be non-zero",
                case.id
            )));
        }
        let error = ((case.predicted - case.reference) / case.reference).abs() * 100.0;
        results.push(BlindCaseResult {
            id: case.id.clone(),
            reference: case.reference,
            predicted: case.predicted,
            absolute_relative_error_percent: error,
            pass: error <= dataset.max_relative_error_percent,
            unit: case.unit.clone(),
        });
    }
    let passed = results.iter().filter(|r| r.pass).count();
    let sum: f64 = results
        .iter()
        .map(|r| r.absolute_relative_error_percent)
        .sum();
    let max = results
        .iter()
        .map(|r| r.absolute_relative_error_percent)
        .fold(0.0_f64, f64::max);
    Ok(BlindReport {
        dataset: dataset.dataset.clone(),
        case_count: results.len(),
        passed,
        mean_absolute_relative_error_percent: sum / results.len() as f64,
        max_absolute_relative_error_percent: max,
        tolerance_percent: dataset.max_relative_error_percent,
        overall_pass: passed == results.len(),
        results,
    })
}

/// Read a JSON dataset, evaluate it, and save JSON and HTML reports.
pub fn evaluate_file(input: &Path, output_dir: &Path) -> Result<BlindReport, BlindError> {
    let dataset: BlindDataset = serde_json::from_slice(&fs::read(input)?)?;
    let report = evaluate(&dataset)?;
    fs::create_dir_all(output_dir)?;
    fs::write(
        output_dir.join("blind_report.json"),
        serde_json::to_vec_pretty(&report)?,
    )?;
    fs::write(output_dir.join("blind_report.html"), to_html(&report))?;
    Ok(report)
}

/// Render a simple standalone HTML report.
pub fn to_html(report: &BlindReport) -> String {
    let status = if report.overall_pass { "PASS" } else { "FAIL" };
    let mut rows = String::new();
    for case in &report.results {
        rows.push_str(&format!(
            "<tr><td>{}</td><td>{:.6}</td><td>{:.6}</td><td>{:.4}%</td><td>{}</td></tr>",
            escape_html(&case.id),
            case.reference,
            case.predicted,
            case.absolute_relative_error_percent,
            if case.pass { "PASS" } else { "FAIL" }
        ));
    }
    format!("<!doctype html><html lang=\"en\"><meta charset=\"utf-8\"><title>AttoFab Blind Validation</title><body><h1>AttoFab Blind Validation</h1><p>Dataset: {}</p><h2>{}</h2><p>Cases passed: {}/{}; mean absolute relative error: {:.4}%; maximum: {:.4}%; tolerance: {:.4}%</p><table border=\"1\" cellpadding=\"6\"><thead><tr><th>Case</th><th>Reference</th><th>Predicted</th><th>Abs. relative error</th><th>Status</th></tr></thead><tbody>{}</tbody></table><p><em>Interpretation is limited to the supplied input records; report generation does not certify calibration.</em></p></body></html>",
        escape_html(&report.dataset), status, report.passed, report.case_count,
        report.mean_absolute_relative_error_percent, report.max_absolute_relative_error_percent,
        report.tolerance_percent, rows)
}

fn escape_html(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&#39;")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evaluates_tolerance_and_html_escapes_names() {
        let d = BlindDataset {
            dataset: "demo <set>".into(),
            max_relative_error_percent: 5.0,
            cases: vec![
                BlindCase {
                    id: "w1".into(),
                    reference: 100.0,
                    predicted: 103.0,
                    unit: "nm".into(),
                },
                BlindCase {
                    id: "w2".into(),
                    reference: 100.0,
                    predicted: 108.0,
                    unit: "nm".into(),
                },
            ],
        };
        let r = evaluate(&d).unwrap();
        assert_eq!(r.passed, 1);
        assert!(!r.overall_pass);
        assert!(to_html(&r).contains("demo &lt;set&gt;"));
    }
}
