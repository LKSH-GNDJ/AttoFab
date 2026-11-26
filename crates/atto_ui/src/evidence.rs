// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The evidence tab: verification cases and blind datasets (PRD G-15).
//!
//! Synthetic data is shown but never counted: a row that reads like validation
//! evidence while being fabricated is worse than no row.

use atto_evidence::blind::BlindDataset;
use atto_evidence::exact::{cases, run, Verdict};

use crate::Colour;

/// One verification case, ready to draw.
#[derive(Debug, Clone, PartialEq)]
pub struct CaseRow {
    pub name: String,
    /// Live, not live, wrong or malformed.
    pub verdict: String,
    pub colour: Colour,
    /// The mutation that must break it.
    pub mutation: String,
}

/// One blind dataset, ready to draw.
#[derive(Debug, Clone, PartialEq)]
pub struct DatasetRow {
    pub name: String,
    pub tolerance_pct: f64,
    pub cases: usize,
    /// Whether the name marks it as fabricated.
    pub synthetic: bool,
    /// What a report may claim about it.
    pub counts_as_evidence: bool,
    pub colour: Colour,
}

/// Build the dataset rows, and count measured against synthetic.
///
/// A synthetic set is shown but never counted: a row that reads like validation
/// evidence while being fabricated is worse than no row.
pub fn dataset_rows(datasets: &[BlindDataset]) -> (Vec<DatasetRow>, usize, usize) {
    let mut rows = Vec::with_capacity(datasets.len());
    let mut measured = 0;
    let mut synthetic = 0;
    for d in datasets {
        if d.synthetic {
            synthetic += 1;
        } else {
            measured += 1;
        }
        rows.push(DatasetRow {
            name: d.name.clone(),
            tolerance_pct: d.tolerance_pct,
            cases: d.cases.len(),
            synthetic: d.synthetic,
            counts_as_evidence: !d.synthetic,
            colour: if d.synthetic { Colour::Grey } else { Colour::Green },
        });
    }
    (rows, measured, synthetic)
}

/// The tab's summary, which must not be inflated by synthetic data.
#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceSummary {
    pub cases_live: usize,
    pub cases_not_live: usize,
    /// Datasets that may be counted as validation evidence.
    pub measured_datasets: usize,
    /// Datasets shown but excluded from the count.
    pub synthetic_datasets: usize,
}

impl EvidenceSummary {
    /// Whether anything in the tab is evidence at all.
    pub fn has_evidence(&self) -> bool {
        self.measured_datasets > 0 && self.cases_not_live == 0
    }

    /// The line a report uses.
    pub fn sentence(&self) -> String {
        format!(
            "{} live cases, {} not live; {} measured datasets, {} synthetic (not counted)",
            self.cases_live, self.cases_not_live, self.measured_datasets, self.synthetic_datasets
        )
    }
}

/// The case rows, with the verdict each case actually earned.
pub fn case_rows<F>(evaluate: F) -> (Vec<CaseRow>, usize, usize)
where
    F: Fn(&str, f64) -> f64,
{
    let mut rows = Vec::new();
    let mut live = 0;
    let mut not_live = 0;
    for case in cases() {
        let verdict = run(&case, |f| evaluate(case.name, f));
        let (word, colour) = match &verdict {
            Verdict::Live => {
                live += 1;
                ("live".to_string(), Colour::Green)
            }
            Verdict::NotLive => {
                not_live += 1;
                ("not live".to_string(), Colour::Amber)
            }
            Verdict::Wrong { got, want } => {
                not_live += 1;
                (format!("wrong: {got:.4} against {want:.4}"), Colour::Red)
            }
            Verdict::Malformed(why) => {
                not_live += 1;
                (format!("malformed: {why}"), Colour::Red)
            }
        };
        rows.push(CaseRow {
            name: case.name.to_string(),
            verdict: word,
            colour,
            mutation: format!("{} x{:.2}", case.mutation.parameter, case.mutation.factor),
        });
    }
    (rows, live, not_live)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn live_model(name: &str, f: f64) -> f64 {
        let all = cases();
        let case = all.iter().find(|c| c.name == name).expect("a known case");
        // a fourth power moves more than any case tolerance, so the model is
        // genuinely sensitive to the mutation it is asked about
        case.reference * f.powi(4)
    }

    #[test]
    fn a_live_model_marks_every_case_live() {
        let (rows, live, not_live) = case_rows(live_model);
        assert_eq!(rows.len(), cases().len());
        assert_eq!(not_live, 0);
        assert_eq!(live, rows.len());
        assert!(rows.iter().all(|r| r.colour == Colour::Green));
    }

    #[test]
    fn a_constant_model_is_reported_not_live() {
        let first = cases()[0].reference;
        let (rows, live, not_live) = case_rows(|_, _| first);
        assert!(not_live > 0, "a constant model cannot be live");
        assert_eq!(live + not_live, rows.len());
        assert!(rows.iter().any(|r| r.verdict == "not live"));
    }

    #[test]
    fn every_row_carries_its_mutation() {
        let (rows, _, _) = case_rows(live_model);
        for row in &rows {
            assert!(row.mutation.contains('x'), "the mutation is shown: {}", row.mutation);
            assert!(!row.name.is_empty());
        }
    }

    fn dataset(name: &str, synthetic: bool) -> BlindDataset {
        BlindDataset {
            name: name.to_string(),
            tolerance_pct: 5.0,
            cases: vec![atto_evidence::blind::BlindCase {
                id: "CD".into(),
                reference: 96.5,
                predicted: 96.0,
                unit: "nm".into(),
            }],
            synthetic,
            generated_from: Some("run 1".into()),
        }
    }

    #[test]
    fn synthetic_datasets_are_shown_but_never_counted() {
        // the rule the tab exists to keep
        let sets = [dataset("FAB-A lot 1", false), dataset("synthetic-1", true)];
        let (rows, measured, synthetic) = dataset_rows(&sets);
        assert_eq!(rows.len(), 2, "both are shown");
        assert_eq!(measured, 1);
        assert_eq!(synthetic, 1);
        let syn = rows.iter().find(|r| r.synthetic).unwrap();
        assert!(!syn.counts_as_evidence, "a synthetic set is not evidence");
        assert_eq!(syn.colour, Colour::Grey, "and it is drawn as such");
        let real = rows.iter().find(|r| !r.synthetic).unwrap();
        assert!(real.counts_as_evidence);
        assert_eq!(real.colour, Colour::Green);
    }

    #[test]
    fn a_dataset_row_carries_its_tolerance_and_case_count() {
        let (rows, _, _) = dataset_rows(&[dataset("FAB-B", false)]);
        assert_eq!(rows[0].tolerance_pct, 5.0);
        assert_eq!(rows[0].cases, 1);
        assert!(!rows[0].name.is_empty());
    }

    #[test]
    fn the_summary_refuses_to_count_synthetic_data() {
        let s = EvidenceSummary {
            cases_live: 3,
            cases_not_live: 0,
            measured_datasets: 0,
            synthetic_datasets: 4,
        };
        assert!(!s.has_evidence(), "four synthetic datasets are not evidence");
        assert!(s.sentence().contains("not counted"));
    }

    #[test]
    fn measured_data_with_live_cases_is_evidence() {
        let s = EvidenceSummary {
            cases_live: 3,
            cases_not_live: 0,
            measured_datasets: 1,
            synthetic_datasets: 2,
        };
        assert!(s.has_evidence(), "{}", s.sentence());
    }
}
