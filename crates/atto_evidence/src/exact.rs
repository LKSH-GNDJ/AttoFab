// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Exact-solution cases and the mutation that proves they are live (PRD R-07).
//!
//! A verification case that passes whatever the code does is decoration. Each
//! case here declares a **named mutation** — a parameter perturbed by a stated
//! factor — and the runner requires two things: the unmutated result matches the
//! reference, and the mutated one does not. A case whose mutation changes nothing
//! is reported as not live rather than passed.

/// The perturbation that must break a case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Mutation {
    /// Which parameter is moved, by name.
    pub parameter: &'static str,
    /// The factor it is multiplied by. Must not be 1.
    pub factor: f64,
}

impl Mutation {
    fn valid(&self) -> bool {
        self.factor.is_finite() && (self.factor - 1.0).abs() > 1e-9
    }
}

/// One exact-solution case.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExactCase {
    pub name: &'static str,
    /// What the model must produce.
    pub reference: f64,
    /// How close it must be.
    pub tolerance: f64,
    /// The perturbation that must move the answer outside the tolerance.
    pub mutation: Mutation,
}

#[derive(Debug, Clone, PartialEq)]
/// What running a case established.
pub enum Verdict {
    /// Matches the reference and the mutation breaks it.
    Live,
    /// Matches the reference, but the mutation leaves it inside the tolerance:
    /// the case cannot detect the change it claims to detect.
    NotLive,
    /// Does not match the reference.
    Wrong { got: f64, want: f64 },
    /// The case or its mutation is malformed.
    Malformed(&'static str),
}

/// Run a case against a model. `evaluate` takes the mutation factor applied to the
/// named parameter and returns the model's answer.
pub fn run<F>(case: &ExactCase, evaluate: F) -> Verdict
where
    F: Fn(f64) -> f64,
{
    if !case.tolerance.is_finite() || case.tolerance <= 0.0 || !case.reference.is_finite() {
        return Verdict::Malformed("the case needs a finite reference and a positive tolerance");
    }
    if !case.mutation.valid() {
        return Verdict::Malformed("the mutation factor must differ from one");
    }
    let plain = evaluate(1.0);
    if !plain.is_finite() || (plain - case.reference).abs() > case.tolerance {
        return Verdict::Wrong { got: plain, want: case.reference };
    }
    let mutated = evaluate(case.mutation.factor);
    if !mutated.is_finite() || (mutated - case.reference).abs() <= case.tolerance {
        return Verdict::NotLive;
    }
    Verdict::Live
}

/// The bundled cases: one per model that has a closed-form answer.
pub fn cases() -> Vec<ExactCase> {
    vec![
        ExactCase {
            name: "deal-grove dry 1000 C, 40 min",
            reference: 37.6,
            tolerance: 2.0,
            mutation: Mutation { parameter: "rate_b", factor: 1.10 },
        },
        ExactCase {
            name: "silicon consumed by the oxide",
            reference: 16.5,
            tolerance: 1.0,
            mutation: Mutation { parameter: "consumption_ratio", factor: 1.05 },
        },
        ExactCase {
            // R-08: the resist CD has an independent route. With no absorption
            // every exposed point bleaches equally, so the half-developed
            // threshold lands exactly on the drawn edge.
            name: "resist CD with no absorption equals the drawn width",
            reference: 200.0,
            tolerance: 1.0,
            mutation: Mutation { parameter: "development_threshold", factor: 1.20 },
        },
        ExactCase {
            name: "implant range at 60 keV",
            reference: 48.0,
            tolerance: 5.0,
            mutation: Mutation { parameter: "energy_kev", factor: 1.10 },
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_live_case_passes_both_requirements() {
        // the model moves with the mutated parameter, as it should
        let c = ExactCase {
            name: "linear",
            reference: 100.0,
            tolerance: 1.0,
            mutation: Mutation { parameter: "slope", factor: 1.10 },
        };
        assert_eq!(run(&c, |f| 100.0 * f), Verdict::Live);
    }

    #[test]
    fn a_case_whose_mutation_changes_nothing_is_not_live() {
        // this is the finding R-07 exists for: a check that cannot fail
        let c = ExactCase {
            name: "insensitive",
            reference: 100.0,
            tolerance: 1.0,
            mutation: Mutation { parameter: "ignored", factor: 2.0 },
        };
        assert_eq!(run(&c, |_| 100.0), Verdict::NotLive);
    }

    #[test]
    fn a_case_that_misses_its_reference_is_wrong_not_live() {
        let c = ExactCase {
            name: "off by ten",
            reference: 100.0,
            tolerance: 1.0,
            mutation: Mutation { parameter: "slope", factor: 1.10 },
        };
        match run(&c, |f| 110.0 * f) {
            Verdict::Wrong { got, want } => {
                assert!((got - 110.0).abs() < 1e-9);
                assert!((want - 100.0).abs() < 1e-9);
            }
            other => panic!("{other:?}"),
        }
    }

    #[test]
    fn a_malformed_case_is_reported_rather_than_run() {
        let mut c = cases()[0];
        c.mutation.factor = 1.0;
        assert!(matches!(run(&c, |_| c.reference), Verdict::Malformed(_)));
        let mut d = cases()[0];
        d.tolerance = 0.0;
        assert!(matches!(run(&d, |_| d.reference), Verdict::Malformed(_)));
    }

    #[test]
    fn every_bundled_case_declares_a_distinct_mutation() {
        let all = cases();
        assert!(!all.is_empty());
        for c in &all {
            assert!(c.mutation.valid(), "{} has an invalid mutation", c.name);
            assert!(!c.mutation.parameter.is_empty(), "{} names no parameter", c.name);
        }
        let mut names: Vec<&str> = all.iter().map(|c| c.name).collect();
        names.sort();
        let n = names.len();
        names.dedup();
        assert_eq!(names.len(), n, "two cases share a name");
    }
}
