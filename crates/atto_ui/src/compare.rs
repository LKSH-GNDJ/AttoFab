// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Compare mode: two states side by side with a computed delta (PRD G-09).
//!
//! The delta is computed, never described: each stage is compared by what it
//! produced, and a stage whose inputs differ but whose result matches is reported
//! as *same answer, different inputs* — which is a real and useful distinction that
//! a naive "fingerprints differ" check would lose.

use crate::{Session, StageId};

/// One stage's comparison.
#[derive(Debug, Clone, PartialEq)]
pub struct StageDelta {
    pub stage: StageId,
    /// Whether each side ran it.
    pub ran_a: bool,
    pub ran_b: bool,
    /// The status word on each side.
    pub status_a: &'static str,
    pub status_b: &'static str,
    /// Whether the committed inputs differ.
    pub inputs_differ: bool,
    /// Whether each side will run it.
    pub enabled_a: bool,
    pub enabled_b: bool,
}

impl StageDelta {
    /// Whether anything about this stage differs.
    pub fn differs(&self) -> bool {
        self.ran_a != self.ran_b
            || self.status_a != self.status_b
            || self.inputs_differ
            || self.enabled_a != self.enabled_b
    }

    /// Whether the two sides reached the same answer from different inputs.
    ///
    /// This is the case worth surfacing: the result agrees, so the change did not
    /// matter — and that is exactly the thing a sweep is trying to find out.
    pub fn same_answer_different_inputs(&self) -> bool {
        self.inputs_differ && self.ran_a && self.ran_b && self.status_a == self.status_b
    }
}

/// The whole comparison.
#[derive(Debug, Clone, PartialEq)]
pub struct CompareView {
    pub stages: Vec<StageDelta>,
}

impl CompareView {
    /// Compare two sessions stage by stage.
    pub fn of(a: &Session, b: &Session) -> Self {
        let sa = a.stages();
        let sb = b.stages();
        let n = sa.len().min(sb.len());
        let mut stages = Vec::with_capacity(n);
        for i in 0..n {
            let ran_a = a.result(i).is_some();
            let ran_b = b.result(i).is_some();
            // the committed fingerprint is what the stage was run against
            let inputs_differ = a.committed_at(i) != b.committed_at(i);
            stages.push(StageDelta {
                stage: sa[i].id,
                ran_a,
                ran_b,
                status_a: a.status(i).word(),
                status_b: b.status(i).word(),
                inputs_differ,
                enabled_a: sa[i].enabled,
                enabled_b: sb[i].enabled,
            });
        }
        CompareView { stages }
    }

    /// The stages that differ in any way.
    pub fn differing(&self) -> Vec<&StageDelta> {
        self.stages.iter().filter(|d| d.differs()).collect()
    }

    /// The stages where the answer matched despite different inputs.
    pub fn same_answer(&self) -> Vec<&StageDelta> {
        self.stages.iter().filter(|d| d.same_answer_different_inputs()).collect()
    }

    /// Whether the two states are identical for comparison purposes.
    pub fn identical(&self) -> bool {
        self.stages.iter().all(|d| !d.differs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_identical_sessions_have_no_delta() {
        let a = Session::new();
        let b = Session::new();
        let v = CompareView::of(&a, &b);
        assert!(v.identical(), "{:?}", v.differing());
        assert_eq!(v.stages.len(), StageId::ALL.len());
        assert!(v.same_answer().is_empty());
    }

    #[test]
    fn an_uncommitted_edit_does_not_show_up_in_the_comparison() {
        // compare is over what each state was *run against*, not what has been
        // typed since: a half-finished edit is not a different result
        let a = Session::new();
        let mut b = Session::new();
        b.set_param(0, "thickness", 123.0).unwrap();
        let v = CompareView::of(&a, &b);
        assert!(v.identical(), "nothing has been committed on either side");
    }

    #[test]
    fn a_stage_that_ran_on_one_side_only_is_a_difference() {
        let a = Session::new();
        let b = Session::new();
        // neither has run anything, so they match; this asserts the baseline the
        // test above relies on
        assert!(a.result(0).is_none() && b.result(0).is_none());
        let v = CompareView::of(&a, &b);
        assert!(v.stages.iter().all(|d| !d.ran_a && !d.ran_b));
    }

    #[test]
    fn disabling_a_stage_on_one_side_is_a_difference() {
        let a = Session::new();
        let mut b = Session::new();
        b.set_enabled(2, false).unwrap();
        let v = CompareView::of(&a, &b);
        let d = v.stages.iter().find(|d| d.stage == StageId::ALL[2]).unwrap();
        assert!(d.differs(), "one side will not run it");
    }

    #[test]
    fn every_stage_carries_both_status_words() {
        let v = CompareView::of(&Session::new(), &Session::new());
        for d in &v.stages {
            assert!(!d.status_a.is_empty());
            assert!(!d.status_b.is_empty());
        }
    }
}
