// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Provenance for whatever the operator is pointing at (PRD G-03).
//!
//! The requirement asks for two routes to the same thing: the pointer and the
//! keyboard. That is only honest if both go through one lookup, so this module is
//! that lookup — a `Focused` value that says what is under the pointer or the
//! caret, and a single `for_focus` that answers for it.

use crate::inspector::InspectorTab;
use crate::Provenance;

/// What the operator is currently on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Focused {
    /// A stage's recipe row.
    RecipeRow(usize),
    /// A measurement row.
    Measurement(usize),
    /// A stack row.
    Stack(usize),
    /// A case in the evidence tab.
    EvidenceCase(usize),
    /// Nothing is focused.
    None,
}

/// One tab's worth of provenance, so the lookup can answer for any of them.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ProvenanceIndex {
    recipe: Vec<Provenance>,
    measurements: Vec<Provenance>,
    stack: Vec<Provenance>,
    evidence: Vec<Provenance>,
}

impl ProvenanceIndex {
    /// An empty index.
    pub fn new() -> Self {
        ProvenanceIndex::default()
    }

    /// Record provenance for a row in a tab.
    pub fn set(&mut self, tab: InspectorTab, row: usize, p: Provenance) {
        let slot = match tab {
            InspectorTab::Stage => &mut self.recipe,
            InspectorTab::Measure => &mut self.measurements,
            InspectorTab::Stack => &mut self.stack,
            InspectorTab::Evidence => &mut self.evidence,
            _ => return,
        };
        if slot.len() <= row {
            // a gap is filled with an empty record, which `complete` reports as a
            // hole rather than as provenance
            slot.resize(row + 1, Provenance {
                model: String::new(),
                parameters: Vec::new(),
                inputs_fingerprint: 0,
                uncertainty: None,
                note: String::new(),
            });
        }
        slot[row] = p;
    }

    /// The provenance for whatever is focused.
    ///
    /// Both routes call this: a pointer hovering a row and a keyboard action on the
    /// same row must produce the same answer, and the only way to guarantee that is
    /// to have one function.
    pub fn for_focus(&self, focus: Focused) -> Option<&Provenance> {
        match focus {
            Focused::RecipeRow(i) => self.recipe.get(i),
            Focused::Measurement(i) => self.measurements.get(i),
            Focused::Stack(i) => self.stack.get(i),
            Focused::EvidenceCase(i) => self.evidence.get(i),
            Focused::None => None,
        }
    }

    /// Whether every row in a tab carries provenance.
    ///
    /// A row without one is a derived number that cannot be explained, which is
    /// what G-03 forbids.
    pub fn complete(&self, tab: InspectorTab) -> bool {
        let slot = match tab {
            InspectorTab::Stage => &self.recipe,
            InspectorTab::Measure => &self.measurements,
            InspectorTab::Stack => &self.stack,
            InspectorTab::Evidence => &self.evidence,
            _ => return true,
        };
        slot.iter().all(|p| !p.model.is_empty())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn prov(model: &str) -> Provenance {
        Provenance {
            model: model.to_string(),
            parameters: vec!["thickness".into()],
            inputs_fingerprint: 7,
            uncertainty: Some(1.0),
            note: String::new(),
        }
    }

    #[test]
    fn the_pointer_and_the_keyboard_get_the_same_answer() {
        // both routes are the same call, so this is what makes the requirement
        // checkable rather than aspirational
        let mut idx = ProvenanceIndex::new();
        idx.set(InspectorTab::Measure, 2, prov("atto_measure::optics"));
        let focus = Focused::Measurement(2);
        let by_pointer = idx.for_focus(focus);
        let by_keyboard = idx.for_focus(focus);
        assert_eq!(by_pointer, by_keyboard);
        assert_eq!(by_pointer.map(|p| p.model.as_str()), Some("atto_measure::optics"));
    }

    #[test]
    fn nothing_focused_has_no_provenance() {
        let idx = ProvenanceIndex::new();
        assert!(idx.for_focus(Focused::None).is_none());
    }

    #[test]
    fn a_gap_left_by_a_later_row_is_reported_as_incomplete() {
        // setting row 2 grows the list to three, and the two gaps are holes the
        // completeness check exists to find
        let mut idx = ProvenanceIndex::new();
        idx.set(InspectorTab::Stage, 2, prov("atto_process::oxidation"));
        assert!(!idx.complete(InspectorTab::Stage), "rows 0 and 1 are holes");
        assert_eq!(idx.for_focus(Focused::RecipeRow(0)).map(|p| p.model.as_str()), Some(""));
        assert!(idx.for_focus(Focused::RecipeRow(5)).is_none(), "beyond the end");
        // fill the gaps and the tab becomes complete
        idx.set(InspectorTab::Stage, 0, prov("atto_process::prep"));
        idx.set(InspectorTab::Stage, 1, prov("atto_process::deposition"));
        assert!(idx.complete(InspectorTab::Stage));
    }

    #[test]
    fn every_tab_is_addressed_separately() {
        let mut idx = ProvenanceIndex::new();
        idx.set(InspectorTab::Stack, 0, prov("atto_geom::column"));
        assert!(idx.for_focus(Focused::Stack(0)).is_some());
        assert!(idx.for_focus(Focused::Measurement(0)).is_none(), "a different tab");
    }
}
