// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Preview values that cannot overwrite committed ones (PRD G-05).
//!
//! A preview is a promise that nothing has changed yet. The way to keep that
//! promise is to make the committed values unreachable from the preview: this type
//! holds its own map keyed by [`Knob`], and the knob set is only written by an
//! explicit `commit`.

use atto_study::{Knob, Knobs};

#[derive(Debug, Clone, Default, PartialEq)]
/// Tentative knob values that have not been committed.
pub struct Preview {
    /// Knobs the operator has changed but not committed.
    pending: Vec<(Knob, f64)>,
}

impl Preview {
    /// An empty preview: nothing pending, nothing committed.
    pub fn new() -> Self {
        Preview { pending: Vec::new() }
    }

    /// Note a tentative value. Does not touch the knob set.
    pub fn set(&mut self, knob: Knob, value: f64) {
        match self.pending.iter_mut().find(|(k, _)| *k == knob) {
            Some(slot) => slot.1 = value,
            None => self.pending.push((knob, value)),
        }
    }

    /// The tentative value for a knob, if it has one.
    pub fn get(&self, knob: Knob) -> Option<f64> {
        self.pending.iter().find(|(k, _)| *k == knob).map(|(_, v)| *v)
    }

    /// How many knobs have a pending value.
    pub fn len(&self) -> usize {
        self.pending.len()
    }

    /// Whether anything is pending.
    pub fn is_empty(&self) -> bool {
        self.pending.is_empty()
    }

    /// What a panel draws for a knob: the preview value where there is one, the
    /// committed value otherwise, and a flag saying which — the visual distinction
    /// the requirement asks for.
    pub fn shown(&self, knobs: &Knobs, knob: Knob) -> (f64, bool) {
        match self.get(knob) {
            Some(v) => (v, true),
            None => (knob.get(knobs), false),
        }
    }

    /// Apply every pending value and clear the preview.
    ///
    /// Returns how many knobs changed. This is the only path from a preview to a
    /// committed value, which is what makes the promise checkable.
    pub fn commit(&mut self, knobs: &mut Knobs) -> usize {
        let n = self.pending.len();
        for (k, v) in self.pending.drain(..) {
            k.set(knobs, v);
        }
        n
    }

    /// Throw the pending values away. The knob set is untouched either way.
    pub fn discard(&mut self) {
        self.pending.clear();
    }
}

/// Whether the bottom dock should collapse (PRD G-02).
///
/// It collapses only when there is nothing to show: no running stage and an empty
/// console. A busy timeline or a console with a line in it keeps the dock open.
pub fn dock_collapsed(timeline_idle: bool, console_empty: bool) -> bool {
    timeline_idle && console_empty
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_preview_never_touches_the_committed_value() {
        let mut knobs = Knobs::default();
        let before = Knob::ImplantEnergy.get(&knobs);
        let mut p = Preview::new();
        p.set(Knob::ImplantEnergy, before + 50.0);
        assert_eq!(Knob::ImplantEnergy.get(&knobs), before, "the knob set moved");
        assert_eq!(p.get(Knob::ImplantEnergy), Some(before + 50.0));
    }

    #[test]
    fn the_shown_value_says_whether_it_is_previewed() {
        let knobs = Knobs::default();
        let mut p = Preview::new();
        let committed = Knob::AnnealTemp.get(&knobs);
        let (v, is_preview) = p.shown(&knobs, Knob::AnnealTemp);
        assert_eq!(v, committed);
        assert!(!is_preview, "nothing pending, so it is committed");

        p.set(Knob::AnnealTemp, committed + 25.0);
        let (v2, is_preview2) = p.shown(&knobs, Knob::AnnealTemp);
        assert_eq!(v2, committed + 25.0);
        assert!(is_preview2, "a pending value must be marked as one");
    }

    #[test]
    fn committing_applies_every_pending_value() {
        let mut knobs = Knobs::default();
        let mut p = Preview::new();
        p.set(Knob::AnnealTemp, 950.0);
        p.set(Knob::OxMinutes, 42.0);
        assert_eq!(p.commit(&mut knobs), 2);
        assert_eq!(Knob::AnnealTemp.get(&knobs), 950.0);
        assert_eq!(Knob::OxMinutes.get(&knobs), 42.0);
        assert!(p.is_empty(), "the preview is consumed by the commit");
    }

    #[test]
    fn discarding_leaves_the_knobs_exactly_as_they_were() {
        let knobs = Knobs::default();
        let before = knobs.clone();
        let mut p = Preview::new();
        p.set(Knob::AnnealTemp, 1234.0);
        p.discard();
        assert_eq!(knobs, before, "discard changed the knobs");
        assert!(p.is_empty());
    }

    #[test]
    fn the_dock_collapses_only_when_there_is_nothing_to_show() {
        assert!(dock_collapsed(true, true));
        assert!(!dock_collapsed(false, true), "a running stage keeps it open");
        assert!(!dock_collapsed(true, false), "a console line keeps it open");
        assert!(!dock_collapsed(false, false));
    }
}
