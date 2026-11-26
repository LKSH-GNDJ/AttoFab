// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The inspector's five tabs (PRD G-11).
//!
//! Each tab keeps its own selection and scroll. Switching tabs and coming back
//! must land where you left, or a user loses their place every time they compare
//! a measurement against the stack — which is the whole reason the tabs exist.

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The five things the inspector can show.
pub enum InspectorTab {
    Stage,
    Stack,
    Measure,
    Device,
    Evidence,
}

impl InspectorTab {
    /// Every tab, in display order.
    pub const ALL: [InspectorTab; 5] = [
        InspectorTab::Stage,
        InspectorTab::Stack,
        InspectorTab::Measure,
        InspectorTab::Device,
        InspectorTab::Evidence,
    ];

    /// The tab name.
    pub fn label(self) -> &'static str {
        match self {
            InspectorTab::Stage => "Stage",
            InspectorTab::Stack => "Stack",
            InspectorTab::Measure => "Measure",
            InspectorTab::Device => "Device",
            InspectorTab::Evidence => "Evidence",
        }
    }

    fn index(self) -> usize {
        InspectorTab::ALL.iter().position(|t| *t == self).unwrap_or(0)
    }
}

/// What a single tab remembers between visits.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct TabState {
    /// Selected row, if the tab has a list.
    pub selection: Option<usize>,
    /// Scroll offset in rows.
    pub scroll: f64,
}

#[derive(Debug, Clone, PartialEq)]
/// Per-tab selection and scroll, so switching tabs does not lose place.
pub struct Inspector {
    current: InspectorTab,
    state: [TabState; 5],
}

impl Default for Inspector {
    fn default() -> Self {
        Inspector::new()
    }
}

impl Inspector {
    /// A fresh inspector on the first tab.
    pub fn new() -> Self {
        Inspector { current: InspectorTab::Stage, state: [TabState::default(); 5] }
    }

    /// The tab being shown.
    pub fn current(&self) -> InspectorTab {
        self.current
    }

    /// Show a tab, keeping each tab's own selection and scroll.
    pub fn select_tab(&mut self, tab: InspectorTab) {
        self.current = tab;
    }

    /// The state of the tab in view.
    pub fn state(&self) -> &TabState {
        &self.state[self.current.index()]
    }

    /// The state of any tab, for restoring a layout.
    pub fn state_of(&self, tab: InspectorTab) -> &TabState {
        &self.state[tab.index()]
    }

    /// Highlight a row, or clear the selection.
    pub fn select_row(&mut self, row: Option<usize>) {
        self.state[self.current.index()].selection = row;
    }

    /// Remember a scroll offset for the current tab.
    pub fn scroll_to(&mut self, rows: f64) {
        self.state[self.current.index()].scroll = rows.max(0.0);
    }

    /// The tab bar, with the current one marked.
    pub fn tabs(&self) -> Vec<(&'static str, bool)> {
        InspectorTab::ALL.iter().map(|t| (t.label(), *t == self.current)).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn each_tab_keeps_its_own_selection_and_scroll() {
        let mut i = Inspector::new();
        i.select_tab(InspectorTab::Measure);
        i.select_row(Some(3));
        i.scroll_to(12.0);

        i.select_tab(InspectorTab::Stack);
        assert_eq!(i.state().selection, None, "a fresh tab starts clean");
        i.select_row(Some(1));

        i.select_tab(InspectorTab::Measure);
        assert_eq!(i.state().selection, Some(3), "the measure tab remembered its row");
        assert!((i.state().scroll - 12.0).abs() < 1e-12, "and its scroll");

        assert_eq!(i.state_of(InspectorTab::Stack).selection, Some(1));
    }

    #[test]
    fn scroll_cannot_go_negative() {
        let mut i = Inspector::new();
        i.scroll_to(-5.0);
        assert!((i.state().scroll - 0.0).abs() < 1e-12);
    }

    #[test]
    fn the_tab_bar_marks_exactly_one_current_tab() {
        let mut i = Inspector::new();
        i.select_tab(InspectorTab::Device);
        let tabs = i.tabs();
        assert_eq!(tabs.len(), 5);
        assert_eq!(tabs.iter().filter(|(_, cur)| *cur).count(), 1);
        assert!(tabs.iter().any(|(l, cur)| *l == "Device" && *cur));
        assert_eq!(i.current(), InspectorTab::Device);
    }

    #[test]
    fn every_tab_has_a_distinct_label() {
        let mut labels: Vec<&str> = InspectorTab::ALL.iter().map(|t| t.label()).collect();
        labels.sort();
        let n = labels.len();
        labels.dedup();
        assert_eq!(labels.len(), n, "two tabs share a label");
    }
}
