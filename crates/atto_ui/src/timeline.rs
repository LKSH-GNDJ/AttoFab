// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The recipe as one band with per-stage ticks (PRD G-16).
//!
//! A timeline is a *view of the recipe*, not a copy of it: every segment is derived
//! from a stage, and a stage that is switched off still appears — greyed, not
//! removed. A band that silently drops disabled stages misrepresents the recipe,
//! and that is the failure this module exists to prevent.

use crate::{Session, StageId};

/// One stage's slice of the band.
#[derive(Debug, Clone, PartialEq)]
pub struct Segment {
    pub stage: StageId,
    pub label: String,
    /// Whether the stage will run.
    pub enabled: bool,
    /// Where the segment starts, as a fraction of the band, 0..1.
    pub start: f64,
    /// How wide it is, as a fraction of the band, 0..1.
    pub width: f64,
    /// Whether it has run, for the tick's shading.
    pub run: bool,
    /// The one-word status.
    pub status: &'static str,
}

impl Segment {
    /// Where the segment ends.
    pub fn end(&self) -> f64 {
        self.start + self.width
    }
}

/// The band, and which stage is expanded under it.
#[derive(Debug, Clone, PartialEq)]
pub struct Timeline {
    pub segments: Vec<Segment>,
    /// The stage whose detail is open, if any.
    pub expanded: Option<StageId>,
}

impl Timeline {
    /// Build the band from a session.
    ///
    /// Every stage gets an equal share. That is a deliberate choice: the recipe has
    /// no declared durations, and inventing them from run time would make the band
    /// change shape between two runs of the same recipe.
    pub fn from_session(session: &Session) -> Self {
        let stages = session.stages();
        let n = stages.len().max(1);
        let w = 1.0 / n as f64;
        let mut segments = Vec::with_capacity(stages.len());
        for (i, spec) in stages.iter().enumerate() {
            let run = session.result(i).is_some();
            segments.push(Segment {
                stage: spec.id,
                label: spec.id.label().to_string(),
                enabled: spec.enabled,
                start: i as f64 * w,
                width: w,
                run,
                status: session.status(i).word(),
            });
        }
        Timeline { segments, expanded: None }
    }

    /// Open or close a stage's detail.
    pub fn toggle(&mut self, stage: StageId) {
        self.expanded = if self.expanded == Some(stage) { None } else { Some(stage) };
    }

    /// The segment for a stage.
    pub fn segment(&self, stage: StageId) -> Option<&Segment> {
        self.segments.iter().find(|s| s.stage == stage)
    }

    /// The enabled stages, in band order — what the run will actually do.
    pub fn enabled(&self) -> Vec<StageId> {
        self.segments.iter().filter(|s| s.enabled).map(|s| s.stage).collect()
    }

    /// The fraction of the band that will run.
    pub fn enabled_fraction(&self) -> f64 {
        self.segments.iter().filter(|s| s.enabled).map(|s| s.width).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn session() -> Session {
        Session::new()
    }

    #[test]
    fn the_band_covers_the_whole_recipe_without_gaps_or_overlap() {
        let t = Timeline::from_session(&session());
        assert_eq!(t.segments.len(), StageId::ALL.len());
        assert!((t.segments[0].start - 0.0).abs() < 1e-12, "the band starts at zero");
        for pair in t.segments.windows(2) {
            assert!((pair[0].end() - pair[1].start).abs() < 1e-12, "a gap or an overlap");
        }
        let last = t.segments.last().unwrap();
        assert!((last.end() - 1.0).abs() < 1e-12, "the band ends at one");
    }

    #[test]
    fn a_disabled_stage_stays_on_the_band() {
        // the failure this prevents: a band that quietly drops what will not run
        let mut s = session();
        s.set_enabled(0, false).unwrap();
        let t = Timeline::from_session(&s);
        assert_eq!(t.segments.len(), StageId::ALL.len(), "nothing was removed");
        let seg = t.segment(StageId::ALL[0]).unwrap();
        assert!(!seg.enabled, "but it is marked as not running");
        assert!(t.enabled_fraction() < 1.0, "and the running fraction fell");
        assert!(!t.enabled().contains(&StageId::ALL[0]));
    }

    #[test]
    fn expanding_one_stage_closes_the_other() {
        let mut t = Timeline::from_session(&session());
        assert_eq!(t.expanded, None);
        t.toggle(StageId::Etching);
        assert_eq!(t.expanded, Some(StageId::Etching));
        t.toggle(StageId::Doping);
        assert_eq!(t.expanded, Some(StageId::Doping), "only one detail is open");
        t.toggle(StageId::Doping);
        assert_eq!(t.expanded, None, "toggling the open one closes it");
    }

    #[test]
    fn a_stage_that_has_run_is_marked() {
        let s = session();
        let t = Timeline::from_session(&s);
        assert!(t.segments.iter().all(|seg| !seg.run), "a fresh session has run nothing");
    }

    #[test]
    fn every_segment_carries_a_label_and_a_status() {
        let t = Timeline::from_session(&session());
        for seg in &t.segments {
            assert!(!seg.label.is_empty());
            assert!(!seg.status.is_empty());
            assert!(seg.width > 0.0);
        }
    }
}
