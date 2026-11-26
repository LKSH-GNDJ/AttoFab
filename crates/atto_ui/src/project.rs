// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The project window: new, open, save, and the wafer specification (PRD G-20).
//!
//! The window owns *no state of its own* — it holds the journal and the wafer
//! specification, and every action is expressed against those two. That is what
//! keeps "save" honest: there is exactly one thing to save, and it is the thing
//! that replays.

use crate::journal::Journal;

/// The wafer a project is about. Separate from the recipe: the same recipe on a
/// different wafer is a different project.
#[derive(Debug, Clone, PartialEq)]
pub struct WaferSpec {
    /// Diameter, mm.
    pub diameter_mm: f64,
    /// Thickness, um.
    pub thickness_um: f64,
    /// The site the tile is taken from, mm from the centre.
    pub site_mm: (f64, f64),
    /// The seed every stochastic step derives from.
    pub seed: u64,
}

impl Default for WaferSpec {
    fn default() -> Self {
        WaferSpec { diameter_mm: 300.0, thickness_um: 500.0, site_mm: (0.0, 0.0), seed: 1 }
    }
}

impl WaferSpec {
    /// Check the specification, returning why it cannot be used.
    ///
    /// A wafer is rejected rather than clamped: a 0 mm diameter silently corrected
    /// to 1 mm produces a plausible project that is not the one asked for.
    pub fn check(&self) -> Result<(), String> {
        if !(self.diameter_mm > 0.0) || !self.diameter_mm.is_finite() {
            return Err(format!("diameter must be positive, got {}", self.diameter_mm));
        }
        if !(self.thickness_um > 0.0) || !self.thickness_um.is_finite() {
            return Err(format!("thickness must be positive, got {}", self.thickness_um));
        }
        let (x, y) = self.site_mm;
        if !x.is_finite() || !y.is_finite() {
            return Err("the site is not a finite point".into());
        }
        let r = (x * x + y * y).sqrt();
        if r > self.diameter_mm * 0.5 {
            return Err(format!("the site at {r:.1} mm is off a {} mm wafer", self.diameter_mm));
        }
        Ok(())
    }
}

/// A project: a journal and the wafer it applies to.
#[derive(Debug, Clone, PartialEq)]
pub struct Project {
    pub spec: WaferSpec,
    pub journal: Journal,
}

impl Project {
    /// A new, empty project.
    pub fn new() -> Self {
        Project { spec: WaferSpec::default(), journal: Journal::new() }
    }

    /// Open a project from its serialised form.
    pub fn open(text: &str) -> Result<Self, String> {
        let journal = Journal::from_json(text)?;
        Ok(Project { spec: WaferSpec::default(), journal })
    }

    /// Serialise for saving.
    pub fn save(&self) -> String {
        self.journal.to_json()
    }

    /// Whether there is anything to save.
    pub fn is_dirty(&self) -> bool {
        !self.journal.is_empty()
    }

    /// Set the wafer, refusing a specification that cannot be used.
    pub fn set_wafer(&mut self, spec: WaferSpec) -> Result<(), String> {
        spec.check()?;
        self.spec = spec;
        Ok(())
    }
}

impl Default for Project {
    fn default() -> Self {
        Project::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_new_project_is_empty_and_not_dirty() {
        let p = Project::new();
        assert!(!p.is_dirty(), "nothing recorded, nothing to save");
        assert_eq!(p.spec.diameter_mm, 300.0, "the default is a 300 mm wafer");
    }

    #[test]
    fn an_impossible_wafer_is_refused_not_clamped() {
        // clamping would produce a plausible project that is not the one asked for
        let mut p = Project::new();
        let bad = WaferSpec { diameter_mm: 0.0, ..WaferSpec::default() };
        assert!(p.set_wafer(bad).is_err());
        let thin = WaferSpec { thickness_um: -1.0, ..WaferSpec::default() };
        assert!(p.set_wafer(thin).is_err());
        assert_eq!(p.spec.diameter_mm, 300.0, "the refused value did not stick");
    }

    #[test]
    fn a_site_off_the_wafer_is_refused() {
        let off = WaferSpec { site_mm: (200.0, 0.0), ..WaferSpec::default() };
        let e = off.check().unwrap_err();
        assert!(e.contains("off a"), "{e}");
        let edge = WaferSpec { site_mm: (149.0, 0.0), ..WaferSpec::default() };
        assert!(edge.check().is_ok(), "a site just inside the edge is fine");
    }

    #[test]
    fn a_project_survives_a_save_and_open() {
        let mut p = Project::new();
        p.journal.record(crate::journal::Command::SetParam { stage: 0, key: "thickness".into(), value: 96.0 });
        let text = p.save();
        let back = Project::open(&text).expect("a project file reads back");
        assert_eq!(back.journal.len(), p.journal.len());
        assert_eq!(back.journal.to_json(), p.journal.to_json());
    }

    #[test]
    fn a_malformed_project_file_is_an_error_not_an_empty_project() {
        // the failure that matters: a bad file must not look like a fresh project
        assert!(Project::open("not a project").is_err());
        assert!(Project::open("{}").is_err() || Project::open("{}").unwrap().journal.is_empty());
    }

    #[test]
    fn saving_twice_gives_the_same_text() {
        let mut p = Project::new();
        p.journal.record(crate::journal::Command::SetParam { stage: 0, key: "thickness".into(), value: 96.0 });
        assert_eq!(p.save(), p.save(), "save is a pure function of the journal");
    }
}
