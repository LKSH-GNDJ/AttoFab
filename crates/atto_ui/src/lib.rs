// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The presentation layer, headless (build PRD Part II, G-01 to G-29).
//!
//! Toolkit-agnostic by construction: every panel is a plain data view model
//! derived from session state, so the whole interface can be driven and tested
//! without a window. Three of the architecture PRD's promises live here rather
//! than in a shell:
//!
//! * **G-02 state is always visible** — a stage's status is *derived* from the
//!   fingerprints of its parameters, never stored as a flag that can go stale.
//! * **G-03 every number is one click from its definition** — each row carries
//!   the model, the parameters and the uncertainty behind it.
//! * **G-13 an empty panel says why** — a measurement that was not taken
//!   carries a reason, never a zero.

pub mod compare;
pub mod console;
pub mod doe;
pub mod docking;
pub mod evidence;
pub mod export;
pub mod gdspreview;
pub mod inspector;
pub mod jobs;
pub mod layout;
pub mod panes;
pub mod plots;
pub mod project;
pub mod provenance;
pub mod reorder;
pub mod preview;
pub mod scaling;
pub mod stack;
pub mod strings;
pub mod timeline;
pub mod journal;
pub mod mini_json;

use atto_measure::{Measurement, Status};

// ---------------------------------------------------------------------------
// Recipe
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
/// The stages a recipe is made of.
pub enum StageId {
    Prep,
    Deposition,
    Lithography,
    Etching,
    Doping,
    Anneal,
    Oxidation,
    Metrology,
}

impl StageId {
    /// Every stage, in process order.
    pub const ALL: [StageId; 8] = [
        StageId::Prep,
        StageId::Deposition,
        StageId::Lithography,
        StageId::Etching,
        StageId::Doping,
        StageId::Anneal,
        StageId::Oxidation,
        StageId::Metrology,
    ];

    /// The name used in a report.
    pub fn label(self) -> &'static str {
        match self {
            StageId::Prep => "Substrate prep",
            StageId::Deposition => "Deposition",
            StageId::Lithography => "Lithography",
            StageId::Etching => "Etching",
            StageId::Doping => "Doping",
            StageId::Anneal => "Anneal",
            StageId::Oxidation => "Oxidation",
            StageId::Metrology => "Metrology",
        }
    }
}

/// One stage's parameters, in a fixed order so a fingerprint is stable.
#[derive(Debug, Clone, PartialEq)]
pub struct StageSpec {
    pub id: StageId,
    pub params: Vec<(String, f64)>,
    pub enabled: bool,
}

impl StageSpec {
    /// A stage specification for an id.
    pub fn new(id: StageId) -> Self {
        StageSpec { id, params: Vec::new(), enabled: true }
    }

    /// A builder form of `set`.
    pub fn with(mut self, key: &str, value: f64) -> Self {
        self.set(key, value);
        self
    }

    /// Set a knob, recording it in the history.
    pub fn set(&mut self, key: &str, value: f64) {
        match self.params.iter_mut().find(|(k, _)| k == key) {
            Some(p) => p.1 = value,
            None => self.params.push((key.to_string(), value)),
        }
        self.params.sort_by(|a, b| a.0.cmp(&b.0));
    }

    /// A knob's value.
    pub fn get(&self, key: &str) -> Option<f64> {
        self.params.iter().find(|(k, _)| k == key).map(|(_, v)| *v)
    }

    /// Everything that can change this stage's result, folded into one number.
    pub fn fingerprint(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |b: &[u8]| {
            for x in b {
                h ^= *x as u64;
                h = h.wrapping_mul(0x0000_0100_0000_01B3);
            }
        };
        feed(&[self.id as u8, self.enabled as u8]);
        for (k, v) in &self.params {
            feed(k.as_bytes());
            feed(&v.to_bits().to_le_bytes());
        }
        h
    }
}

/// Derived, never stored.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StageStatus {
    NotRun,
    /// Cannot run because something upstream has not run (build PRD G-25).
    /// A blocked stage is not a failed one: nothing is wrong with it.
    Blocked,
    Committed,
    /// Its own parameters changed, or an upstream stage's did.
    Stale,
    Failed,
}

impl StageStatus {
    /// The one-word status.
    pub fn word(self) -> &'static str {
        // G-27: the words live in the resource table, not here
        match self {
            StageStatus::NotRun => crate::strings::STATUS_NOT_RUN,
            StageStatus::Blocked => crate::strings::STATUS_BLOCKED,
            StageStatus::Committed => crate::strings::STATUS_COMMITTED,
            StageStatus::Stale => crate::strings::STATUS_STALE,
            StageStatus::Failed => crate::strings::STATUS_FAILED,
        }
    }

    /// Semantic colour, so a shell maps it without inventing one (PRD G-04).
    /// A word is always shown alongside: colour alone fails for the roughly one
    /// man in twelve with a colour-vision deficiency.
    pub fn colour(self) -> Colour {
        match self {
            StageStatus::NotRun => Colour::Grey,
            StageStatus::Blocked => Colour::Grey,
            StageStatus::Committed => Colour::Green,
            StageStatus::Stale => Colour::Amber,
            StageStatus::Failed => Colour::Red,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// The colour a status is drawn in.
pub enum Colour {
    Grey,
    Green,
    Amber,
    Red,
}

/// What a stage produced, for the inspector and for compare mode.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct StageResult {
    pub measurements: Vec<Measurement>,
    /// Model identity, for the provenance popover.
    pub model: String,
    /// Parameter ids the model read.
    pub parameters: Vec<String>,
    pub inputs_fingerprint: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Why a session operation was refused.
pub enum UiError {
    NoSuchStage,
    NotEnabled,
    UndoEmpty,
    RedoEmpty,
}

/// The evaluator the shell supplies. Tests substitute a fake.
pub type Evaluator<'a> = &'a dyn Fn(&StageSpec, u64) -> Result<StageResult, String>;

// ---------------------------------------------------------------------------
// Session
// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct Snapshot {
    stages: Vec<StageSpec>,
    committed: Vec<Option<u64>>,
}

#[derive(Debug, Clone)]
/// The interface state: recipe, selection, history and fingerprint.
pub struct Session {
    stages: Vec<StageSpec>,
    /// The fingerprint each stage was committed at. `None` means never run.
    committed: Vec<Option<u64>>,
    results: Vec<Option<StageResult>>,
    failed: Vec<bool>,
    /// The message from the last failed commit, per stage (PRD G-07).
    failure: Vec<Option<String>>,
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
    selected: usize,
    /// Maximum number of undo steps retained. A session that runs for hours
    /// would otherwise keep every snapshot alive (build PRD X-06).
    history_limit: usize,
}

impl Default for Session {
    fn default() -> Self {
        Session::new()
    }
}

impl Session {
    /// A session with the default recipe.
    pub fn new() -> Self {
        let stages = StageId::ALL.iter().map(|id| StageSpec::new(*id)).collect::<Vec<_>>();
        let n = stages.len();
        Session {
            stages,
            committed: vec![None; n],
            results: vec![None; n],
            failed: vec![false; n],
            failure: vec![None; n],
            undo: Vec::new(),
            redo: Vec::new(),
            selected: 0,
            history_limit: 64,
        }
    }

    /// The recipe's stages.
    pub fn stages(&self) -> &[StageSpec] {
        &self.stages
    }

    /// The selected stage index.
    pub fn selected(&self) -> usize {
        self.selected
    }

    /// Select a stage by index.
    pub fn select(&mut self, i: usize) -> Result<(), UiError> {
        if i >= self.stages.len() {
            return Err(UiError::NoSuchStage);
        }
        self.selected = i;
        Ok(())
    }

    /// Status is computed from the fingerprints, so it cannot go stale.
    pub fn status(&self, i: usize) -> StageStatus {
        if i >= self.stages.len() {
            return StageStatus::NotRun;
        }
        if self.failed[i] {
            return StageStatus::Failed;
        }
        if self.committed[i].is_none() {
            // never run itself: blocked if an upstream stage that would feed it
            // has not run either, otherwise simply not run
            let upstream_missing = (0..i).any(|j| self.committed[j].is_none() && self.stages[j].enabled);
            return if upstream_missing { StageStatus::Blocked } else { StageStatus::NotRun };
        }
        // this stage, or anything upstream that actually ran, has changed
        // since its commit. A stage that was never run contributed nothing, so
        // it cannot make anything stale.
        for j in 0..=i {
            if let Some(c) = self.committed[j] {
                if c != self.stages[j].fingerprint() {
                    return StageStatus::Stale;
                }
            }
        }
        StageStatus::Committed
    }

    /// Why the last commit of a stage failed, if it did.
    pub fn failure_reason(&self, i: usize) -> Option<&str> {
        self.failure.get(i).and_then(|f| f.as_deref())
    }

    /// The result of a stage, when it has one.
    /// The fingerprint a stage was last committed at, if it ever was.
    ///
    /// Two sessions with different values here were run against different inputs,
    /// which is what compare mode reports (PRD G-09).
    pub fn committed_at(&self, i: usize) -> Option<u64> {
        self.committed.get(i).copied().flatten()
    }

    /// The result of a stage, when it has one.
    pub fn result(&self, i: usize) -> Option<&StageResult> {
        self.results.get(i).and_then(|r| r.as_ref())
    }

    /// Change a parameter. The status of this stage and everything downstream
    /// follows from the fingerprints; nothing is marked dirty.
    pub fn set_param(&mut self, i: usize, key: &str, value: f64) -> Result<(), UiError> {
        if i >= self.stages.len() {
            return Err(UiError::NoSuchStage);
        }
        self.push_undo();
        self.stages[i].set(key, value);
        Ok(())
    }

    /// Turn a stage on or off, refusing if it cannot run.
    pub fn set_enabled(&mut self, i: usize, on: bool) -> Result<(), UiError> {
        if i >= self.stages.len() {
            return Err(UiError::NoSuchStage);
        }
        self.push_undo();
        self.stages[i].enabled = on;
        Ok(())
    }

    /// Commit a stage. `Ok(None)` means the model ran and failed, which is a
    /// visible state; `Err` means the request itself was invalid.
    pub fn commit(&mut self, i: usize, eval: Evaluator<'_>) -> Result<Option<&StageResult>, UiError> {
        if i >= self.stages.len() {
            return Err(UiError::NoSuchStage);
        }
        if !self.stages[i].enabled {
            return Err(UiError::NotEnabled);
        }
        self.push_undo();
        let fp = self.stages[i].fingerprint();
        match eval(&self.stages[i], fp) {
            Ok(r) => {
                self.results[i] = Some(r);
                self.committed[i] = Some(fp);
                self.failed[i] = false;
                self.failure[i] = None;
            }
            Err(e) => {
                // the message is kept, not discarded: the recipe row shows it
                self.failure[i] = Some(e);
                // a failure is visible state, not a silent no-op
                self.failed[i] = true;
                self.committed[i] = Some(fp);
                self.results[i] = None;
            }
        }
        Ok(self.results[i].as_ref())
    }

    fn push_undo(&mut self) {
        self.undo.push(Snapshot { stages: self.stages.clone(), committed: self.committed.clone() });
        // drop the oldest step once the limit is reached, so memory is bounded
        // and the oldest history is what goes
        while self.undo.len() > self.history_limit {
            self.undo.remove(0);
        }
        self.redo.clear();
    }

    /// Set how many undo steps are kept. Zero disables undo.
    pub fn set_history_limit(&mut self, n: usize) {
        self.history_limit = n;
        while self.undo.len() > n {
            self.undo.remove(0);
        }
        while self.redo.len() > n {
            self.redo.remove(0);
        }
    }

    /// The most undo steps the session keeps.
    pub fn history_limit(&self) -> usize {
        self.history_limit
    }

    /// How many undo steps are held.
    pub fn history_depth(&self) -> usize {
        self.undo.len()
    }

    /// Undo the last change.
    pub fn undo(&mut self) -> Result<(), UiError> {
        let s = self.undo.pop().ok_or(UiError::UndoEmpty)?;
        self.redo.push(Snapshot { stages: self.stages.clone(), committed: self.committed.clone() });
        self.stages = s.stages;
        self.committed = s.committed;
        Ok(())
    }

    /// Reapply the last undone change.
    pub fn redo(&mut self) -> Result<(), UiError> {
        let s = self.redo.pop().ok_or(UiError::RedoEmpty)?;
        self.undo.push(Snapshot { stages: self.stages.clone(), committed: self.committed.clone() });
        self.stages = s.stages;
        self.committed = s.committed;
        Ok(())
    }

    /// Whether anything can be undone.
    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }

    /// Whether anything can be redone.
    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }

    /// Everything that determines what the panels show, folded into one number.
    /// A replay that lands on the same id produced the same session.
    pub fn fingerprint(&self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut feed = |b: &[u8]| {
            for x in b {
                h ^= *x as u64;
                h = h.wrapping_mul(0x0000_0100_0000_01B3);
            }
        };
        for (i, s) in self.stages.iter().enumerate() {
            feed(&s.fingerprint().to_le_bytes());
            feed(&self.committed[i].unwrap_or(0).to_le_bytes());
            feed(&[self.failed[i] as u8, self.stages[i].enabled as u8]);
        }
        feed(&(self.selected as u64).to_le_bytes());
        h
    }

    /// The recipe panel's rows.
    pub fn recipe_rows(&self) -> Vec<RecipeRow> {
        self.stages
            .iter()
            .enumerate()
            .map(|(i, s)| RecipeRow {
                index: i,
                label: s.id.label(),
                status: self.status(i),
                selected: i == self.selected,
                enabled: s.enabled,
            })
            .collect()
    }
}

#[derive(Debug, Clone, PartialEq)]
/// One recipe row, ready to draw.
pub struct RecipeRow {
    pub index: usize,
    pub label: &'static str,
    pub status: StageStatus,
    pub selected: bool,
    pub enabled: bool,
}

// ---------------------------------------------------------------------------
// Provenance and the inspector
// ---------------------------------------------------------------------------

/// Everything behind one displayed number (G-03).
#[derive(Debug, Clone, PartialEq)]
pub struct Provenance {
    pub model: String,
    pub parameters: Vec<String>,
    pub inputs_fingerprint: u64,
    /// One sigma, in the row's unit, when the model declares one.
    pub uncertainty: Option<f64>,
    pub note: String,
}

impl Provenance {
    /// A one-line description of the session.
    pub fn summary(&self) -> String {
        let p = if self.parameters.is_empty() {
            "no parameters".to_string()
        } else {
            self.parameters.join(", ")
        };
        format!("{} · {} · fp {:016x}", self.model, p, self.inputs_fingerprint)
    }
}

/// One row of the measurement card.
#[derive(Debug, Clone, PartialEq)]
pub struct MeasureRow {
    pub label: String,
    pub unit: String,
    pub value: Option<f64>,
    pub nominal: f64,
    pub tolerance: f64,
    pub status: Status,
    /// Present when the measurement was not taken.
    pub reason: Option<String>,
    pub provenance: Provenance,
}

impl MeasureRow {
    /// The value as a report shows it, with its unit.
    pub fn display(&self) -> String {
        match self.value {
            Some(v) => format!("{v:.2} {}", self.unit),
            None => format!("not measured — {}", self.reason.clone().unwrap_or_else(|| "no reason given".into())),
        }
    }
}

/// Build the measurement card from a stage's result.
pub fn measure_rows(result: &StageResult) -> Vec<MeasureRow> {
    result
        .measurements
        .iter()
        .map(|m| MeasureRow {
            label: m.name.clone(),
            unit: m.unit.clone(),
            value: m.value,
            nominal: m.nominal,
            tolerance: m.tolerance,
            status: m.status,
            reason: m.reason.clone(),
            provenance: Provenance {
                model: result.model.clone(),
                parameters: result.parameters.clone(),
                inputs_fingerprint: result.inputs_fingerprint,
                uncertainty: if m.sigma > 0.0 { Some(m.sigma) } else { None },
                // a measured row has no reason to give; an unmeasured one always
                // does, and the absence is a real empty note rather than a
                // swallowed error
                note: match &m.reason {
                    Some(r) => r.clone(),
                    None => String::new(),
                },
            },
        })
        .collect()
}

/// A not-measured row must carry a reason: this is the rule, checkable.
pub fn rows_are_honest(rows: &[MeasureRow]) -> bool {
    rows.iter().all(|r| match r.value {
        Some(_) => true,
        None => r.reason.as_ref().map(|s| !s.trim().is_empty()).unwrap_or(false),
    })
}

// ---------------------------------------------------------------------------
// Compare
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, PartialEq)]
/// The difference between two sessions.
pub struct Delta {
    pub label: String,
    pub a: Option<f64>,
    pub b: Option<f64>,
    pub delta: Option<f64>,
    pub unit: String,
}

impl Delta {
    /// The value as a report shows it, with its unit.
    pub fn display(&self) -> String {
        match (self.a, self.b) {
            (Some(a), Some(b)) => format!("{a:.2} → {b:.2} ({:+.2}) {}", b - a, self.unit),
            (Some(a), None) => format!("{a:.2} → not measured {}", self.unit),
            (None, Some(b)) => format!("not measured → {b:.2} {}", self.unit),
            (None, None) => "not measured in either state".to_string(),
        }
    }
}

/// Two states side by side. Rows are matched by label.
pub fn compare(a: &[MeasureRow], b: &[MeasureRow]) -> Vec<Delta> {
    let mut out = Vec::new();
    for r in a {
        let other = b.iter().find(|x| x.label == r.label);
        out.push(Delta {
            label: r.label.clone(),
            a: r.value,
            b: other.and_then(|x| x.value),
            delta: match (r.value, other.and_then(|x| x.value)) {
                (Some(x), Some(y)) => Some(y - x),
                _ => None,
            },
            unit: r.unit.clone(),
        });
    }
    for r in b {
        if !a.iter().any(|x| x.label == r.label) {
            out.push(Delta { label: r.label.clone(), a: None, b: r.value, delta: None, unit: r.unit.clone() });
        }
    }
    out
}

/// Named, quantified progress (build PRD G-23).
///
/// "Working..." is not progress. A step carries a name and a count, so a shell
/// can show `Etching — 3 of 8` and a user can tell a slow step from a stuck one.
#[derive(Debug, Clone, PartialEq)]
pub struct Progress {
    pub label: String,
    pub done: usize,
    pub total: usize,
}

impl Progress {
    /// Progress for a named stage of a known length.
    pub fn new(label: &str, total: usize) -> Self {
        Progress { label: label.to_string(), done: 0, total }
    }

    /// Advance the progress by one step.
    pub fn step(&mut self) {
        if self.done < self.total {
            self.done += 1;
        }
    }

    /// How far through, 0 to 1.
    pub fn fraction(&self) -> f64 {
        if self.total == 0 {
            1.0
        } else {
            self.done as f64 / self.total as f64
        }
    }

    /// Whether the progress has reached its total.
    pub fn is_done(&self) -> bool {
        self.done >= self.total
    }

    /// The text a shell shows: never just a percentage, always the step name.
    pub fn text(&self) -> String {
        format!("{} — {} of {}", self.label, self.done, self.total)
    }
}

/// One entry in the command palette (build PRD G-19).
#[derive(Debug, Clone, PartialEq)]
pub struct PaletteEntry {
    pub id: String,
    pub label: String,
    /// Key binding, empty when the command has none.
    pub binding: String,
}

/// Every command the interface offers, with its binding.
///
/// Derived from the same table the help sheet uses, so a command cannot exist
/// without a binding being listed or vice versa.
pub fn palette() -> Vec<PaletteEntry> {
    let mut out: Vec<PaletteEntry> = SHORTCUTS
        .iter()
        .map(|(k, d)| PaletteEntry { id: d.to_lowercase().replace(' ', "-"), label: d.to_string(), binding: k.to_string() })
        .collect();
    for (i, id) in StageId::ALL.iter().enumerate() {
        out.push(PaletteEntry {
            id: format!("commit-{}", id.label().to_lowercase().replace(' ', "-")),
            label: format!("Commit {}", id.label()),
            // G-22: every action is reachable from the keyboard, so every commit
            // carries a binding rather than being palette-only
            binding: format!("Alt+{}", i + 1),
        });
    }
    out
}

/// The full binding sheet: every action that has a key.
pub fn binding_sheet() -> Vec<PaletteEntry> {
    palette().into_iter().filter(|e| !e.binding.is_empty()).collect()
}

/// Actions a keyboard cannot reach. Empty when G-22 holds.
pub fn unreachable_actions() -> Vec<String> {
    palette().into_iter().filter(|e| e.binding.is_empty()).map(|e| e.id).collect()
}

/// The keyboard map, as data so a shell and the help sheet cannot disagree.
pub const SHORTCUTS: &[(&str, &str)] = &[
    ("Ctrl+Z", "Undo"),
    ("Ctrl+Shift+Z", "Redo"),
    ("Space", "Commit the selected stage"),
    ("Escape", "Cancel the running job"),
    ("1..8", "Select a stage"),
    ("Ctrl+K", "Command palette"),
    ("Ctrl+D", "Compare with the previous commit"),
    ("Ctrl+Enter", "Run the whole recipe"),
    ("Alt+1..8", "Commit a stage by number"),
];

#[cfg(test)]
mod tests {
    use super::*;
    use atto_measure::Measurement;

    fn fake(ok: bool) -> impl Fn(&StageSpec, u64) -> Result<StageResult, String> {
        move |spec: &StageSpec, fp: u64| {
            if !ok {
                return Err("the model failed".into());
            }
            let t = spec.get("thickness").unwrap_or(100.0);
            Ok(StageResult {
                measurements: vec![
                    Measurement::against("Film thickness", "nm", t, 100.0, 10.0, 1.0),
                    Measurement::missing("Junction depth", "nm", "no dopant in this window"),
                ],
                model: "deposition.directional".into(),
                parameters: spec.params.iter().map(|(k, _)| k.clone()).collect(),
                inputs_fingerprint: fp,
            })
        }
    }

    fn committed() -> Session {
        let mut s = Session::new();
        s.set_param(1, "thickness", 96.0).unwrap();
        s.commit(1, &fake(true)).unwrap();
        s
    }

    fn committed_result(s: &Session) -> &StageResult {
        s.result(1).expect("stage 1 has a result")
    }

    #[test]
    fn a_fresh_stage_is_not_run_and_a_committed_one_is_committed() {
        let s = committed();
        assert_eq!(s.status(1), StageStatus::Committed);
        assert_eq!(s.status(0), StageStatus::NotRun, "prep was never run");
    }

    #[test]
    fn changing_a_parameter_makes_it_stale_and_changing_it_back_restores_it() {
        let mut s = committed();
        s.set_param(1, "thickness", 120.0).unwrap();
        assert_eq!(s.status(1), StageStatus::Stale);
        s.set_param(1, "thickness", 96.0).unwrap();
        assert_eq!(s.status(1), StageStatus::Committed, "staleness is derived, so it clears");
    }

    #[test]
    fn an_upstream_change_stales_everything_downstream() {
        let mut s = Session::new();
        for i in 0..4 {
            s.commit(i, &fake(true)).unwrap();
        }
        assert_eq!(s.status(3), StageStatus::Committed);
        s.set_param(1, "thickness", 50.0).unwrap();
        assert_eq!(s.status(0), StageStatus::Committed, "upstream is untouched");
        assert_eq!(s.status(1), StageStatus::Stale);
        assert_eq!(s.status(2), StageStatus::Stale, "downstream follows");
        assert_eq!(s.status(3), StageStatus::Stale);
    }

    #[test]
    fn a_failure_is_visible_not_silent() {
        let mut s = Session::new();
        s.set_param(1, "thickness", 96.0).unwrap();
        let out = s.commit(1, &fake(false)).unwrap();
        assert!(out.is_none(), "a failed evaluation returns no result");
        assert_eq!(s.status(1), StageStatus::Failed);
        assert!(s.result(1).is_none());
    }

    #[test]
    fn undo_and_redo_restore_parameters_and_status() {
        let mut s = committed();
        s.set_param(1, "thickness", 140.0).unwrap();
        assert_eq!(s.status(1), StageStatus::Stale);
        s.undo().unwrap();
        assert_eq!(s.status(1), StageStatus::Committed);
        assert_eq!(s.stages()[1].get("thickness"), Some(96.0));
        s.redo().unwrap();
        assert_eq!(s.stages()[1].get("thickness"), Some(140.0));
        assert_eq!(s.status(1), StageStatus::Stale);
    }

    #[test]
    fn undo_on_an_empty_history_is_an_error() {
        let mut s = Session::new();
        assert_eq!(s.undo(), Err(UiError::UndoEmpty));
        assert_eq!(s.redo(), Err(UiError::RedoEmpty));
        assert!(!s.can_undo());
    }

    #[test]
    fn the_measurement_card_carries_provenance_and_a_reason() {
        let s = committed();
        let rows = measure_rows(committed_result(&s));
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].value, Some(96.0));
        assert_eq!(rows[0].provenance.model, "deposition.directional");
        assert!(rows[0].provenance.parameters.contains(&"thickness".to_string()));
        // the missing one explains itself
        assert!(rows[1].value.is_none());
        assert!(rows[1].display().contains("no dopant"));
        assert!(rows_are_honest(&rows));
    }

    #[test]
    fn the_card_shows_nominal_tolerance_uncertainty_and_status_together() {
        // G-13 asks for all four at once, not a subset: a number without its
        // tolerance cannot be judged, and a tolerance without its uncertainty
        // cannot be trusted.
        let res = StageResult {
            measurements: vec![atto_measure::Measurement::against("Film thickness", "nm", 96.0, 100.0, 10.0, 1.0)],
            model: "atto_geom::grow_film".into(),
            parameters: vec!["thickness".into()],
            inputs_fingerprint: 1,
        };
        let rows = measure_rows(&res);
        let r = &rows[0];
        assert_eq!(r.nominal, 100.0, "nominal");
        assert_eq!(r.tolerance, 10.0, "tolerance");
        assert_eq!(r.provenance.uncertainty, Some(1.0), "uncertainty");
        assert!(!r.status.label().is_empty(), "status");
        // and all four are present on the same row, not spread across the card
        assert!(r.value.is_some(), "the value itself");
    }

    #[test]
    fn a_row_without_a_reason_is_dishonest() {
        let rows = vec![MeasureRow {
            label: "x".into(),
            unit: "nm".into(),
            value: None,
            nominal: 0.0,
            tolerance: 0.0,
            status: Status::NotMeasured,
            reason: None,
            provenance: Provenance { model: "m".into(), parameters: vec![], inputs_fingerprint: 0, uncertainty: None, note: String::new() },
        }];
        assert!(!rows_are_honest(&rows), "a missing value must say why");
    }

    #[test]
    fn compare_reports_the_delta_and_handles_missing_sides() {
        let mut s = committed();
        let before = measure_rows(committed_result(&s));
        s.set_param(1, "thickness", 120.0).unwrap();
        s.commit(1, &fake(true)).unwrap();
        let after = measure_rows(committed_result(&s));
        let d = compare(&before, &after);
        let film = d.iter().find(|x| x.label == "Film thickness").unwrap();
        assert_eq!(film.delta, Some(24.0));
        assert!(film.display().contains("+24.00"));
        let jd = d.iter().find(|x| x.label == "Junction depth").unwrap();
        assert_eq!(jd.delta, None);
        assert!(jd.display().contains("not measured"));
    }

    #[test]
    fn the_undo_history_is_bounded() {
        let mut s = Session::new();
        s.set_history_limit(5);
        for i in 0..50 {
            s.set_param(1, "thickness", 100.0 + i as f64).unwrap();
        }
        assert_eq!(s.history_depth(), 5, "only five steps are kept");
        let mut undos = 0;
        while s.undo().is_ok() {
            undos += 1;
        }
        assert_eq!(undos, 5, "and only five can be undone");
    }

    #[test]
    fn a_zero_limit_disables_undo() {
        let mut s = Session::new();
        s.set_history_limit(0);
        s.set_param(1, "thickness", 120.0).unwrap();
        assert_eq!(s.history_depth(), 0);
        assert_eq!(s.undo(), Err(UiError::UndoEmpty));
    }

    #[test]
    fn a_stage_with_no_upstream_run_is_blocked_not_merely_unrun() {
        let mut s = Session::new();
        // stage 0 has no upstream, so it is simply not run
        assert_eq!(s.status(0), StageStatus::NotRun);
        // stage 4 cannot run until 0..3 have
        assert_eq!(s.status(4), StageStatus::Blocked);
        s.commit(0, &fake(true)).unwrap();
        assert_eq!(s.status(1), StageStatus::NotRun, "one upstream is now done");
    }

    #[test]
    fn status_carries_a_word_and_a_colour() {
        // G-04: never a colour alone
        assert_eq!(StageStatus::Failed.word(), "Failed");
        assert_eq!(StageStatus::Failed.colour(), Colour::Red);
        assert_eq!(StageStatus::Stale.colour(), Colour::Amber);
        assert_eq!(StageStatus::Committed.colour(), Colour::Green);
        assert_eq!(StageStatus::Blocked.colour(), Colour::Grey);
        assert!(!StageStatus::NotRun.word().is_empty());
    }

    #[test]
    fn a_failed_stage_keeps_its_message() {
        let mut s = Session::new();
        s.set_param(1, "thickness", 96.0).unwrap();
        assert_eq!(s.failure_reason(1), None);
        s.commit(1, &fake(false)).unwrap();
        assert_eq!(s.status(1), StageStatus::Failed);
        assert_eq!(s.failure_reason(1), Some("the model failed"), "G-07: the message survives");
        // a successful commit clears it
        s.commit(1, &fake(true)).unwrap();
        assert_eq!(s.failure_reason(1), None);
    }

    #[test]
    fn progress_is_named_and_quantified() {
        let mut p = Progress::new("Etching", 8);
        assert_eq!(p.text(), "Etching — 0 of 8");
        assert!((p.fraction() - 0.0).abs() < 1e-12);
        for _ in 0..3 {
            p.step();
        }
        assert_eq!(p.text(), "Etching — 3 of 8");
        assert!((p.fraction() - 0.375).abs() < 1e-12);
        for _ in 0..10 {
            p.step();
        }
        assert!(p.is_done(), "it cannot overrun its total");
        assert_eq!(p.done, 8);
    }

    #[test]
    fn the_palette_covers_every_shortcut_and_every_stage() {
        let p = palette();
        for (binding, desc) in SHORTCUTS {
            assert!(p.iter().any(|e| e.binding == *binding && e.label == *desc), "missing {binding}");
        }
        for id in StageId::ALL {
            let want = format!("Commit {}", id.label());
            assert!(p.iter().any(|e| e.label == want), "missing {want}");
        }
        let mut ids: Vec<&str> = p.iter().map(|e| e.id.as_str()).collect();
        ids.sort();
        let n = ids.len();
        ids.dedup();
        assert_eq!(ids.len(), n, "command ids must be unique");
    }

    #[test]
    fn the_recipe_panel_reflects_selection_and_status() {
        let mut s = committed();
        s.select(3).unwrap();
        let rows = s.recipe_rows();
        assert_eq!(rows.len(), 8);
        assert!(rows[3].selected);
        assert_eq!(rows[1].status, StageStatus::Committed);
        assert_eq!(rows[0].status, StageStatus::NotRun);
    }

    #[test]
    fn a_disabled_stage_cannot_commit() {
        let mut s = Session::new();
        s.set_enabled(1, false).unwrap();
        assert_eq!(s.commit(1, &fake(true)).err(), Some(UiError::NotEnabled));
    }

    #[test]
    fn every_action_is_reachable_from_the_keyboard() {
        // G-22: an action in the palette with no key is an action half the
        // interface cannot reach
        let missing = unreachable_actions();
        assert!(missing.is_empty(), "actions without a key: {missing:?}");
    }

    #[test]
    fn the_binding_sheet_has_no_duplicate_keys() {
        let sheet = binding_sheet();
        let mut keys: Vec<String> = sheet.iter().map(|e| e.binding.clone()).collect();
        let n = keys.len();
        keys.sort();
        keys.dedup();
        assert_eq!(keys.len(), n, "two actions share a key");
        assert!(!sheet.is_empty());
    }

    #[test]
    fn every_shortcut_is_documented_and_unique() {
        let mut keys: Vec<&str> = SHORTCUTS.iter().map(|(k, _)| *k).collect();
        keys.sort();
        let n = keys.len();
        keys.dedup();
        assert_eq!(keys.len(), n, "a binding is listed twice");
        assert!(SHORTCUTS.iter().all(|(_, d)| !d.is_empty()));
    }
}
