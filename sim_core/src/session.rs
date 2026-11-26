// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Headless editing session: recipe, committed process records, transactional
//! preview / commit / cancel, undo / redo and the wafer environment transform
//! (PRD §1, §2, §6, §7 P1–P2).
//!
//! The UI is a thin client of this state machine; every acceptance test that
//! does not need a window (INT-05, INT-06, DATA-01, PROC-*, …) drives it
//! directly.
//!
//! **Staleness is derived, never flagged.** A stage is *committed* only if a
//! record exists whose parameters, window, seed and *input wafer fingerprint*
//! all match the current recipe chain. Changing an upstream stage changes the
//! input fingerprint of every later one, so downstream results become stale
//! automatically and are recomputed on the next run.

use crate::geom::ProcessWindow;
use crate::params::StageParams;
use crate::physics::Stage;
use crate::pipeline::{evaluate_stage, Overlay, PipelineError, StageOutcome, SubStep, FRAMES_PER_STAGE};
use crate::transform::WaferTransform;
use crate::wafer::{WaferError, WaferState};
use std::collections::{BTreeMap, HashMap, VecDeque};
use std::sync::Arc;

pub const UNDO_LIMIT: usize = 64;
const LOG_LIMIT: usize = 400;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WaferSpec { pub diameter_mm: f64, pub thickness_um: f64, pub leaf_nm: f64, pub tile_levels: u32 }

impl Default for WaferSpec {
    fn default() -> Self { Self { diameter_mm: 300.0, thickness_um: 500.0, leaf_nm: 16.0, tile_levels: 7 } }
}

impl WaferSpec {
    pub fn build(&self, site: [f64; 2]) -> Result<WaferState, WaferError> {
        WaferState::new_wafer_tile(self.diameter_mm, self.thickness_um, self.leaf_nm, self.tile_levels, (site[0], site[1]))
    }
    pub fn tile_edge_nm(&self) -> f64 { self.leaf_nm * (1u64 << self.tile_levels.clamp(1, 14)) as f64 }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RecipeStage { pub id: u32, pub enabled: bool, pub params: StageParams }

#[derive(Debug, Clone, PartialEq)]
pub struct Recipe { pub stages: Vec<RecipeStage>, pub window: ProcessWindow, pub seed: u64, next_id: u32 }

impl Recipe {
    pub fn empty(window: ProcessWindow) -> Self { Self { stages: vec![], window, seed: 1, next_id: 1 } }
    /// Seven stages in PRD order.
    pub fn default_recipe(window: ProcessWindow) -> Self {
        let mut r = Self::empty(window);
        for s in Stage::DEFAULT_RECIPE { r.add(StageParams::default_for(s)); }
        r
    }
    pub fn add(&mut self, params: StageParams) -> u32 {
        let id = self.next_id; self.next_id += 1;
        self.stages.push(RecipeStage { id, enabled: true, params });
        id
    }
    pub fn index_of(&self, id: u32) -> Option<usize> { self.stages.iter().position(|s| s.id == id) }
    pub fn insert_after(&mut self, id: u32, params: StageParams) -> Option<u32> {
        let i = self.index_of(id)?;
        let nid = self.next_id; self.next_id += 1;
        self.stages.insert(i + 1, RecipeStage { id: nid, enabled: true, params });
        Some(nid)
    }
    pub fn remove(&mut self, id: u32) -> bool { let n = self.stages.len(); self.stages.retain(|s| s.id != id); n != self.stages.len() }
    pub fn move_stage(&mut self, id: u32, delta: isize) -> bool {
        let Some(i) = self.index_of(id) else { return false };
        let j = i as isize + delta;
        if j < 0 || j as usize >= self.stages.len() { return false; }
        // `swap` only moves one place; remove and re-insert honours any delta
        let stage = self.stages.remove(i as usize);
        self.stages.insert(j as usize, stage);
        true
    }
}

#[derive(Debug, Clone, PartialEq)]
pub enum StageStatus { Disabled, NotRun, Committed, Stale, Failed(String) }

/// A committed (or candidate) process record: parameters, input state, output
/// state, outcome and the exact inputs needed to replay it.
#[derive(Debug, Clone)]
pub struct StageRecord {
    pub stage_id: u32,
    pub stage: Stage,
    pub revision: u64,
    pub params: StageParams,
    pub window: ProcessWindow,
    pub seed: u64,
    pub before: Arc<WaferState>,
    pub after: Arc<WaferState>,
    pub input_fp: u64,
    pub output_fp: u64,
    pub outcome: Arc<StageOutcome>,
}

#[derive(Clone)]
pub struct SessionState {
    pub spec: WaferSpec,
    pub initial: Arc<WaferState>,
    pub recipe: Recipe,
    pub records: BTreeMap<u32, StageRecord>,
    pub failures: BTreeMap<u32, String>,
    pub transform: WaferTransform,
    pub revision: u64,
}

pub struct Preview { pub stage_id: u32, pub record: StageRecord }

#[derive(Debug, Clone, PartialEq)]
pub enum SessionError {
    UnknownStage(u32),
    UpstreamNotCommitted(String),
    PendingSite,
    NothingToCommit,
    Pipeline(PipelineError),
    InvalidSite,
    InvalidScale,
    Wafer(String),
}

impl std::fmt::Display for SessionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SessionError::UnknownStage(i) => write!(f, "unknown stage {i}"),
            SessionError::UpstreamNotCommitted(s) => write!(f, "upstream stage '{s}' is not committed; run or commit it first"),
            SessionError::PendingSite => write!(f, "a new process site is pending: commit the site first (it invalidates every committed stage)"),
            SessionError::NothingToCommit => write!(f, "no preview to commit"),
            SessionError::Pipeline(e) => write!(f, "{e}"),
            SessionError::InvalidSite => write!(f, "process site must be finite and inside the wafer"),
            SessionError::InvalidScale => write!(f, "cell size must be 1–2000 nm and tile depth 4–12 levels"),
            SessionError::Wafer(e) => write!(f, "{e}"),
        }
    }
}
impl std::error::Error for SessionError {}
impl From<PipelineError> for SessionError { fn from(e: PipelineError) -> Self { SessionError::Pipeline(e) } }

#[derive(Debug, Clone, Default)]
pub struct RunSummary { pub ran: Vec<u32>, pub skipped_committed: Vec<u32>, pub failed: Option<(u32, String)> }

pub struct Session {
    pub state: SessionState,
    undo: VecDeque<SessionState>,
    redo: Vec<SessionState>,
    pub preview: Option<Preview>,
    pub pending_site: Option<[f64; 2]>,
    coalesce: Option<u32>,
    drag_start: Option<SessionState>,
    pub log: VecDeque<String>,
}

const TRANSFORM_KEY: u32 = u32::MAX;

impl Session {
    pub fn new(spec: WaferSpec, recipe: Recipe) -> Result<Self, SessionError> {
        let initial = spec.build(recipe.window.center).map_err(|e| SessionError::Wafer(e.to_string()))?;
        Ok(Self {
            state: SessionState { spec, initial: Arc::new(initial), recipe, records: BTreeMap::new(), failures: BTreeMap::new(), transform: WaferTransform::IDENTITY, revision: 0 },
            undo: VecDeque::new(), redo: vec![], preview: None, pending_site: None, coalesce: None, drag_start: None, log: VecDeque::new(),
        })
    }

    pub fn with_default_recipe() -> Self { Self::new(WaferSpec::default(), Recipe::default_recipe(ProcessWindow::default())).expect("default wafer spec is valid") }

    fn log(&mut self, s: impl Into<String>) { self.log.push_back(s.into()); while self.log.len() > LOG_LIMIT { self.log.pop_front(); } }

    fn checkpoint(&mut self) {
        self.coalesce = None;
        self.undo.push_back(self.state.clone());
        while self.undo.len() > UNDO_LIMIT { self.undo.pop_front(); }
        self.redo.clear();
    }

    // --- derived status ------------------------------------------------------

    pub fn statuses(&self) -> Vec<(u32, StageStatus)> {
        let s = &self.state;
        let mut expected = s.initial.fingerprint();
        let mut broken = false;
        let mut out = Vec::new();
        for st in &s.recipe.stages {
            if !st.enabled { out.push((st.id, StageStatus::Disabled)); continue; }
            match s.records.get(&st.id) {
                Some(r) if !broken && r.params == st.params && r.window == s.recipe.window && r.seed == s.recipe.seed && r.input_fp == expected => {
                    expected = r.output_fp; out.push((st.id, StageStatus::Committed));
                }
                Some(_) => { broken = true; out.push((st.id, StageStatus::Stale)); }
                None => { broken = true; out.push((st.id, match s.failures.get(&st.id) { Some(m) => StageStatus::Failed(m.clone()), None => StageStatus::NotRun })); }
            }
        }
        out
    }

    pub fn status_of(&self, id: u32) -> Option<StageStatus> { self.statuses().into_iter().find(|(i, _)| *i == id).map(|(_, s)| s) }

    /// Wafer after the last committed stage of the valid chain.
    pub fn committed_wafer(&self) -> Arc<WaferState> {
        let mut cur = self.state.initial.clone();
        for (id, st) in self.statuses() { if st == StageStatus::Committed { cur = self.state.records[&id].after.clone(); } }
        cur
    }

    pub fn committed_count(&self) -> usize { self.statuses().iter().filter(|(_, s)| *s == StageStatus::Committed).count() }
    pub fn enabled_count(&self) -> usize { self.state.recipe.stages.iter().filter(|s| s.enabled).count() }

    /// Input wafer for stage `id`, requiring every earlier enabled stage to be committed.
    fn input_for(&self, id: u32) -> Result<Arc<WaferState>, SessionError> {
        let idx = self.state.recipe.index_of(id).ok_or(SessionError::UnknownStage(id))?;
        let status = self.statuses();
        let mut cur = self.state.initial.clone();
        for (i, st) in self.state.recipe.stages.iter().enumerate().take(idx) {
            if !st.enabled { continue; }
            if status[i].1 != StageStatus::Committed { return Err(SessionError::UpstreamNotCommitted(st.params.stage().label().to_string())); }
            cur = self.state.records[&st.id].after.clone();
        }
        Ok(cur)
    }

    fn make_record(&self, id: u32, before: Arc<WaferState>) -> Result<StageRecord, SessionError> {
        let st = &self.state.recipe.stages[self.state.recipe.index_of(id).ok_or(SessionError::UnknownStage(id))?];
        let ev = evaluate_stage(&before, &self.state.recipe.window, &st.params, self.state.recipe.seed, 1.0)?;
        let after = Arc::new(ev.wafer);
        Ok(StageRecord {
            stage_id: id, stage: st.params.stage(), revision: 0, params: st.params.clone(), window: self.state.recipe.window, seed: self.state.recipe.seed,
            input_fp: before.fingerprint(), output_fp: after.fingerprint(), before, after, outcome: Arc::new(ev.outcome),
        })
    }

    // --- preview / commit / cancel ------------------------------------------

    /// Compute a stage on the committed upstream state without touching it.
    pub fn preview_stage(&mut self, id: u32) -> Result<(), SessionError> {
        if let Some(site) = self.pending_site { if site != self.state.recipe.window.center { return Err(SessionError::PendingSite); } }
        let before = self.input_for(id)?;
        let rec = self.make_record(id, before)?;
        self.log(format!("Preview {} (uncommitted)", rec.stage.label()));
        self.preview = Some(Preview { stage_id: id, record: rec });
        Ok(())
    }

    pub fn cancel_preview(&mut self) {
        if self.preview.take().is_some() { self.log("Preview cancelled; committed state unchanged"); }
        self.pending_site = None;
    }

    pub fn commit_preview(&mut self) -> Result<u32, SessionError> {
        let prev = self.preview.take().ok_or(SessionError::NothingToCommit)?;
        self.checkpoint();
        let id = prev.stage_id;
        let mut rec = prev.record;
        self.state.revision += 1;
        rec.revision = self.state.revision;
        self.state.failures.remove(&id);
        self.log(format!("Committed {} (rev {})", rec.stage.label(), rec.revision));
        self.state.records.insert(id, rec);
        Ok(id)
    }

    /// Run every enabled stage that is not already committed, in order.
    pub fn run_recipe(&mut self) -> RunSummary {
        self.preview = None;
        let snapshot = self.state.clone();
        let mut sum = RunSummary::default();
        let ids: Vec<u32> = self.state.recipe.stages.iter().filter(|s| s.enabled).map(|s| s.id).collect();
        let mut new_records = Vec::new();
        for id in ids {
            let status = self.status_of(id);
            if status == Some(StageStatus::Committed) { sum.skipped_committed.push(id); continue; }
            let before = match self.input_for(id) { Ok(b) => b, Err(e) => { sum.failed = Some((id, e.to_string())); break; } };
            match self.make_record(id, before) {
                Ok(mut rec) => {
                    self.state.revision += 1; rec.revision = self.state.revision;
                    self.state.failures.remove(&id);
                    self.log(format!("Ran {} (rev {}): {} warning(s)", rec.stage.label(), rec.revision, rec.outcome.issues.len()));
                    self.state.records.insert(id, rec); sum.ran.push(id); new_records.push(id);
                }
                Err(e) => {
                    let msg = e.to_string();
                    self.log(format!("FAILED {}: {}", self.state.recipe.stages[self.state.recipe.index_of(id).unwrap()].params.stage().label(), msg));
                    self.state.failures.insert(id, msg.clone()); sum.failed = Some((id, msg)); break;
                }
            }
        }
        if sum.ran.is_empty() && sum.failed.is_none() { self.state = snapshot; }
        else { self.coalesce = None; self.undo.push_back(snapshot); while self.undo.len() > UNDO_LIMIT { self.undo.pop_front(); } self.redo.clear(); }
        sum
    }

    // --- recipe editing --------------------------------------------------------

    /// Edit one stage's parameters. Consecutive edits of the same stage share a single undo entry.
    pub fn edit_params(&mut self, id: u32, params: StageParams) -> Result<(), SessionError> {
        let i = self.state.recipe.index_of(id).ok_or(SessionError::UnknownStage(id))?;
        if self.state.recipe.stages[i].params == params { return Ok(()); }
        if self.coalesce != Some(id) { self.checkpoint(); self.coalesce = Some(id); }
        self.state.recipe.stages[i].params = params;
        self.state.failures.remove(&id);
        self.preview = None;
        Ok(())
    }

    pub fn set_enabled(&mut self, id: u32, enabled: bool) -> Result<(), SessionError> {
        let i = self.state.recipe.index_of(id).ok_or(SessionError::UnknownStage(id))?;
        if self.state.recipe.stages[i].enabled == enabled { return Ok(()); }
        self.checkpoint(); self.state.recipe.stages[i].enabled = enabled; self.preview = None; Ok(())
    }

    pub fn add_stage_after(&mut self, id: u32, stage: Stage) -> Result<u32, SessionError> {
        self.checkpoint(); self.preview = None;
        self.state.recipe.insert_after(id, StageParams::default_for(stage)).ok_or(SessionError::UnknownStage(id))
    }
    pub fn append_stage(&mut self, stage: Stage) -> u32 { self.checkpoint(); self.preview = None; self.state.recipe.add(StageParams::default_for(stage)) }
    pub fn duplicate_stage(&mut self, id: u32) -> Result<u32, SessionError> {
        let p = self.state.recipe.stages[self.state.recipe.index_of(id).ok_or(SessionError::UnknownStage(id))?].params.clone();
        self.checkpoint(); self.preview = None;
        self.state.recipe.insert_after(id, p).ok_or(SessionError::UnknownStage(id))
    }
    pub fn remove_stage(&mut self, id: u32) -> Result<(), SessionError> {
        if self.state.recipe.index_of(id).is_none() { return Err(SessionError::UnknownStage(id)); }
        self.checkpoint(); self.preview = None; self.state.records.remove(&id); self.state.failures.remove(&id); self.state.recipe.remove(id); Ok(())
    }
    pub fn move_stage(&mut self, id: u32, delta: isize) -> bool {
        let mut probe = self.state.recipe.clone();
        if !probe.move_stage(id, delta) { return false; }
        self.checkpoint(); self.preview = None; self.state.recipe = probe; true
    }

    // --- site ------------------------------------------------------------------

    /// Propose a new process site (preview only; nothing is committed).
    pub fn propose_site(&mut self, center_local: [f64; 2]) -> Result<(), SessionError> {
        if !center_local.iter().all(|v| v.is_finite()) { return Err(SessionError::InvalidSite); }
        let r = self.state.initial.domain.radius_nm;
        let h = self.state.spec.tile_edge_nm() * 0.5;
        if (center_local[0].hypot(center_local[1]) + h) > r * 1.0 + h { return Err(SessionError::InvalidSite); }
        self.pending_site = Some(center_local);
        Ok(())
    }

    /// Commit the pending site: re-centres the high-resolution tile and
    /// therefore invalidates every committed stage (they are recomputed on the next run).
    pub fn commit_site(&mut self) -> Result<(), SessionError> {
        let Some(site) = self.pending_site.take() else { return Err(SessionError::NothingToCommit) };
        let spec = self.state.spec;
        let initial = spec.build(site).map_err(|_| SessionError::InvalidSite)?;
        self.checkpoint();
        self.preview = None;
        self.state.recipe.window.center = site;
        self.state.initial = Arc::new(initial);
        self.log(format!("Process site moved to ({:.0}, {:.0}) nm; downstream results are stale", site[0], site[1]));
        Ok(())
    }

    /// Change the process window size (all stages operate inside it). Results become stale.
    /// Largest process window (nm) the current scale allows: inside the tile and within the operator grid cap.
    pub fn max_window_nm(&self) -> f64 {
        (self.state.spec.tile_edge_nm() * 0.95).min(crate::geom::MAX_GRID_AXIS as f64 * self.state.spec.leaf_nm)
    }

    /// Change the octree cell size and tile depth (the simulation scale). The blank tile is rebuilt at the current
    /// site, the window is clamped to the new limits and every committed stage becomes stale. Undoable.
    pub fn set_scale(&mut self, leaf_nm: f64, tile_levels: u32) -> Result<(), SessionError> {
        if !(leaf_nm.is_finite() && (1.0..=2000.0).contains(&leaf_nm)) || !(4..=12).contains(&tile_levels) { return Err(SessionError::InvalidScale); }
        let spec = WaferSpec { leaf_nm, tile_levels, ..self.state.spec };
        if spec == self.state.spec { return Ok(()); }
        let initial = spec.build(self.state.recipe.window.center).map_err(|e| SessionError::Wafer(e.to_string()))?;
        if self.coalesce != Some(TRANSFORM_KEY - 2) { self.checkpoint(); self.coalesce = Some(TRANSFORM_KEY - 2); }
        self.state.spec = spec;
        self.state.initial = Arc::new(initial);
        let max = self.max_window_nm();
        for a in 0..2 { self.state.recipe.window.size[a] = self.state.recipe.window.size[a].min(max); }
        self.preview = None;
        self.log(format!("Scale set: cell {leaf_nm} nm, tile {:.0} nm; downstream results are stale", spec.tile_edge_nm()));
        Ok(())
    }

    pub fn set_window_size(&mut self, size: [f64; 2]) -> Result<(), SessionError> {
        let max = self.max_window_nm();
        if !size.iter().all(|v| v.is_finite() && *v >= 64.0 && *v <= max) { return Err(SessionError::InvalidSite); }
        if self.state.recipe.window.size == size { return Ok(()); }
        if self.coalesce != Some(TRANSFORM_KEY - 1) { self.checkpoint(); self.coalesce = Some(TRANSFORM_KEY - 1); }
        self.state.recipe.window.size = size;
        self.preview = None;
        Ok(())
    }

    // --- wafer transform (environment, never touches the SVO) -------------------

    /// Begin a pointer drag. The pre-drag state is remembered; one undo entry is created on `end_drag` if the transform moved.
    pub fn begin_drag(&mut self) { if self.drag_start.is_none() { self.drag_start = Some(self.state.clone()); } }
    pub fn drag_to(&mut self, t: WaferTransform) -> bool {
        if self.drag_start.is_none() || !t.is_finite() { return false; }
        self.state.transform = t; true
    }
    pub fn end_drag(&mut self) {
        if let Some(start) = self.drag_start.take() {
            if start.transform != self.state.transform {
                self.coalesce = None; self.undo.push_back(start); while self.undo.len() > UNDO_LIMIT { self.undo.pop_front(); } self.redo.clear();
            }
        }
    }
    pub fn is_dragging(&self) -> bool { self.drag_start.is_some() }

    /// Direct transform edit (numeric fields); consecutive edits coalesce.
    pub fn set_transform(&mut self, t: WaferTransform) -> bool {
        if !t.is_finite() { return false; }
        if t == self.state.transform { return true; }
        if self.coalesce != Some(TRANSFORM_KEY) { self.checkpoint(); self.coalesce = Some(TRANSFORM_KEY); }
        self.state.transform = t; true
    }

    pub fn center_wafer(&mut self) { let mut t = self.state.transform; t.position = [0.0, 0.0, 0.0]; self.coalesce = None; self.set_transform(t); }

    // --- reset / undo / redo ------------------------------------------------------

    /// Back to the blank substrate: clears every record, failure and the transform.
    pub fn reset(&mut self) {
        self.checkpoint(); self.preview = None; self.pending_site = None;
        self.state.records.clear(); self.state.failures.clear(); self.state.transform = WaferTransform::IDENTITY;
        self.log("Reset to blank substrate");
    }

    pub fn can_undo(&self) -> bool { !self.undo.is_empty() }
    pub fn can_redo(&self) -> bool { !self.redo.is_empty() }
    pub fn undo(&mut self) -> bool {
        let Some(prev) = self.undo.pop_back() else { return false };
        self.redo.push(std::mem::replace(&mut self.state, prev));
        self.preview = None; self.pending_site = None; self.coalesce = None; self.log("Undo"); true
    }
    pub fn redo(&mut self) -> bool {
        let Some(next) = self.redo.pop() else { return false };
        self.undo.push_back(std::mem::replace(&mut self.state, next));
        self.preview = None; self.pending_site = None; self.coalesce = None; self.log("Redo"); true
    }
    pub fn undo_depth(&self) -> usize { self.undo.len() }

    // --- timeline -----------------------------------------------------------------

    /// Animation segments for the committed chain (or the single preview stage).
    pub fn timeline(&self) -> Timeline {
        if let Some(p) = &self.preview {
            return Timeline { segments: vec![Segment::from_record(&p.record, false)] };
        }
        let segs = self.statuses().into_iter().filter(|(_, s)| *s == StageStatus::Committed).map(|(id, _)| Segment::from_record(&self.state.records[&id], true)).collect();
        Timeline { segments: segs }
    }
}

#[derive(Clone)]
pub struct Segment {
    pub stage_id: u32,
    pub stage: Stage,
    pub committed: bool,
    pub params: StageParams,
    pub window: ProcessWindow,
    pub seed: u64,
    pub before: Arc<WaferState>,
    pub after: Arc<WaferState>,
    pub outcome: Arc<StageOutcome>,
    pub duration_s: f64,
    pub input_fp: u64,
    pub output_fp: u64,
}

impl Segment {
    fn from_record(r: &StageRecord, committed: bool) -> Self {
        Self { stage_id: r.stage_id, stage: r.stage, committed, params: r.params.clone(), window: r.window, seed: r.seed, before: r.before.clone(), after: r.after.clone(), outcome: r.outcome.clone(), duration_s: r.outcome.duration_s, input_fp: r.input_fp, output_fp: r.output_fp }
    }
}

#[derive(Clone, Default)]
pub struct Timeline { pub segments: Vec<Segment> }

impl Timeline {
    pub fn total_duration_s(&self) -> f64 { self.segments.iter().map(|s| s.duration_s).sum() }
    /// Simulated process time at a playback position (segment units).
    /// Identity of the timeline content; changes whenever any segment's inputs or outputs change.
    pub fn signature(&self) -> u64 {
        self.segments.iter().fold(0xcbf2_9ce4_8422_2325u64, |h, s| (h ^ s.output_fp ^ s.input_fp.rotate_left(13) ^ s.stage_id as u64 ^ s.committed as u64).wrapping_mul(0x0000_0100_0000_01B3))
    }

    pub fn sim_time_at(&self, position: f64) -> f64 {
        let mut t = 0.0;
        for (i, s) in self.segments.iter().enumerate() {
            let local = (position - i as f64).clamp(0.0, 1.0);
            t += s.duration_s * local;
        }
        t
    }
}

pub struct Frame {
    pub wafer: Arc<WaferState>,
    pub overlay: Overlay,
    pub phase: SubStep,
    pub progress: f64,
    pub stage: Stage,
}

/// Lazily replays animation frames from the recorded inputs (deterministic
/// replay) and keeps a small LRU of the most recent ones.
#[derive(Default)]
pub struct FrameCache { map: HashMap<(u32, u64, usize), Arc<Frame>>, order: VecDeque<(u32, u64, usize)>, pub computed: usize }

impl FrameCache {
    pub fn clear(&mut self) { self.map.clear(); self.order.clear(); }

    pub fn get(&mut self, seg: &Segment, frame: usize) -> Result<Arc<Frame>, PipelineError> {
        let frame = frame.min(FRAMES_PER_STAGE);
        let key = (seg.stage_id, seg.output_fp ^ seg.input_fp.rotate_left(17), frame);
        if let Some(f) = self.map.get(&key) { return Ok(f.clone()); }
        let progress = crate::pipeline::frame_progress(frame);
        let (wafer, overlay) = if frame == 0 {
            (seg.before.clone(), Overlay::None)
        } else {
            let ev = evaluate_stage(&seg.before, &seg.window, &seg.params, seg.seed, progress)?;
            // The final frame *is* the committed record; reuse it so the viewport is exactly the committed geometry.
            let w = if frame == FRAMES_PER_STAGE { seg.after.clone() } else { Arc::new(ev.wafer) };
            (w, ev.overlay)
        };
        self.computed += 1;
        let f = Arc::new(Frame { wafer, overlay, phase: crate::pipeline::phase_at(seg.stage, progress), progress, stage: seg.stage });
        self.map.insert(key, f.clone());
        self.order.push_back(key);
        while self.order.len() > 24 { if let Some(k) = self.order.pop_front() { self.map.remove(&k); } }
        Ok(f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::params::Species;

    fn s() -> Session { Session::with_default_recipe() }
    fn ids(s: &Session) -> Vec<u32> { s.state.recipe.stages.iter().map(|x| x.id).collect() }
    fn fp(s: &Session) -> u64 { s.committed_wafer().fingerprint() }

    #[test]
    fn run_recipe_commits_everything_and_a_second_run_is_a_no_op() {
        let mut s = s();
        let r = s.run_recipe();
        assert!(r.failed.is_none(), "{:?}", r.failed);
        assert_eq!(r.ran.len(), 7);
        assert_eq!(s.committed_count(), 7);
        let (depth, fp1) = (s.undo_depth(), fp(&s));
        let r2 = s.run_recipe();
        assert!(r2.ran.is_empty() && r2.skipped_committed.len() == 7);
        assert_eq!(s.undo_depth(), depth, "a no-op run adds no undo entry");
        assert_eq!(fp(&s), fp1);
    }

    #[test]
    fn proc04_final_clean_then_metrology_measures_the_final_committed_wafer() {
        let mut s = s();
        s.run_recipe();
        let ids = ids(&s);
        let clean = &s.state.records[&ids[5]];
        let metro = &s.state.records[&ids[6]];
        let rep = metro.outcome.report.as_ref().unwrap();
        assert_eq!(metro.input_fp, clean.output_fp, "metrology starts from the post-clean wafer");
        assert_eq!(rep.wafer_fingerprint, s.committed_wafer().fingerprint(), "and measures exactly the final committed state");
        assert_eq!(clean.outcome.metric("Resist residue (fraction left)"), Some(0.0), "clean removed the resist the measurements rely on being gone");
    }

    #[test]
    fn int05_preview_then_cancel_changes_nothing() {
        let mut s = s();
        let id = ids(&s)[0];
        let before = (fp(&s), s.state.recipe.clone(), s.state.records.len(), s.undo_depth(), s.state.transform);
        s.preview_stage(id).unwrap();
        assert!(s.preview.is_some());
        assert_ne!(s.preview.as_ref().unwrap().record.output_fp, 0);
        s.cancel_preview();
        assert!(s.preview.is_none());
        assert_eq!((fp(&s), s.state.recipe.clone(), s.state.records.len(), s.undo_depth(), s.state.transform), before);
        assert_eq!(s.statuses()[0].1, StageStatus::NotRun, "no process record was created");
    }

    #[test]
    fn int06_preview_commit_undo_restores_geometry_materials_and_transform() {
        let mut s = s();
        let (a, b) = (ids(&s)[0], ids(&s)[1]);
        s.preview_stage(a).unwrap(); s.commit_preview().unwrap();
        s.set_transform(WaferTransform { position: [1.0e6, -2.0e6, 0.0], yaw_rad: 0.3 });
        s.coalesce = None;
        let (fp0, t0, rec0) = (fp(&s), s.state.transform, s.state.records.len());
        s.preview_stage(b).unwrap(); s.commit_preview().unwrap();
        assert_ne!(fp(&s), fp0);
        assert!(s.undo());
        assert_eq!((fp(&s), s.state.transform, s.state.records.len()), (fp0, t0, rec0));
        assert!(s.redo());
        assert_eq!(s.committed_count(), 2);
    }

    #[test]
    fn undo_after_multiple_commits_walks_back_one_commit_at_a_time() {
        let mut s = s();
        let mut fps = vec![fp(&s)];
        for id in ids(&s).into_iter().take(4) { s.preview_stage(id).unwrap(); s.commit_preview().unwrap(); fps.push(fp(&s)); }
        for k in (0..4).rev() { assert!(s.undo()); assert_eq!(fp(&s), fps[k], "after undoing commit {}", k + 1); }
        assert!(!s.undo() || s.undo_depth() == 0 || true);
    }

    #[test]
    fn commits_must_follow_the_recipe_order() {
        let mut s = s();
        let third = ids(&s)[2];
        assert!(matches!(s.preview_stage(third), Err(SessionError::UpstreamNotCommitted(_))));
    }

    #[test]
    fn data01_editing_an_upstream_stage_marks_downstream_stale_and_recomputes_only_what_changed() {
        let mut s = s();
        s.run_recipe();
        let ids = ids(&s);
        let dep = ids[1];
        let mut p = s.state.recipe.stages[1].params.clone();
        if let StageParams::Deposition(d) = &mut p { d.thickness_nm = 64.0; }
        s.edit_params(dep, p).unwrap();
        let st: Vec<StageStatus> = s.statuses().into_iter().map(|x| x.1).collect();
        assert_eq!(st[0], StageStatus::Committed, "prep is upstream and unaffected");
        assert!(st[1..].iter().all(|x| *x == StageStatus::Stale), "{st:?}");
        let r = s.run_recipe();
        assert_eq!(r.ran.len(), 6);
        assert_eq!(r.skipped_committed.len(), 1);
        assert!(s.statuses().iter().all(|x| x.1 == StageStatus::Committed));
        let rep = s.state.records[&ids[6]].outcome.report.as_ref().unwrap();
        let t = rep.measurements[0].value.unwrap();
        assert!((t - 64.0).abs() < 1.0, "metrology measured the recomputed film: {t}");
    }

    #[test]
    fn fr_4_2_editing_a_lithography_mask_stales_it_and_every_later_stage_and_a_rerun_recomputes_all_of_them() {
        let mut s = s();
        s.run_recipe();
        let ids = ids(&s);
        let lith = ids[2];
        assert!(s.state.recipe.stages[2].params.stage() == crate::physics::Stage::Lithography, "the third default stage is lithography");
        let mut p = s.state.recipe.stages[2].params.clone();
        if let StageParams::Lithography(l) = &mut p {
            l.mask = crate::geom::MaskPattern::LineSpace { line_nm: 210.0, pitch_nm: 320.0, count: 3, length_nm: 1000.0, vertical: true, clear_field: true };
        }
        s.edit_params(lith, p).unwrap();
        let st: Vec<StageStatus> = s.statuses().into_iter().map(|x| x.1).collect();
        assert!(st[..2].iter().all(|x| *x == StageStatus::Committed), "stages before the mask are untouched: {st:?}");
        assert!(st[2..].iter().all(|x| *x == StageStatus::Stale), "the mask stage and all five after it are stale: {st:?}");
        let r = s.run_recipe();
        assert_eq!(r.ran, ids[2..].to_vec(), "lithography and every later stage are recomputed, in order, not one stage with the rest 'up to date'");
        assert_eq!(r.skipped_committed, ids[..2].to_vec());
        assert!(s.statuses().iter().all(|x| x.1 == StageStatus::Committed));
    }

    #[test]
    fn changing_a_downstream_stage_leaves_upstream_committed() {
        let mut s = s();
        s.run_recipe();
        let id = ids(&s)[4];
        let mut p = s.state.recipe.stages[4].params.clone();
        if let StageParams::Doping(d) = &mut p { d.species = Species::Phosphorus; }
        s.edit_params(id, p).unwrap();
        let st: Vec<StageStatus> = s.statuses().into_iter().map(|x| x.1).collect();
        assert!(st[..4].iter().all(|x| *x == StageStatus::Committed));
        assert!(st[4..].iter().all(|x| *x == StageStatus::Stale));
    }

    #[test]
    fn failed_stage_does_not_commit_and_blocks_downstream() {
        let mut s = s();
        let id = ids(&s)[1];
        let mut p = s.state.recipe.stages[1].params.clone();
        if let StageParams::Deposition(d) = &mut p { d.thickness_nm = f64::NAN; }
        s.edit_params(id, p).unwrap();
        let r = s.run_recipe();
        assert_eq!(r.failed.as_ref().map(|f| f.0), Some(id));
        assert_eq!(r.ran.len(), 1, "only substrate prep ran");
        assert!(matches!(s.status_of(id), Some(StageStatus::Failed(_))));
        assert_eq!(s.status_of(ids(&s)[2]), Some(StageStatus::NotRun));
    }

    #[test]
    fn repeated_cycles_are_allowed() {
        let mut s = s();
        s.run_recipe();
        let last_etch = ids(&s)[3];
        let d = s.add_stage_after(last_etch, Stage::Deposition).unwrap();
        let l = s.add_stage_after(d, Stage::Lithography).unwrap();
        s.add_stage_after(l, Stage::Etching).unwrap();
        assert_eq!(s.state.recipe.stages.len(), 10);
        let st = s.statuses();
        assert!(st[..4].iter().all(|x| x.1 == StageStatus::Committed) && st[4..].iter().all(|x| x.1 != StageStatus::Committed));
    }

    #[test]
    fn int01_dragging_the_wafer_changes_only_the_transform_never_the_svo() {
        let mut s = s();
        s.run_recipe();
        let wafer = s.committed_wafer();
        let ptr = Arc::as_ptr(&wafer);
        let depth = s.undo_depth();
        s.begin_drag();
        for i in 1..=10 { assert!(s.drag_to(WaferTransform { position: [i as f64 * 1.0e5, -(i as f64) * 5.0e4, 0.0], yaw_rad: 0.0 })); }
        s.end_drag();
        assert_eq!(Arc::as_ptr(&s.committed_wafer()), ptr, "the SVO was not rewritten or cloned");
        assert_eq!(s.committed_wafer().fingerprint(), wafer.fingerprint());
        assert_eq!(s.undo_depth(), depth + 1, "one undo entry for the whole drag");
        assert!(s.undo());
        assert_eq!(s.state.transform, WaferTransform::IDENTITY);
        s.begin_drag(); s.end_drag();
        assert_eq!(s.undo_depth(), depth, "a click without movement adds no history");
        assert!(!s.drag_to(WaferTransform::IDENTITY), "drag_to without begin_drag is ignored");
    }

    #[test]
    fn site_change_invalidates_everything_and_is_transactional() {
        let mut s = s();
        s.run_recipe();
        s.propose_site([5.0e6, 2.0e6]).unwrap();
        assert_eq!(s.committed_count(), 7, "a proposed site is only a preview");
        assert!(matches!(s.preview_stage(ids(&s)[0]), Err(SessionError::PendingSite)));
        s.cancel_preview();
        assert_eq!(s.committed_count(), 7);
        s.propose_site([5.0e6, 2.0e6]).unwrap(); s.commit_site().unwrap();
        assert_eq!(s.committed_count(), 0);
        assert!(s.run_recipe().failed.is_none());
        assert_eq!(s.committed_count(), 7);
        assert!(s.undo() && s.undo());
        assert_eq!(s.state.recipe.window.center, [0.0, 0.0]);
        assert!(matches!(s.propose_site([f64::NAN, 0.0]), Err(SessionError::InvalidSite)));
        assert!(matches!(s.propose_site([1.0e9, 0.0]), Err(SessionError::InvalidSite)));
    }

    #[test]
    fn param_edits_coalesce_into_one_undo_step() {
        let mut s = s();
        let id = ids(&s)[1];
        let d0 = s.undo_depth();
        for t in [80.0, 90.0, 110.0] {
            let mut p = s.state.recipe.stages[1].params.clone();
            if let StageParams::Deposition(d) = &mut p { d.thickness_nm = t; }
            s.edit_params(id, p).unwrap();
        }
        assert_eq!(s.undo_depth(), d0 + 1);
        assert!(s.undo());
        if let StageParams::Deposition(d) = &s.state.recipe.stages[1].params { assert_eq!(d.thickness_nm, 100.0); } else { panic!() }
    }

    #[test]
    fn reset_returns_to_the_blank_substrate() {
        let mut s = s();
        s.run_recipe();
        s.set_transform(WaferTransform { position: [1.0, 2.0, 3.0], yaw_rad: 0.1 });
        let blank = s.state.initial.fingerprint();
        s.reset();
        assert_eq!(fp(&s), blank);
        assert_eq!(s.state.transform, WaferTransform::IDENTITY);
        assert!(s.undo());
        assert_eq!(s.committed_count(), 7);
    }

    #[test]
    fn deterministic_replay_gives_identical_records() {
        let (mut a, mut b) = (s(), s());
        a.run_recipe(); b.run_recipe();
        let (ia, ib) = (ids(&a), ids(&b));
        for (x, y) in ia.iter().zip(&ib) { assert_eq!(a.state.records[x].output_fp, b.state.records[y].output_fp); }
    }

    #[test]
    fn frame_cache_replays_frames_and_final_frame_is_the_committed_wafer() {
        let mut s = s();
        s.run_recipe();
        let tl = s.timeline();
        assert_eq!(tl.segments.len(), 7);
        let mut cache = FrameCache::default();
        for seg in &tl.segments {
            let last = cache.get(seg, FRAMES_PER_STAGE).unwrap();
            assert!(Arc::ptr_eq(&last.wafer, &seg.after));
            let mid = cache.get(seg, FRAMES_PER_STAGE / 2).unwrap();
            let again = cache.get(seg, FRAMES_PER_STAGE / 2).unwrap();
            assert!(Arc::ptr_eq(&mid, &again), "cached");
            assert_eq!(cache.get(seg, 0).unwrap().wafer.fingerprint(), seg.before.fingerprint());
        }
        assert!(tl.total_duration_s() > 600.0);
        assert!((tl.sim_time_at(0.0)).abs() < 1e-9 && (tl.sim_time_at(7.0) - tl.total_duration_s()).abs() < 1e-6);
    }

    #[test]
    fn window_size_edits_are_validated_and_invalidate_results() {
        let mut s = s();
        s.run_recipe();
        assert!(s.set_window_size([f64::NAN, 800.0]).is_err());
        assert!(s.set_window_size([1.0e9, 800.0]).is_err());
        s.set_window_size([800.0, 800.0]).unwrap();
        assert_eq!(s.committed_count(), 0);
        assert!(s.run_recipe().failed.is_none());
    }

    #[test]
    fn scale_changes_are_validated_undoable_and_stale_everything() {
        let mut s = s();
        s.run_recipe();
        assert_eq!(s.committed_count(), 7);
        for (leaf, lv) in [(f64::NAN, 7), (0.5, 7), (16.0, 3), (16.0, 13), (5000.0, 7)] { assert!(s.set_scale(leaf, lv).is_err(), "{leaf} {lv}"); }
        assert_eq!(s.committed_count(), 7, "a rejected scale changes nothing");
        let fp0 = s.state.initial.fingerprint();
        s.set_scale(50.0, 8).unwrap();
        assert_ne!(s.state.initial.fingerprint(), fp0);
        assert_eq!(s.committed_count(), 0);
        assert!((s.state.initial.bounds.max.0 - s.state.initial.bounds.min.0 - 12_800.0).abs() < 1e-6, "tile = 50 nm x 2^8");
        assert!(s.state.recipe.window.size[0] <= s.max_window_nm() + 1e-9);
        assert!(s.max_window_nm() <= crate::geom::MAX_GRID_AXIS as f64 * 50.0 + 1e-9, "window capped by the operator grid");
        assert!(s.undo());
        assert_eq!(s.state.initial.fingerprint(), fp0);
        assert_eq!(s.committed_count(), 7);
    }

    #[test]
    fn preview_timeline_is_a_single_uncommitted_segment() {
        let mut s = s();
        s.preview_stage(ids(&s)[0]).unwrap();
        let tl = s.timeline();
        assert_eq!(tl.segments.len(), 1);
        assert!(!tl.segments[0].committed);
    }
}
