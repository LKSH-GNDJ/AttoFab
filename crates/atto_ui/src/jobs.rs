// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Jobs: the interface never blocks, and cancelling discards partial work
//! (PRD G-06, G-24).
//!
//! There is no threading here and no engine — the caller drives the work and the
//! queue only records what state each job is in. That is the point: a panel asks
//! the queue what is running and draws it, instead of waiting on a computation.
//! A job that is cancelled leaves nothing behind, because a half-applied stage is
//! worse than no stage at all.

use crate::{Progress, StageId};

#[derive(Debug, Clone, PartialEq)]
/// Where a job has got to.
pub enum JobState {
    Queued,
    Running(Progress),
    Done,
    Failed(String),
    Cancelled,
}

impl JobState {
    /// Whether this state is final.
    pub fn is_finished(&self) -> bool {
        matches!(self, JobState::Done | JobState::Failed(_) | JobState::Cancelled)
    }

    /// The one-word status a panel draws.
    pub fn word(&self) -> &'static str {
        match self {
            JobState::Queued => "Queued",
            JobState::Running(_) => "Running",
            JobState::Done => "Done",
            JobState::Failed(_) => "Failed",
            JobState::Cancelled => "Cancelled",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
/// One submitted stage.
pub struct Job {
    pub stage: StageId,
    pub state: JobState,
}

/// At most one job runs at a time. Starting a second one queues it rather than
/// interleaving two stages over the same wafer.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct JobQueue {
    jobs: Vec<Job>,
}

impl JobQueue {
    /// An empty queue.
    pub fn new() -> Self {
        JobQueue { jobs: Vec::new() }
    }

    /// How many jobs have been submitted.
    pub fn len(&self) -> usize {
        self.jobs.len()
    }

    /// Whether anything has been submitted.
    pub fn is_empty(&self) -> bool {
        self.jobs.is_empty()
    }

    /// What is running, if anything. A panel draws this; it never waits on it.
    pub fn running(&self) -> Option<&Job> {
        self.jobs.iter().find(|j| matches!(j.state, JobState::Running(_)))
    }

    /// Submit a stage. Returns the index for later polling.
    pub fn submit(&mut self, stage: StageId, total_steps: usize) -> usize {
        let idx = self.jobs.len();
        self.jobs.push(Job { stage, state: JobState::Queued });
        // nothing else running: start immediately, but do not run it
        if self.running().is_none() {
            self.jobs[idx].state = JobState::Running(Progress::new(stage.label(), total_steps));
        }
        idx
    }

    /// Advance a running job by one step.
    pub fn step(&mut self, idx: usize) -> bool {
        match self.jobs.get_mut(idx) {
            Some(job) => match &mut job.state {
                JobState::Running(p) => {
                    p.step();
                    true
                }
                _ => false,
            },
            None => false,
        }
    }

    /// Mark a job done and let the queue move on.
    pub fn finish(&mut self, idx: usize) {
        if let Some(job) = self.jobs.get_mut(idx) {
            if matches!(job.state, JobState::Running(_)) {
                job.state = JobState::Done;
            }
        }
        self.promote();
    }

    /// Mark a job failed, keeping the reason.
    pub fn fail(&mut self, idx: usize, why: &str) {
        if let Some(job) = self.jobs.get_mut(idx) {
            if matches!(job.state, JobState::Running(_)) {
                job.state = JobState::Failed(why.to_string());
            }
        }
        self.promote();
    }

    /// Cancel a job. A cancelled job reports nothing: no partial result is kept.
    pub fn cancel(&mut self, idx: usize) -> bool {
        match self.jobs.get_mut(idx) {
            Some(job) if !job.state.is_finished() => {
                job.state = JobState::Cancelled;
                self.promote();
                true
            }
            _ => false,
        }
    }

    /// The state of a submitted job.
    pub fn state(&self, idx: usize) -> Option<&JobState> {
        self.jobs.get(idx).map(|j| &j.state)
    }

    /// The next queued job becomes the running one.
    fn promote(&mut self) {
        if self.running().is_some() {
            return;
        }
        if let Some(next) = self.jobs.iter_mut().find(|j| matches!(j.state, JobState::Queued)) {
            next.state = JobState::Running(Progress::new(next.stage.label(), 1));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn submitting_does_not_run_the_work() {
        // G-06: the call returns immediately; the job is only marked running
        let mut q = JobQueue::new();
        let i = q.submit(StageId::Etching, 8);
        assert!(matches!(q.state(i), Some(JobState::Running(_))));
        match q.state(i) {
            Some(JobState::Running(p)) => assert_eq!(p.done, 0, "no step has happened yet"),
            other => panic!("{other:?}"),
        }
        assert!(q.running().is_some());
    }

    #[test]
    fn a_second_job_waits_for_the_first() {
        let mut q = JobQueue::new();
        let a = q.submit(StageId::Etching, 4);
        let b = q.submit(StageId::Doping, 4);
        assert!(matches!(q.state(b), Some(JobState::Queued)), "one job at a time");
        q.finish(a);
        assert!(matches!(q.state(b), Some(JobState::Running(_))), "the queue promotes");
    }

    #[test]
    fn cancelling_discards_the_partial_work() {
        // G-24: a cancelled job reports nothing, not a half-applied stage
        let mut q = JobQueue::new();
        let i = q.submit(StageId::Oxidation, 10);
        for _ in 0..4 {
            q.step(i);
        }
        assert!(q.cancel(i));
        assert_eq!(q.state(i), Some(&JobState::Cancelled));
        assert!(q.state(i).unwrap().is_finished());
        assert!(!q.cancel(i), "cancelling twice is not an error, just a no-op");
        assert!(q.running().is_none(), "nothing runs after a cancel");
    }

    #[test]
    fn a_failure_keeps_its_message_and_frees_the_queue() {
        let mut q = JobQueue::new();
        let a = q.submit(StageId::Etching, 2);
        let b = q.submit(StageId::Anneal, 2);
        q.fail(a, "the plasma model rejected the pressure");
        match q.state(a) {
            Some(JobState::Failed(why)) => assert!(why.contains("pressure"), "{why}"),
            other => panic!("{other:?}"),
        }
        assert!(matches!(q.state(b), Some(JobState::Running(_))));
    }

    #[test]
    fn a_finished_job_cannot_be_cancelled_or_advanced() {
        let mut q = JobQueue::new();
        let i = q.submit(StageId::Metrology, 2);
        q.finish(i);
        assert_eq!(q.state(i), Some(&JobState::Done));
        assert!(!q.cancel(i));
        assert!(!q.step(i));
    }

    #[test]
    fn stepping_advances_the_progress() {
        let mut q = JobQueue::new();
        let i = q.submit(StageId::Etching, 4);
        assert!(q.step(i));
        assert!(q.step(i));
        match q.state(i) {
            Some(JobState::Running(p)) => {
                assert_eq!(p.done, 2);
                assert!((p.fraction() - 0.5).abs() < 1e-12);
            }
            other => panic!("{other:?}"),
        }
    }
}
