//! Reactive dependency graph management + asynchronous non-blocking
//! execution.
//!
//! The six fabrication stages are not independent: an etch's thermal
//! profile can perturb an upstream deposition's stress state, and doping
//! diffusion is sensitive to thermal budget from later anneals. Rather than
//! re-running the entire linear pipeline whenever one parameter changes,
//! the engine tracks a directed dependency graph between stages and only
//! marks the affected downstream (and, where wired, upstream) stages dirty.
//!
//! Heavy stage computations are dispatched to a background worker thread via
//! a crossbeam channel so the UI thread (egui, 60+ FPS target) is never
//! blocked waiting on a Monte Carlo or SVO raymarch job.

use crate::physics::Stage;
use crossbeam::channel::{unbounded, Receiver, Sender};
use std::collections::{HashMap, HashSet};

/// A job dispatched to the background worker: which stage to recompute, and
/// an opaque serialized parameter payload (JSON) for that stage's `run()`.
#[derive(Debug, Clone)]
pub struct StageJob {
    pub stage: Stage,
    pub params_json: String,
}

/// The result of a completed job: the stage that ran, and its serialized
/// output (or an error message).
#[derive(Debug, Clone)]
pub struct StageResult {
    pub stage: Stage,
    pub output_json: Result<String, String>,
}

/// Manages which stages are "dirty" (need recomputation) based on a
/// directed dependency graph, and owns the channels used to hand work off
/// to a background thread pool.
pub struct Engine {
    /// stage -> stages that depend on it (i.e. must be marked dirty when it changes)
    downstream: HashMap<Stage, Vec<Stage>>,
    dirty: HashSet<Stage>,
    job_tx: Sender<StageJob>,
    job_rx: Receiver<StageJob>,
    result_tx: Sender<StageResult>,
    result_rx: Receiver<StageResult>,
}

impl Engine {
    /// Build an engine wired with the canonical linear pipeline dependency
    /// graph (each stage feeds the next), plus the two known cross-stage
    /// feedback edges called out in the PRD: Etching -> Deposition (etch
    /// heat perturbs film stress) and Etching -> Doping (etch heat drives
    /// dopant diffusion).
    pub fn new_default() -> Self {
        let mut downstream: HashMap<Stage, Vec<Stage>> = HashMap::new();
        downstream.insert(Stage::SubstratePrep, vec![Stage::Deposition]);
        downstream.insert(Stage::Deposition, vec![Stage::Lithography]);
        downstream.insert(Stage::Lithography, vec![Stage::Etching]);
        downstream.insert(Stage::Etching, vec![Stage::Doping, Stage::Deposition]); // feedback edge
        downstream.insert(Stage::Doping, vec![Stage::FinalClean]);
        downstream.insert(Stage::Cmp, vec![Stage::FinalClean, Stage::Metrology]);
downstream.insert(Stage::FinalClean, vec![Stage::Metrology]);
        downstream.insert(Stage::Metrology, vec![]);

        let (job_tx, job_rx) = unbounded();
        let (result_tx, result_rx) = unbounded();

        Self {
            downstream,
            dirty: HashSet::new(),
            job_tx,
            job_rx,
            result_tx,
            result_rx,
        }
    }

    /// Mark a stage dirty because its parameters changed, and propagate
    /// dirtiness to every downstream-dependent stage (transitively, with
    /// cycle protection since Etching<->Deposition forms a feedback loop).
    pub fn mark_dirty(&mut self, stage: Stage) {
        let mut stack = vec![stage];
        let mut visited: HashSet<Stage> = HashSet::new();
        while let Some(s) = stack.pop() {
            if !visited.insert(s) {
                continue;
            }
            self.dirty.insert(s);
            if let Some(deps) = self.downstream.get(&s) {
                for &d in deps {
                    stack.push(d);
                }
            }
        }
    }

    pub fn dirty_stages(&self) -> Vec<Stage> {
        self.dirty.iter().copied().collect()
    }

    pub fn is_dirty(&self, stage: Stage) -> bool {
        self.dirty.contains(&stage)
    }

    /// Submit a recompute job for a stage to the background worker channel
    /// and clear its dirty flag optimistically (it will be re-marked dirty
    /// if a later parameter change arrives before the result comes back).
    pub fn submit(&mut self, job: StageJob) {
        self.dirty.remove(&job.stage);
        // A production build spawns a small pool of OS threads that `recv()`
        // from `job_rx` in a loop, run the matching `physics::<stage>::run`,
        // and `send` a `StageResult` back on `result_tx`. That pool lives in
        // the platform shell (app_desktop/app_web) so it can be shut down
        // cleanly with the window; sim_core only owns the channel plumbing.
        let _ = self.job_tx.send(job);
    }

    pub fn job_receiver(&self) -> Receiver<StageJob> {
        self.job_rx.clone()
    }

    pub fn result_sender(&self) -> Sender<StageResult> {
        self.result_tx.clone()
    }

    /// Non-blocking drain of any results that have arrived, for the UI
    /// thread to call once per frame.
    pub fn poll_results(&self) -> Vec<StageResult> {
        self.result_rx.try_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dirtying_lithography_propagates_downstream() {
        let mut engine = Engine::new_default();
        engine.mark_dirty(Stage::Lithography);
        assert!(engine.is_dirty(Stage::Lithography));
        assert!(engine.is_dirty(Stage::Etching));
        assert!(engine.is_dirty(Stage::Doping));
        assert!(engine.is_dirty(Stage::Metrology));
        // Upstream, independent stages should not be touched.
        assert!(!engine.is_dirty(Stage::SubstratePrep));
    }

    #[test]
    fn etch_deposition_feedback_loop_does_not_hang() {
        let mut engine = Engine::new_default();
        // This would infinite-loop without cycle protection, since
        // Etching -> Deposition -> Lithography -> Etching is a cycle.
        engine.mark_dirty(Stage::Etching);
        assert!(engine.is_dirty(Stage::Deposition));
        assert!(engine.is_dirty(Stage::Lithography));
    }

    #[test]
    fn submitted_job_is_received_by_worker_channel() {
        let mut engine = Engine::new_default();
        engine.mark_dirty(Stage::Doping);
        let rx = engine.job_receiver();
        engine.submit(StageJob { stage: Stage::Doping, params_json: "{}".into() });
        assert!(!engine.is_dirty(Stage::Doping));
        let job = rx.try_recv().unwrap();
        assert!(matches!(job.stage, Stage::Doping));
    }
}

/// A running pool of worker threads draining the engine's job channel.
///
/// R-23 asked for this: the engine owned the channels but nothing drained them, so a
/// submitted job went nowhere. The pool takes the handler as an argument because the
/// physics belongs to the caller; what lives here is the threading, which is the part
/// that must be shut down cleanly.
///
/// Dropping the pool signals every worker to stop and joins it.
pub struct WorkerPool {
    stop: std::sync::Arc<std::sync::atomic::AtomicBool>,
    handles: Vec<std::thread::JoinHandle<()>>,
}

impl WorkerPool {
    /// How many workers are running.
    pub fn workers(&self) -> usize { self.handles.len() }
}

impl Drop for WorkerPool {
    fn drop(&mut self) {
        self.stop.store(true, std::sync::atomic::Ordering::Relaxed);
        for h in self.handles.drain(..) {
            let _ = h.join();
        }
    }
}

impl Engine {
    /// Start `workers` threads that drain the job channel and answer on the result channel.
    ///
    /// `run` is the physics: one job in, one result out. It must be `Send + Sync` because
    /// every worker calls it. Returns the pool, which stops the workers when dropped.
    pub fn start_pool<F>(&self, workers: usize, run: std::sync::Arc<F>) -> WorkerPool
    where
        F: Fn(StageJob) -> StageResult + Send + Sync + 'static,
    {
        let stop_flag = std::sync::Arc::new(std::sync::atomic::AtomicBool::new(false));
        let mut handles = Vec::new();
        for _ in 0..workers.max(1) {
            let jobs = self.job_rx.clone();
            let results = self.result_tx.clone();
            let run = std::sync::Arc::clone(&run);
            let stop_flag = std::sync::Arc::clone(&stop_flag);
            handles.push(std::thread::spawn(move || loop {
                if stop_flag.load(std::sync::atomic::Ordering::Relaxed) {
                    break;
                }
                match jobs.recv_timeout(std::time::Duration::from_millis(20)) {
                    Ok(j) => {
                        let _ = results.send(run(j));
                    }
                    Err(_) => continue,
                }
            }));
        }
        WorkerPool { stop: stop_flag, handles }
    }
}

#[cfg(test)]
mod pool_tests {
    use super::*;
    use std::sync::atomic::{AtomicUsize, Ordering};

    fn job(stage: Stage) -> StageJob {
        StageJob { stage, params_json: String::new() }
    }

    #[test]
    fn the_pool_drains_every_submitted_job() {
        // the defect R-23 names: jobs went into a channel nothing read
        let mut e = Engine::new_default();
        let seen = std::sync::Arc::new(AtomicUsize::new(0));
        let s2 = std::sync::Arc::clone(&seen);
        let _pool = e.start_pool(3, std::sync::Arc::new(move |j: StageJob| {
            s2.fetch_add(1, Ordering::SeqCst);
            StageResult { stage: j.stage, output_json: Ok(String::new()) }
        }));
        for st in [Stage::Lithography, Stage::Etching, Stage::Doping] {
            e.submit(job(st));
        }
        let mut got = 0;
        for _ in 0..200 {
            got += e.poll_results().len();
            if got == 3 { break; }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }
        assert_eq!(got, 3, "every job must come back");
        assert_eq!(seen.load(Ordering::SeqCst), 3);
    }

    #[test]
    fn dropping_the_pool_stops_its_workers() {
        // a pool that outlives its engine is a thread leak the window would inherit
        let e = Engine::new_default();
        let pool = e.start_pool(4, std::sync::Arc::new(|j: StageJob| StageResult {
            stage: j.stage, output_json: Ok(String::new()),
        }));
        assert_eq!(pool.workers(), 4);
        drop(pool); // must join, not hang
    }
}
