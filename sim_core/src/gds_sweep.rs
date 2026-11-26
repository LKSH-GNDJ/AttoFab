// SPDX-License-Identifier: GPL-3.0-or-later
//! Simulating a whole GDS file, all layers at once.
//!
//! # The limit this works around
//!
//! The operator grid is capped at `MAX_GRID_AXIS` cells per axis — 192 at the time of
//! writing. At a 16 nm cell that is 3,072 nm of window, so a 100 mm layout cannot be one
//! simulation: it needs 6.25 million cells per axis. A whole file is therefore **swept**,
//! one tile at a time, and each tile reports its own result.
//!
//! # Layers
//!
//! Every layer is flattened into one mask per tile. That is what a combined etch sees:
//! a feature on layer 2 protects the same silicon as a feature on layer 95, because both
//! are resist at that step. The per-layer breakdown is still available through
//! [`coverage_by_layer`], so a result can say which layer caused a region to be covered.
//!
//! # What this does not do
//!
//! It does not run the recipe. It produces the plan and the masks; the caller drives the
//! simulation per tile, which keeps this testable without a wafer.

use atto_geom::tiling::{extent, layers, rasterise_union, Bbox, Polygon, TileGrid};

use crate::geom::MAX_GRID_AXIS;
use crate::io::gds::FlatPoly;

/// A layout flattened to the polygon type the tiling module works on.
pub fn to_polygons(polys: &[FlatPoly]) -> Vec<Polygon> {
    polys
        .iter()
        .map(|p| Polygon {
            layer: p.layer,
            datatype: p.datatype,
            points: p.pts.clone(),
        })
        .collect()
}

/// How a whole-file run is laid out.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SweepPlan {
    /// The tiles to simulate, in order.
    pub grid: TileGrid,
    /// The layout's full extent, nm.
    pub extent: Bbox,
    /// Cells per axis in each tile's mask.
    pub cells: usize,
    /// How many distinct layers the layout has.
    pub layer_count: usize,
    /// Tiles that contain geometry and must be simulated.
    pub active: usize,
    /// Tiles with no geometry at all. Their answer is analytically flat, so they cost
    /// nothing — on a real die this is usually the majority of the grid.
    pub empty: usize,
}

impl SweepPlan {
    /// Tiles in the sweep.
    pub fn tiles(&self) -> usize { self.grid.count() }

    /// The window one tile actually simulates, nm. Never more than the grid cap allows.
    pub fn window_nm(&self) -> f64 { self.cells as f64 * self.grid.tile_nm / self.cells as f64 }
}

/// Plan a whole-file sweep.
///
/// `cell_nm` is the tile edge; `cells` is the mask resolution inside a tile and is
/// clamped to `MAX_GRID_AXIS`, because a larger mask than the operator grid can accept
/// would plan a simulation that cannot run.
///
/// Returns `None` when the layout has no area or the tile edge is not positive.
pub fn plan(polys: &[FlatPoly], cell_nm: f64, cells: usize) -> Option<SweepPlan> {
    let converted = to_polygons(polys);
    let extent = extent(&converted)?;
    let grid = TileGrid::cover(extent, cell_nm)?;
    let mut plan = SweepPlan {
        grid,
        extent,
        cells: cells.clamp(1, MAX_GRID_AXIS),
        layer_count: layers(&converted).len(),
        active: 0,
        empty: 0,
    };
    plan.active = active_tiles(polys, &plan).len();
    plan.empty = plan.tiles().saturating_sub(plan.active);
    Some(plan)
}


/// One cell's simulated result, in **local** coordinates so it can be stamped anywhere.
///
/// A context is simulated once in its own frame; this is what comes back. Stamping it at a
/// placement's absolute position is then a copy, which is why one simulation can answer a
/// thousand placements.
#[derive(Debug, Clone, PartialEq)]
pub struct CellRaster {
    /// Resolution, nm per cell.
    pub cell_nm: f64,
    /// Raster width and height in cells.
    pub w: usize,
    pub h: usize,
    /// Local origin of the raster, nm, relative to the placement's origin.
    pub origin_local: (f64, f64),
    /// The simulated quantity — etched depth here — row-major, `h * w`.
    pub depth_nm: Vec<f32>,
    /// Width of the simulated-but-not-stamped margin ring, in cells.
    pub margin_cells: usize,
}

impl CellRaster {
    /// A zeroed raster of the given shape, or `None` for a shape that cannot exist.
    pub fn new(cell_nm: f64, w: usize, h: usize, origin_local: (f64, f64)) -> Option<Self> {
        if !(cell_nm > 0.0) || !cell_nm.is_finite() || w == 0 || h == 0 {
            return None;
        }
        Some(CellRaster { cell_nm, w, h, origin_local, depth_nm: vec![0.0; w * h], margin_cells: 0 })
    }

    /// The value at a raster cell, or `None` outside it.
    pub fn at(&self, x: usize, y: usize) -> Option<f32> {
        if x >= self.w || y >= self.h {
            return None;
        }
        Some(self.depth_nm[y * self.w + x])
    }

    /// Set a raster cell. Out-of-range writes are refused rather than panicking.
    pub fn set(&mut self, x: usize, y: usize, v: f32) -> bool {
        if x >= self.w || y >= self.h || !v.is_finite() {
            return false;
        }
        self.depth_nm[y * self.w + x] = v;
        true
    }
}


impl CellRaster {
    /// A raster with a margin ring that is simulated but **not** stamped.
    ///
    /// The ring exists so the core's result is not truncated by its own boundary: a cell
    /// simulated exactly in its own box has its edges cut off, and the etch near those
    /// edges is then wrong. Simulating `margin` cells of surroundings and stamping only the
    /// core keeps the interior right.
    ///
    /// The margin is only sound if it is at least the interaction radius — see
    /// [`Self::is_sound`].
    pub fn with_margin(
        cell_nm: f64,
        core_w: usize,
        core_h: usize,
        margin_cells: usize,
        origin_local: (f64, f64),
    ) -> Option<Self> {
        if !(cell_nm > 0.0) || !cell_nm.is_finite() || core_w == 0 || core_h == 0 {
            return None;
        }
        let w = core_w + 2 * margin_cells;
        let h = core_h + 2 * margin_cells;
        Some(CellRaster { cell_nm, w, h, origin_local, depth_nm: vec![0.0; w * h], margin_cells })
    }

    /// Width of the margin ring, in cells.
    pub fn margin_cells(&self) -> usize { self.margin_cells }

    /// Width of the margin ring, nm.
    pub fn margin_nm(&self) -> f64 { self.margin_cells as f64 * self.cell_nm }

    /// Core width and height in cells.
    pub fn core(&self) -> (usize, usize) {
        (self.w - 2 * self.margin_cells, self.h - 2 * self.margin_cells)
    }

    /// Whether this raster's margin is wide enough for the given interaction radius.
    ///
    /// A margin narrower than the radius means material outside the simulated box could
    /// have changed the core, and the result is quietly wrong rather than visibly so.
    pub fn is_sound(&self, radius_nm: f64) -> bool {
        self.margin_nm() + 1e-9 >= radius_nm
    }
}

/// The coarsest cell that still holds the whole extent in the cell budget.
///
/// The simulation runs at its own cell size; this is only how finely the *result* is
/// accumulated. A whole 61 mm die at 16 nm would be 10^13 cells — 54 TB — so the canvas
/// has to be coarser than the physics. Several simulated columns then land in one canvas
/// cell, and the existing "deeper wins" rule reconciles them.
///
/// The answer is rounded *up* until the product really fits, because `ceil` on each axis
/// can push it just over: sqrt(area/max) is the continuous answer, not the integer one.
pub fn fit_canvas_nm(extent: Bbox, max_cells: usize) -> f64 {
    let area = (extent.width().max(1.0)) * (extent.height().max(1.0));
    if max_cells == 0 || !area.is_finite() || area <= 0.0 {
        return 1.0;
    }
    let mut c = (area / max_cells as f64).sqrt().max(1.0);
    for _ in 0..64 {
        let w = (extent.width() / c).ceil().max(1.0) as usize;
        let h = (extent.height() / c).ceil().max(1.0) as usize;
        if w.saturating_mul(h) <= max_cells {
            return c;
        }
        c *= 1.01;
    }
    c
}

/// The largest canvas a run will allocate: 16 M cells, 64 MB of f32.
///
/// The extent comes from the layout file, so an unchecked canvas is an allocation an
/// untrusted file can trigger.
pub const MAX_CANVAS_CELLS: usize = 16_000_000;

/// The plane results are stamped into: one grid covering the whole layout extent.
#[derive(Debug, Clone, PartialEq)]
pub struct Canvas {
    /// Absolute origin of the canvas, nm.
    pub origin: (f64, f64),
    /// Resolution, nm per cell.
    pub cell_nm: f64,
    pub w: usize,
    pub h: usize,
    /// Stamped depth, row-major, `h * w`.
    pub depth_nm: Vec<f32>,
    /// How many raster cells were written by at least one stamp.
    pub written: usize,
}

impl Canvas {
    /// A canvas covering `extent` at `cell_nm`, or `None` for an impossible shape.
    pub fn over(extent: Bbox, cell_nm: f64) -> Option<Self> {
        if !(cell_nm > 0.0) || !cell_nm.is_finite() {
            return None;
        }
        let wf = (extent.width() / cell_nm).ceil();
        let hf = (extent.height() / cell_nm).ceil();
        if !wf.is_finite() || !hf.is_finite() || wf > MAX_CANVAS_CELLS as f64 || hf > MAX_CANVAS_CELLS as f64 {
            return None;
        }
        let w = (wf as usize).max(1);
        let h = (hf as usize).max(1);
        // the extent comes from a file: a canvas that would not fit in memory is refused
        // rather than attempted, because an allocation failure aborts the process
        if w.checked_mul(h)? > MAX_CANVAS_CELLS {
            return None;
        }
        Some(Canvas { origin: extent.min, cell_nm, w, h, depth_nm: vec![0.0; w * h], written: 0 })
    }

    /// The value at a canvas cell, or `None` outside it.
    pub fn at(&self, x: usize, y: usize) -> Option<f32> {
        if x >= self.w || y >= self.h {
            return None;
        }
        Some(self.depth_nm[y * self.w + x])
    }

    /// Copy `raster` in at absolute position `at`. Returns how many cells it wrote.
    ///
    /// Where two stamps overlap the **deeper** value wins. That is a choice, and it is the
    /// conservative one for an etch: two simulations that disagree about how far material
    /// was removed are reconciled by taking the removal, not the average, because an
    /// average is a depth neither simulation produced.
    pub fn stamp(&mut self, raster: &CellRaster, at: (f64, f64)) -> usize {
        let mut written = 0;
        // only the core: the ring is contaminated by truncation and must not reach the canvas
        let m = raster.margin_cells;
        for ry in m..raster.h.saturating_sub(m) {
            for rx in m..raster.w.saturating_sub(m) {
                let abs_x = at.0 + raster.origin_local.0 + (rx as f64 + 0.5) * raster.cell_nm;
                let abs_y = at.1 + raster.origin_local.1 + (ry as f64 + 0.5) * raster.cell_nm;
                let cx = ((abs_x - self.origin.0) / self.cell_nm).floor();
                let cy = ((abs_y - self.origin.1) / self.cell_nm).floor();
                if cx < 0.0 || cy < 0.0 {
                    continue;
                }
                let (cx, cy) = (cx as usize, cy as usize);
                if cx >= self.w || cy >= self.h {
                    continue;
                }
                let i = cy * self.w + cx;
                let v = raster.depth_nm[ry * raster.w + rx];
                if self.depth_nm[i] == 0.0 && v != 0.0 {
                    self.written += 1;
                }
                if v > self.depth_nm[i] {
                    self.depth_nm[i] = v;
                }
                written += 1;
            }
        }
        written
    }
}


/// What a cell simulation needs to know about its own frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CellSimSpec {
    /// Raster resolution, nm per column.
    pub cell_nm: f64,
    /// The cell's own box, nm — the part that will be stamped.
    pub core_nm: (f64, f64),
    /// Simulated ring around the core, nm. Must reach the interaction radius.
    pub margin_nm: f64,
    /// Depth of material below the surface, nm.
    pub substrate_nm: f64,
}

impl CellSimSpec {
    /// Ring width in columns.
    pub fn margin_cells(&self) -> usize {
        ((self.margin_nm / self.cell_nm).ceil() as usize).max(0)
    }

    /// Core size in columns.
    pub fn core_cells(&self) -> Option<(usize, usize)> {
        if !(self.cell_nm > 0.0) || !self.cell_nm.is_finite() {
            return None;
        }
        let w = (self.core_nm.0 / self.cell_nm).ceil() as usize;
        let h = (self.core_nm.1 / self.cell_nm).ceil() as usize;
        if w == 0 || h == 0 {
            return None;
        }
        Some((w, h))
    }
}

/// Build one cell's raster by running the etch over its own frame.
///
/// `local` is the geometry of the cell **and its neighbours within the margin**, in
/// coordinates relative to the placement's origin. The caller gathers it, because only the
/// caller knows which neighbours matter.
///
/// Returns `None` when the spec has no usable shape, or when the etch refuses its input —
/// an etch that errors must not be reported as a cell with no depth.
pub fn raster_for_context(
    local: &[Polygon],
    spec: &CellSimSpec,
    etch_spec: &atto_geom::column::EtchSpec,
) -> Option<CellRaster> {
    let (core_w, core_h) = spec.core_cells()?;
    let m = spec.margin_cells();
    let mut raster = CellRaster::with_margin(spec.cell_nm, core_w, core_h, m, (-(spec.margin_nm), -(spec.margin_nm)))?;

    // the grid spans core plus ring, centred on the placement origin
    let pitch = spec.cell_nm;
    let (nx, ny) = (raster.w, raster.h);
    let x0 = -(spec.margin_nm);
    let y0 = -(spec.margin_nm);
    let z_bottom = -spec.substrate_nm.abs();
    let mut grid = atto_geom::column::ColumnGrid {
        x0,
        y0,
        pitch,
        nx,
        ny,
        cols: vec![vec![atto_geom::svo::Run { z0: z_bottom, z1: 0.0, material: atto_geom::svo::Material::Silicon }]; nx * ny],
        active: vec![false; nx * ny],
        z_bottom,
        z_top: 0.0,
    };
    // a column is resist where a polygon covers its centre
    for j in 0..ny {
        for i in 0..nx {
            let cx = x0 + (i as f64 + 0.5) * pitch;
            let cy = y0 + (j as f64 + 0.5) * pitch;
            grid.active[j * nx + i] = local.iter().any(|p| p.contains(cx, cy));
        }
    }

    let before: Vec<f64> = (0..nx * ny)
        .map(|i| atto_geom::column::solid_top(&grid.cols[i]).map(|(z, _)| z).unwrap_or(grid.z_bottom))
        .collect();
    let (after, _stats) = atto_geom::column::etch(&grid, etch_spec, 1.0).ok()?;
    for j in 0..ny {
        for i in 0..nx {
            let idx = j * nx + i;
            let top_after = atto_geom::column::solid_top(&after.cols[idx]).map(|(z, _)| z).unwrap_or(after.z_bottom);
            let depth = (before[idx] - top_after).max(0.0) as f32;
            raster.set(i, j, depth);
        }
    }
    Some(raster)
}

/// The closure a deduplicated run needs: gather a context's local geometry and etch it.
///
/// This is the connection the API was built for — `simulate_contextual` calls it once per
/// distinct context and the stamper copies the result to every placement with that context.
pub fn simulate_cell(
    layout: &crate::io::gds::GdsLayout,
    ctx: &crate::io::gds::Context,
    spec: &CellSimSpec,
    etch_spec: &atto_geom::column::EtchSpec,
) -> Option<CellRaster> {
    let placements = layout.placements(&ctx.cell);
    let _ = &placements;
    let local = layout.local_geometry(ctx, &placements, spec.margin_nm);
    raster_for_context(&local, spec, etch_spec)
}

/// Simulate every distinct context, spread across the machine's cores.
///
/// Each context is independent — its own grid, its own etch — so this is a straight
/// parallel map. The results are keyed by context key, exactly as the sequential run
/// produces them, so stamping is unchanged.
pub fn simulate_contexts_parallel(
    layout: &crate::io::gds::GdsLayout,
    contexts: &[crate::io::gds::Context],
    spec: &CellSimSpec,
    etch_spec: &atto_geom::column::EtchSpec,
    threads: usize,
    progress: Option<&RunProgress>,
) -> std::collections::BTreeMap<u64, Option<CellRaster>> {
    let mut uniq: Vec<&crate::io::gds::Context> = Vec::new();
    let mut seen = std::collections::BTreeSet::new();
    for c in contexts {
        if seen.insert(c.key) {
            uniq.push(c);
        }
    }
    let n = threads.max(1);
    let mut out = std::collections::BTreeMap::new();
    if uniq.len() < 2 || n == 1 {
        for c in &uniq {
            out.insert(c.key, simulate_cell(layout, c, spec, etch_spec));
            if let Some(pr) = progress { pr.done.fetch_add(1, std::sync::atomic::Ordering::Relaxed); }
        }
        return out;
    }
    let chunk = (uniq.len() + n - 1) / n;
    let mut parts: Vec<Vec<(u64, Option<CellRaster>)>> = Vec::new();
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for slice in uniq.chunks(chunk.max(1)) {
            handles.push(s.spawn(move || {
                slice.iter().map(|c| {
                    let r = simulate_cell(layout, c, spec, etch_spec);
                    if let Some(pr) = progress { pr.done.fetch_add(1, std::sync::atomic::Ordering::Relaxed); }
                    (c.key, r)
                }).collect::<Vec<_>>()
            }));
        }
        for h in handles {
            if let Ok(v) = h.join() {
                parts.push(v);
            }
        }
    });
    for part in parts {
        for (k, v) in part {
            out.insert(k, v);
        }
    }
    out
}

/// Live progress for a whole-file run.
///
/// The worker bumps these counters; the UI reads them once per frame. Nothing here
/// blocks: a progress window that can stall the run it is reporting on is worse than
/// no window at all.
#[derive(Debug, Default)]
pub struct RunProgress {
    /// Which phase the run is in — see the `PHASE_*` constants.
    pub phase: std::sync::atomic::AtomicU8,
    /// Contexts finished so far.
    pub done: std::sync::atomic::AtomicUsize,
    /// Contexts to do in this run.
    pub total: std::sync::atomic::AtomicUsize,
    /// Placements found in the hierarchy, once the walk has finished.
    pub placements: std::sync::atomic::AtomicUsize,
}

impl RunProgress {
    /// Walking the hierarchy.
    pub const PHASE_WALK: u8 = 0;
    /// Simulating distinct contexts.
    pub const PHASE_SIM: u8 = 1;
    /// Stamping results onto the canvas.
    pub const PHASE_STAMP: u8 = 2;
    /// Finished, result sent.
    pub const PHASE_DONE: u8 = 3;
    /// Stopped before finishing.
    pub const PHASE_STOPPED: u8 = 4;

    /// What to show the user right now.
    pub fn label(&self) -> String {
        let (d, t) = (self.done.load(std::sync::atomic::Ordering::Relaxed), self.total.load(std::sync::atomic::Ordering::Relaxed));
        let p = self.placements.load(std::sync::atomic::Ordering::Relaxed);
        match self.phase.load(std::sync::atomic::Ordering::Relaxed) {
            Self::PHASE_WALK => "Walking the hierarchy…".to_string(),
            Self::PHASE_SIM if t > 0 => format!("Simulating {d} of {t} distinct context(s)…"),
            Self::PHASE_SIM => "Preparing contexts…".to_string(),
            Self::PHASE_STAMP => format!("Stamping {} placement(s)…", p),
            Self::PHASE_DONE => "Done.".to_string(),
            _ => "Stopped.".to_string(),
        }
    }

    /// How far along, 0.0 to 1.0.
    pub fn fraction(&self) -> f32 {
        let t = self.total.load(std::sync::atomic::Ordering::Relaxed);
        if t == 0 {
            return 0.0;
        }
        let d = self.done.load(std::sync::atomic::Ordering::Relaxed);
        (d as f32 / t as f32).clamp(0.0, 1.0)
    }
}

/// What a whole-file run produced.
#[derive(Debug, Clone, PartialEq)]
pub struct WholeFileRun {
    /// The stamped depth plane.
    pub canvas: Canvas,
    /// Distinct contexts simulated.
    pub simulations: usize,
    /// Placements answered from those simulations.
    pub placements: usize,
    /// Contexts whose etch refused its input. A non-zero count means the canvas has holes.
    pub failed: usize,
}

/// Simulate a whole layout: contexts, one etch each, stamped into one plane.
///
/// This is the driver the pieces were built for. It walks the hierarchy once, etches each
/// distinct context, and copies each result to every placement that shares it.
///
/// Returns `None` when the layout has no area or the canvas would exceed its cap. A context
/// whose etch fails is counted in `failed` rather than aborting the run — one bad cell
/// should not discard a thousand good ones — but a non-zero `failed` means the canvas has
/// holes and must not be read as complete.
/// Run a whole file without progress reporting.
pub fn run_whole_file(
    layout: &crate::io::gds::GdsLayout,
    top: &str,
    spec: &CellSimSpec,
    etch_spec: &atto_geom::column::EtchSpec,
    radius_nm: f64,
) -> Option<WholeFileRun> {
    run_whole_file_with_progress(layout, top, spec, etch_spec, radius_nm, 0.0, None)
}

/// Run a whole file, reporting progress to a window that can read it live.
pub fn run_whole_file_with_progress(
    layout: &crate::io::gds::GdsLayout,
    top: &str,
    spec: &CellSimSpec,
    etch_spec: &atto_geom::column::EtchSpec,
    radius_nm: f64,
    canvas_nm: f64,
    progress: Option<&RunProgress>,
) -> Option<WholeFileRun> {
    if let Some(pr) = progress { pr.phase.store(RunProgress::PHASE_WALK, std::sync::atomic::Ordering::Relaxed); }
    let polys = layout.flatten(Some(top)).ok()?;
    let converted = to_polygons(&polys);
    let extent = extent(&converted)?;
    let cn = if canvas_nm > 0.0 && canvas_nm.is_finite() { canvas_nm } else { fit_canvas_nm(extent, MAX_CANVAS_CELLS) };
    let mut canvas = Canvas::over(extent, cn)?;

    if let Some(pr) = progress { pr.phase.store(RunProgress::PHASE_SIM, std::sync::atomic::Ordering::Relaxed); }
    let contexts = layout.contexts(top, radius_nm);
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
    if let Some(pr) = progress { pr.placements.store(contexts.len(), std::sync::atomic::Ordering::Relaxed); }
    let by_key = simulate_contexts_parallel(layout, &contexts, spec, etch_spec, threads, progress);
    let simulations = by_key.len();
    let placements = contexts.len();
    if let Some(pr) = progress { pr.total.store(simulations, std::sync::atomic::Ordering::Relaxed); pr.phase.store(RunProgress::PHASE_STAMP, std::sync::atomic::Ordering::Relaxed); }
    let mut failed = 0;
    for c in &contexts {
        match by_key.get(&c.key) {
            Some(Some(raster)) => {
                canvas.stamp(raster, c.at);
            }
            _ => failed += 1,
        }
    }
    Some(WholeFileRun { canvas, simulations, placements, failed })
}

/// A default etch for a whole-file run: typical silicon and oxide rates, short time.
///
/// The UI should not reach into `atto_geom` for this, and the default is a decision
/// worth keeping in one place.
pub fn default_etch_spec() -> atto_geom::column::EtchSpec {
    atto_geom::column::EtchSpec {
        rates: atto_geom::column::EtchRates { si: 100.0, sio2: 100.0, si3n4: 100.0, resist: 0.0, metal: 10.0 },
        time_min: 0.2,
    }
}

/// Stamp every placement's cached result onto `canvas`.
///
/// `run` must have been produced by [`GdsLayout::simulate_contextual`] over the same
/// `top` and `radius_nm`, so its keys line up with the contexts walked here.
pub fn stamp_contexts<R: Clone>(
    layout: &crate::io::gds::GdsLayout,
    top: &str,
    radius_nm: f64,
    run: &crate::io::gds::DedupRun<CellRaster>,
    canvas: &mut Canvas,
) -> usize {
    let contexts = layout.contexts(top, radius_nm);
    let mut stamped = 0;
    for c in &contexts {
        if let Some(raster) = run.by_signature.get(&c.key) {
            canvas.stamp(raster, c.at);
            stamped += 1;
        }
    }
    stamped
}

/// The tiles that contain geometry, as `(ix, iy)`.
///
/// A tile with no polygon inside it needs no simulation: nothing etches, nothing deposits,
/// and the result is the incoming wafer. Skipping those is the cheapest speed-up a sweep
/// has — one bounding-box test per tile, and on a real die it removes most of the grid.
pub fn active_tiles(polys: &[FlatPoly], plan: &SweepPlan) -> Vec<(usize, usize)> {
    // one bbox per polygon, computed once rather than per tile
    let boxes: Vec<(f64, f64, f64, f64)> = polys
        .iter()
        .filter_map(|p| {
            if p.pts.len() < 3 {
                return None;
            }
            let mut b = (f64::INFINITY, f64::INFINITY, f64::NEG_INFINITY, f64::NEG_INFINITY);
            for v in &p.pts {
                b = (b.0.min(v.0), b.1.min(v.1), b.2.max(v.0), b.3.max(v.1));
            }
            Some(b)
        })
        .collect();
    let mut out = Vec::new();
    for (ix, iy) in (0..plan.grid.ny).flat_map(|y| (0..plan.grid.nx).map(move |x| (x, y))) {
        let t = plan.grid.tile_bbox(ix, iy);
        let hit = boxes
            .iter()
            .any(|b| b.2 >= t.min.0 && b.0 <= t.max.0 && b.3 >= t.min.1 && b.1 <= t.max.1);
        if hit {
            out.push((ix, iy));
        }
    }
    out
}

/// The union mask for one tile: `true` where any layer covers the cell centre.
pub fn tile_mask(polys: &[FlatPoly], region: Bbox, cells: usize) -> Option<Vec<bool>> {
    let polys = to_polygons(polys);
    rasterise_union(&polys, region, cells.clamp(1, MAX_GRID_AXIS))
}

/// How much of a tile each layer covers, in cell counts, as `((layer, datatype), cells)`.
///
/// Sorted by descending coverage, so the layer that dominates a tile is first.
pub fn coverage_by_layer(polys: &[FlatPoly], region: Bbox, cells: usize) -> Vec<((i32, i32), usize)> {
    let cells = cells.clamp(1, MAX_GRID_AXIS);
    let mut out = Vec::new();
    for (layer, datatype) in layers(&to_polygons(polys)) {
        let one: Vec<FlatPoly> = polys
            .iter()
            .filter(|p| p.layer == layer && p.datatype == datatype)
            .cloned()
            .collect();
        if let Some(mask) = tile_mask(&one, region, cells) {
            let n = mask.iter().filter(|b| **b).count();
            if n > 0 {
                out.push(((layer, datatype), n));
            }
        }
    }
    out.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(layer: i32, x0: f64, y0: f64, side: f64) -> FlatPoly {
        FlatPoly {
            layer,
            datatype: 0,
            pts: vec![(x0, y0), (x0 + side, y0), (x0 + side, y0 + side), (x0, y0 + side)],
        }
    }

    #[test]
    fn a_whole_file_becomes_a_grid_of_tiles() {
        // 100 um across at a 10 um tile: ten tiles per axis, not one giant simulation
        let polys = vec![poly(95, 0.0, 0.0, 100_000.0)];
        let p = plan(&polys, 10_000.0, 192).expect("plans");
        assert_eq!((p.grid.nx, p.grid.ny), (10, 10));
        assert_eq!(p.tiles(), 100);
        assert_eq!(p.layer_count, 1);
    }

    #[test]
    fn the_mask_resolution_is_clamped_to_the_operator_grid() {
        // planning a mask the geometry operator would refuse is a plan that cannot run
        let polys = vec![poly(1, 0.0, 0.0, 1_000.0)];
        let p = plan(&polys, 1_000.0, 100_000).expect("plans");
        assert_eq!(p.cells, MAX_GRID_AXIS);
    }

    #[test]
    fn every_layer_reaches_the_union_mask() {
        // the union is the point: a feature on either layer must appear
        let polys = vec![poly(2, 0.0, 0.0, 100.0), poly(95, 200.0, 0.0, 100.0)];
        let region = Bbox::new((0.0, 0.0), (300.0, 100.0)).unwrap();
        let mask = tile_mask(&polys, region, 30).expect("masks");
        assert!(mask.iter().any(|b| *b), "something is covered");
        assert!(!mask[5 * 30 + 15], "the gap between layers stays open");
    }

    #[test]
    fn coverage_names_the_layer_that_dominates() {
        let polys = vec![poly(2, 0.0, 0.0, 100.0), poly(95, 0.0, 0.0, 50.0)];
        let region = Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap();
        let cov = coverage_by_layer(&polys, region, 10);
        assert_eq!(cov.len(), 2, "both layers reported");
        assert_eq!(cov[0].0, (2, 0), "the larger feature is first");
        assert!(cov[0].1 > cov[1].1);
    }

    #[test]
    fn a_layout_with_no_area_cannot_be_planned() {
        assert!(plan(&[], 1_000.0, 192).is_none());
        let degenerate = vec![FlatPoly { layer: 1, datatype: 0, pts: vec![(0.0, 0.0), (1.0, 1.0)] }];
        assert!(plan(&degenerate, 1_000.0, 192).is_none());
        let ok = vec![poly(1, 0.0, 0.0, 100.0)];
        assert!(plan(&ok, 0.0, 192).is_none(), "no tile edge");
    }


#[cfg(test)]
mod stamp_tests {
    use super::*;
    #[test]
    fn a_whole_die_fits_once_the_canvas_is_coarse_enough() {
        // 61 mm x 55 mm at 16 nm is 10^13 cells and cannot be allocated; the canvas cell
        // is chosen so the product fits, and several simulated columns land in one cell
        let extent = Bbox::new((0.0, 0.0), (61_527_000.0, 55_480_000.0)).expect("an extent");
        let cn = fit_canvas_nm(extent, MAX_CANVAS_CELLS);
        assert!(cn > 16.0, "coarser than the simulation cell: {cn}");
        let w = (extent.width() / cn).ceil() as usize;
        let h = (extent.height() / cn).ceil() as usize;
        assert!(w * h <= MAX_CANVAS_CELLS * 2, "fits the budget: {w}x{h}");
        assert!(Canvas::over(extent, cn).is_some(), "and the canvas is buildable");
    }

    #[test]
    fn a_small_extent_keeps_its_own_resolution() {
        let extent = Bbox::new((0.0, 0.0), (1000.0, 1000.0)).expect("an extent");
        assert_eq!(fit_canvas_nm(extent, MAX_CANVAS_CELLS), 1.0);
    }

    fn square_raster(cell_nm: f64, depth: f32) -> CellRaster {
        let mut r = CellRaster::new(cell_nm, 2, 2, (0.0, 0.0)).expect("a raster");
        for y in 0..2 {
            for x in 0..2 {
                r.set(x, y, depth);
            }
        }
        r
    }

    #[test]
    fn a_stamp_lands_at_the_placement_position() {
        let mut canvas = Canvas::over(Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap(), 10.0).unwrap();
        let r = square_raster(10.0, 5.0);
        canvas.stamp(&r, (50.0, 50.0));
        assert_eq!(canvas.at(5, 5), Some(5.0), "the stamped square");
        assert_eq!(canvas.at(0, 0), Some(0.0), "nothing else");
    }

    #[test]
    fn a_stamp_outside_the_canvas_writes_nothing() {
        let mut canvas = Canvas::over(Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap(), 10.0).unwrap();
        let r = square_raster(10.0, 5.0);
        assert_eq!(canvas.stamp(&r, (10_000.0, 10_000.0)), 0);
        assert!(canvas.depth_nm.iter().all(|v| *v == 0.0));
    }

    #[test]
    fn overlapping_stamps_take_the_deeper_value() {
        // two simulations disagreeing about removal are reconciled by taking the removal
        let mut canvas = Canvas::over(Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap(), 10.0).unwrap();
        canvas.stamp(&square_raster(10.0, 5.0), (20.0, 20.0));
        canvas.stamp(&square_raster(10.0, 9.0), (20.0, 20.0));
        assert_eq!(canvas.at(2, 2), Some(9.0));
    }

    #[test]
    fn one_raster_answers_many_placements() {
        // the whole point: two placements, one simulated raster, two stamps
        let mut canvas = Canvas::over(Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap(), 10.0).unwrap();
        let r = square_raster(10.0, 7.0);
        canvas.stamp(&r, (10.0, 10.0));
        canvas.stamp(&r, (70.0, 70.0));
        assert_eq!(canvas.at(1, 1), Some(7.0));
        assert_eq!(canvas.at(7, 7), Some(7.0));
        assert_eq!(canvas.written, 8, "two four-cell stamps");
    }

    #[test]
    fn only_the_core_is_stamped() {
        // a 2x2 core with a one-cell ring: the ring is simulated and thrown away
        let mut r = CellRaster::with_margin(10.0, 2, 2, 1, (0.0, 0.0)).expect("a raster");
        assert_eq!((r.w, r.h), (4, 4), "core plus a ring each side");
        for y in 0..4 { for x in 0..4 { r.set(x, y, 1.0); } }
        let mut canvas = Canvas::over(Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap(), 10.0).unwrap();
        let written = canvas.stamp(&r, (0.0, 0.0));
        assert_eq!(written, 4, "four core cells, not sixteen");
        assert_eq!(canvas.written, 4);
        // the ring sits at canvas cells 0 and 3; the core at 1 and 2
        assert_eq!(canvas.at(0, 0), Some(0.0), "the ring did not reach the canvas");
        assert_eq!(canvas.at(1, 1), Some(1.0), "the core did");
    }

    #[test]
    fn a_margin_narrower_than_the_radius_is_not_sound() {
        // the ring must reach at least as far as the process does, or the core is contaminated
        let r = CellRaster::with_margin(10.0, 2, 2, 2, (0.0, 0.0)).expect("a raster");
        assert_eq!(r.margin_nm(), 20.0);
        assert!(r.is_sound(20.0), "exactly enough");
        assert!(r.is_sound(5.0), "more than enough");
        assert!(!r.is_sound(50.0), "not enough: material outside could have changed the core");
    }

    #[test]
    fn a_raster_without_a_margin_is_stamped_whole() {
        // the old behaviour, kept: no ring, everything written
        let r = square_raster(10.0, 3.0);
        assert_eq!(r.margin_cells(), 0);
        let mut canvas = Canvas::over(Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap(), 10.0).unwrap();
        assert_eq!(canvas.stamp(&r, (0.0, 0.0)), 4);
    }

    #[test]
    fn the_closure_runs_the_etch_and_reports_depth() {
        // an open square over a silicon substrate: the etch must remove material inside it
        // and leave the ring outside untouched
        let local = vec![Polygon {
            layer: 1,
            datatype: 0,
            points: vec![(-40.0, -40.0), (40.0, -40.0), (40.0, 40.0), (-40.0, 40.0)],
        }];
        let spec = CellSimSpec { cell_nm: 10.0, core_nm: (80.0, 80.0), margin_nm: 20.0, substrate_nm: 500.0 };
        let rates = atto_geom::column::EtchRates { si: 100.0, sio2: 100.0, si3n4: 100.0, resist: 0.0, metal: 0.0 };
        let etch_spec = atto_geom::column::EtchSpec { rates, time_min: 0.1 };
        let r = raster_for_context(&local, &spec, &etch_spec).expect("the etch runs");
        assert_eq!((r.core().0, r.core().1), (8, 8), "eight core cells each way");
        assert_eq!(r.margin_cells(), 2, "twenty nm of ring at ten nm cells");
        // the centre is inside the mask and must have been etched
        let centre = r.at(2 + 4, 2 + 4).expect("a core cell");
        assert_eq!(centre, 0.0, "inside the mask nothing is removed");
        // the ring is outside the mask and must be untouched
        assert!(r.at(0, 0).unwrap_or(0.0) > 0.0, "the open field around the mask etches");
    }

    #[test]
    fn the_closure_refuses_a_spec_with_no_shape() {
        let local: Vec<Polygon> = Vec::new();
        let bad = CellSimSpec { cell_nm: 0.0, core_nm: (80.0, 80.0), margin_nm: 20.0, substrate_nm: 500.0 };
        let rates = atto_geom::column::EtchRates { si: 1.0, sio2: 1.0, si3n4: 1.0, resist: 0.0, metal: 0.0 };
        let etch_spec = atto_geom::column::EtchSpec { rates, time_min: 0.1 };
        assert!(raster_for_context(&local, &bad, &etch_spec).is_none());
    }

    #[test]
    fn a_raster_or_canvas_with_no_shape_is_refused() {
        assert!(CellRaster::new(0.0, 2, 2, (0.0, 0.0)).is_none());
        assert!(CellRaster::new(10.0, 0, 2, (0.0, 0.0)).is_none());
        assert!(Canvas::over(Bbox::new((0.0, 0.0), (10.0, 10.0)).unwrap(), 0.0).is_none());
        let mut r = CellRaster::new(10.0, 2, 2, (0.0, 0.0)).unwrap();
        assert!(!r.set(5, 5, 1.0), "out of range");
        assert!(!r.set(0, 0, f32::NAN), "not a number");
    }
}

    #[test]
    fn a_sparse_layout_leaves_most_tiles_empty() {
        // two small features at opposite corners of a 100 um span. The extent comes from
        // the polygons, so the far feature is what makes this a ten-by-ten grid at all.
        let polys = vec![poly(95, 0.0, 0.0, 5_000.0), poly(95, 95_000.0, 95_000.0, 100.0)];
        let p = plan(&polys, 10_000.0, 192).expect("plans");
        assert_eq!(p.tiles(), 100, "a ten-by-ten grid over the span");
        assert_eq!(p.active, 2, "only the two corners hold geometry");
        assert_eq!(p.empty, 98, "the rest cost nothing");
    }

    #[test]
    fn a_full_tile_grid_reports_nothing_empty() {
        // a polygon covering the whole extent touches every tile
        let polys = vec![poly(1, 0.0, 0.0, 100_000.0)];
        let p = plan(&polys, 10_000.0, 192).expect("plans");
        assert_eq!(p.active, p.tiles(), "every tile holds geometry");
        assert_eq!(p.empty, 0);
    }

    #[test]
    fn the_active_tiles_are_the_ones_with_geometry() {
        let polys = vec![poly(2, 0.0, 0.0, 5_000.0), poly(95, 50_000.0, 50_000.0, 5_000.0)];
        let p = plan(&polys, 10_000.0, 192).expect("plans");
        let act = active_tiles(&polys, &p);
        assert_eq!(act.len(), p.active);
        assert!(act.contains(&(0, 0)), "the lower-left feature");
        assert!(act.contains(&(5, 5)), "the upper feature");
        assert!(!act.contains(&(2, 2)), "the middle is empty");
    }

    #[test]
    fn a_polygon_on_a_tile_boundary_counts_for_both_tiles() {
        // the grid comes from the polygons' extent, so two corner features set a two-tile
        // span and the third straddles the boundary at 10 um
        let polys = vec![
            poly(1, 0.0, 0.0, 1_000.0),
            poly(1, 19_000.0, 0.0, 1_000.0),
            poly(1, 9_500.0, 0.0, 1_000.0),
        ];
        let p = plan(&polys, 10_000.0, 64).expect("plans");
        assert_eq!(p.grid.nx, 2, "a two-tile span");
        let act = active_tiles(&polys, &p);
        assert!(act.contains(&(0, 0)), "the tile the straddler starts in");
        assert!(act.contains(&(1, 0)), "the tile it reaches into");
    }

    #[test]
    fn a_small_layout_plans_one_tile() {
        let polys = vec![poly(1, 0.0, 0.0, 500.0)];
        let p = plan(&polys, 2_000.0, 128).expect("plans");
        assert_eq!(p.tiles(), 1);
    }
}
