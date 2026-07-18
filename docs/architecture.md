# AttoFab architecture

## Guiding principle: start dense, prove correctness, then add sparsity

Every mesh/solver decision in this codebase follows one rule: get a plain,
dense, obviously-correct implementation working and regression-tested
*first*. Only after that baseline is validated do we introduce a memory or
performance optimization (RLE, chunked sparsity, GPU, etc.) on top of it,
re-running the same tests to prove the optimization didn't change the
physics.

This is why:
- `Wafer1d` and `Wafer2d` both use plain `Vec<f64>` / flattened arrays, not
  an octree or sparse structure, even though the original design brief
  suggested octrees as "highly recommended." True adaptive-resolution
  octrees break the uniform-grid-spacing assumption that explicit
  finite-difference diffusion depends on (`C[i+1] - 2C[i] + C[i-1]`
  assumes constant `dx`); reconciling that is a well-known source of bugs
  in real TCAD tools and not worth the complexity for an educational v1.
- The recommended path to wafer-scale (Phase 3) is: bulk silicon stored as
  run-length-encoded columns (diffusion has a finite, analytically-known
  penetration depth per anneal, so most of the bulk never needs dense
  storage) + a chunked sparse grid of dense, uniform-spacing chunks near
  the active surface. Octrees and sparse voxel DAGs are demoted to
  rendering/export-time optimizations, not simulation-time data
  structures.

## Phases

- **Phase 1 (done)**: 1D column (`wafer1d.rs`) - Deal-Grove oxidation,
  Gaussian implant, Fick's-Law FDM diffusion. Ported from, and
  regression-tested against, a working JS reference implementation.
- **Phase 2 (done)**: 2D cross-section grid (`grid2d.rs`) - masked
  (LOCOS-style) oxidation, masked implant, full 2D Fick's-Law diffusion,
  lithography mask expose/develop (`lithography.rs`), and CA-based
  anisotropic/isotropic etch + deposition (`etch.rs`). A canvas renderer
  (`web/index.html`) visualizes the JSON-exported grid.
- **Phase 3 (planned)**: whole-wafer 3D via GDSII-parsed masks, RLE bulk +
  chunked sparse active zone, WebGL topography viewer.

## Crate layout

- `crates/core` (`attofab-core`) - the physics engine, pure Rust, no FFI.
  This is the crate all correctness work happens in.
- `crates/core-wasm`, `crates/core-py` - planned FFI bindings (wasm-bindgen
  / PyO3) that wrap `attofab-core`. Not yet scaffolded; add back to the
  workspace `members` list in the root `Cargo.toml` once they have a
  manifest.
- `backend/` - FastAPI backend, subprocess-bridges to the compiled
  `recipe_runner` binary (see `crates/core/src/bin/recipe_runner.rs`).
- `frontend_react/` (recommended UI) and `frontend/` (Streamlit,
  alternative UI) both talk to `backend/` over HTTP. `web/index.html`
  remains as a dependency-free, backend-free reference renderer (loads a
  JSON export directly, no server required) and the visual spec
  `frontend_react/src/components/WaferCanvas.jsx` reproduces.
- `bot.py` - CLI entrypoint that calls the engine directly (via
  `backend/core/engine_bridge.py`), no server needed.

## Numerical stability rules (do not violate)

- Diffusion is explicit FDM. The CFL-like stability condition
  `D * dt / dx^2 <= 0.5` (1D) / `D * dt * (1/dx^2 + 1/dy^2) <= 0.5` (2D)
  is enforced in code via an auto-computed, sub-stepped `dt` - never accept
  a caller-supplied `dt` that bypasses this.
- Oxidation (Deal-Grove) is closed-form, so it's unconditionally stable -
  no iteration, no timestep to tune.
- Etch/deposit (cellular automaton) is stable by construction as long as
  it stays nearest-neighbor-per-step; do not "speed it up" with multi-cell
  jumps per iteration, since that changes the physical undercut profile.

## Regression testing

`wafer1d.rs`'s test suite is pinned to the exact output of the original
validated JS Phase 1 engine (`node demo.js` in the JS prototype), not to a
hand-derived textbook estimate - this caught a real bug once already
(a stale hardcoded oxide-thickness expectation that didn't match either
engine). When adding new physics, always check tests against a second
independent computation (closed-form formula, the JS engine, or a known
analytical solution) rather than asserting whatever the code currently
outputs.

## IP compliance and export control rationale

AttoFab's exclusive reliance on public-domain physics formulas (Deal-Grove,
Fick's Laws, Pearson-IV, Dill/Mack, Level-Set - all textbook/peer-reviewed
science, none of it proprietary foundry data) is a deliberate legal
strategy, not just a scientific one:

- Universal physical laws and published constants cannot be copyrighted or
  patented, so an original implementation of them is a new, non-derivative
  work - no clean-room reverse-engineering process is needed because there
  is no proprietary source being cloned.
- Permissive licensing (MIT here; Apache 2.0 is the SkyWater PDK/DEVSIM/IHP
  precedent) avoids copyleft contamination that would make the engine
  unusable inside commercial or proprietary workflows.
- Software that is genuinely public-domain-physics-based, openly published,
  and free of NDAs or proprietary process data generally qualifies for
  EAR99 classification under U.S. export control rules, unlike commercial
  TCAD/EDA tools (which fall under ECCN 3D001/3E001/3E991 and face
  deemed-export licensing requirements for foreign national collaborators).

Practical implication for contributors: never add code that encodes a
specific real foundry's proprietary process recipe, calibrated fitting
constants from a non-public source, or reverse-engineered behavior from a
commercial tool's output. Constants must trace to a citable public source
(see `docs/math_references.md` for the standard this repo holds itself
to) - this is what keeps the "public domain physics only" claim legally
meaningful rather than aspirational.

## Development-language layer boundaries (why each piece is what it is)

- **Rust core** (`crates/core`): all physics/numerics live here for memory
  safety without GC overhead, and so results are reproducible independent
  of any scripting layer. This is the only crate correctness work should
  target directly.
- **`recipe_runner` binary + Python `backend/`**: the backend shells out to
  a compiled binary rather than embedding via PyO3, trading a small amount
  of subprocess overhead for zero compiled-extension/toolchain-version
  coupling - any Python environment that can spawn a process can drive the
  real engine. If profiling ever shows this boundary is a bottleneck (it
  shouldn't be, for reasonable recipe/grid sizes), PyO3 bindings
  (`crates/core-py`, currently unscaffolded) are the documented upgrade
  path.
- **`web/index.html`**: a dependency-free canvas renderer that reads the
  engine's JSON export directly - useful standalone and as the reference
  implementation for what `frontend_react/` should visually reproduce.
