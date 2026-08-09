# AttoFab

**An open-source, education-first electronics fabrication simulator.**

AttoFab simulates the core physics of semiconductor manufacturing — thermal
oxidation, ion implantation, thermal diffusion, photolithography, and
etching — using only classical, public-domain textbook physics (Deal-Grove,
Fick's Laws, Gaussian implant profiles). No proprietary foundry data, no
NDAs, no black boxes: every constant in this engine traces back to a
citable public source (see [`docs/math_references.md`](docs/math_references.md)).

```
[ Photoresist ]         growing field oxide, masked active area,
[ SiO2 ][Metal][ SiO2 ]  metal contact trench — a real AttoFab render
[     Silicon      ]
```

## Why

Commercial TCAD/EDA tools cost anywhere from ~$5M (130nm) to $500M+ (5nm)
in licensing, are closed-source, and are entangled in export-control
restrictions (ECCN 3D001/3E001/3E991) that create real friction for
international students and researchers. AttoFab sidesteps all of that: it
implements universal physical laws (which can't be copyrighted or
patented), under a permissive MIT license, with no proprietary process
data — the kind of software that generally qualifies for EAR99
classification. See [`docs/architecture.md`](docs/architecture.md#ip-compliance-and-export-control-rationale)
for the full rationale.

## What's implemented

| Phase | Status | What it does |
|---|---|---|
| **Phase 1 — 1D engine** | ✅ Done | Single vertical column: Deal-Grove oxidation, Gaussian implant, Fick's-Law diffusion (explicit FDM), junction-depth finding |
| **Phase 2 — 2D engine** | ✅ Done | Full 2D wafer cross-section: masked (LOCOS-style) oxidation/implant, 2D diffusion, mask expose/develop, cellular-automaton anisotropic & isotropic etch + deposition, canvas renderer |
| **Phase 3 — GDSII + 3D** | 🟡 Core done, integration partial | GDSII parser (hand-rolled, tested against known byte patterns), dense 3D wafer mesh (masked oxidation/implant, 3D diffusion, 3D anisotropic etch), WebGL heightfield topography viewer. **Not yet done:** 3D lithography expose/develop, a GDS-aware CLI recipe runner (3D flows currently only run via a Rust example), RLE/chunked-sparse memory, true Level-Set/HRLE topography evolution |

Also implemented: crystal-orientation-dependent oxidation ((111) vs (100)),
a full-stack backend + two frontends + CLI, and a physics regression bench
that runs on every backend startup.

**Known limitations, documented honestly rather than hidden:**
- The Pearson-IV implant distribution (for ion channeling tails) has correct, tested ODE-integration machinery, but its moment-matching formula isn't yet producing correct output — explicitly gated off (Gaussian remains the default), with the two failing tests marked `#[ignore]` and explained rather than silently shipped.
- `web/topography_3d.html` (the WebGL viewer) hasn't been visually verified in an actual browser in this development environment — the underlying computed topology has been verified independently (Rust tests + a from-scratch canvas render), but the Three.js file itself should be treated as unverified until someone opens it.
- Oxidation currently substitutes material in place rather than modeling the volume expansion that would push the surface upward — so oxide growth alone doesn't create 3D topography; only etch does.

## Project structure

```
open-fab-sim/
├── attofab_logo.svg / attofab_mark.svg   Logo (wordmark+mark / icon-only)
│
├── crates/
│   ├── core/                  attofab-core — the Rust physics engine
│   │   ├── src/
│   │   │   ├── materials.rs       Public-domain constants (Deal-Grove, diffusivities)
│   │   │   ├── oxidation.rs       Deal-Grove closed-form solver
│   │   │   ├── implant.rs         Gaussian (+ experimental Pearson-IV) implant profiles
│   │   │   ├── diffusion.rs       Fick's 2nd Law, explicit FDM, 1D
│   │   │   ├── wafer1d.rs         1D wafer column state (Phase 1)
│   │   │   ├── grid2d.rs          2D wafer grid state (Phase 2), 2D diffusion
│   │   │   ├── grid3d.rs          3D wafer mesh state (Phase 3), 3D diffusion + etch
│   │   │   ├── gds.rs              GDSII binary parser + polygon-to-mask rasterizer
│   │   │   ├── lithography.rs     Mask expose / develop
│   │   │   ├── etch.rs            Anisotropic / isotropic etch + deposition (2D)
│   │   │   └── bin/recipe_runner.rs   CLI: JSON recipe in, wafer JSON out (2D only)
│   │   └── examples/
│   │       ├── locos_2d_demo.rs
│   │       └── gds_3d_demo.rs      GDSII-driven 3D LOCOS + contact trench demo
│   ├── core-wasm/             (stub — wasm-bindgen target, not yet built)
│   └── core-py/               (stub — PyO3 bindings, not yet built)
│
├── backend/                   FastAPI backend
│   ├── api.py                     All routes, middleware, startup lifecycle
│   ├── core/
│   │   ├── auth.py                X-API-Key middleware (timing-safe)
│   │   ├── physics_bench.py       Startup regression check vs. known-good physics
│   │   ├── database.py            SQLAlchemy async ORM + auto-migrations
│   │   ├── logging_config.py      Structured JSON logging
│   │   ├── engine_bridge.py       Subprocess bridge to recipe_runner
│   │   └── recipe_pipeline.py     3-stage: validate → execute → summarize
│   └── models/schemas.py          Pydantic v2 request/response models
│
├── frontend/                  Streamlit UI (alternative)
│   ├── app.py
│   ├── components/sidebar.py
│   ├── pages/1_Run_History.py
│   └── utils/api_client.py · visualizer.py
│
├── frontend_react/            React UI (recommended)
│   ├── index.html · package.json · vite.config.js · start_react.sh/.bat
│   └── src/
│       ├── App.jsx · main.jsx
│       ├── api/client.js
│       ├── components/{Logo,Layout,WaferCanvas,ResultCard}.jsx
│       ├── pages/{Simulate,History,Analytics}.jsx
│       └── styles/globals.css
│
├── web/index.html             Standalone canvas renderer (no backend needed)
├── web/topography_3d.html     WebGL heightfield viewer for Wafer3d exports (Phase 3)
├── data/ · logs/ · reports/ · temp_uploads/    Runtime dirs (auto-created)
│
├── bot.py                     CLI recipe runner (no server required)
├── requirements.txt · .env.example · .gitignore · LICENSE
├── start_backend.sh/.bat · start_frontend.sh/.bat · start_bot.sh/.bat
│
├── docs/
│   ├── architecture.md        Mesh strategy, layer boundaries, IP rationale
│   └── math_references.md     Every formula, cited, with provenance notes
└── tests/                     (integration test scaffolding)
```

## Quick start

**Rust engine (required by everything else):**
```bash
cargo build -p attofab-core --bin recipe_runner
cargo test -p attofab-core
```

**CLI — no server needed:**
```bash
pip install -r requirements.txt
python3 bot.py recipe.json
# or: cat recipe.json | python3 bot.py -
```

**Backend + React UI (recommended):**
```bash
./start_backend.sh        # http://127.0.0.1:8000
./frontend_react/start_react.sh   # http://127.0.0.1:5173
```

**Backend + Streamlit UI (alternative):**
```bash
./start_backend.sh
./start_frontend.sh       # http://127.0.0.1:8501
```

**Standalone renderer, no backend at all:**
```bash
# Open web/index.html in a browser, drop in a wafer JSON export
# (e.g. from `cargo run --example locos_2d_demo -p attofab-core`)
```

**Phase 3 — GDSII-driven 3D (Rust example only, not yet wired into the CLI/backend):**
```bash
cargo run --example gds_3d_demo -p attofab-core
# writes output/gds_3d_demo.json - open web/topography_3d.html and drop it in
```

### Example recipe

```json
{
  "nx": 60, "ny": 80, "dx_um": 0.01, "dy_um": 0.01,
  "substrate": { "dopant": "Boron", "concentration_cm3": 1e15 },
  "steps": [
    { "op": "oxidize", "temperature_c": 1000, "time_hours": 0.75, "ambient": "Dry" },
    { "op": "implant", "dopant": "Phosphorus", "dose_cm2": 1e15, "energy_kev": 80 },
    { "op": "anneal", "temperature_c": 1000, "time_minutes": 20 }
  ]
}
```

See `crates/core/src/bin/recipe_runner.rs` for the full recipe format,
including masked (LOCOS-style) steps, lithography (`spin_photoresist` /
`expose` / `develop`), and etch/deposit.

## Documentation

- [`docs/architecture.md`](docs/architecture.md) — mesh/memory strategy ("start dense, prove correctness"), why each layer is built the way it is, IP/export-control rationale
- [`docs/math_references.md`](docs/math_references.md) — every physics formula used, with citations and provenance notes (including honest documentation of the Pearson-IV limitation)

## Contributing

Physics changes should cite a public source (see the standard
`docs/math_references.md` holds itself to) and include a test — ideally
one that validates against an independent computation (closed-form
formula, a second engine, or a known analytical solution), not just
whatever the code currently outputs. Benchmarking contributions against
the [SkyWater 130nm PDK](https://github.com/google/skywater-pdk) are
especially welcome.

## License

MIT — see [`LICENSE`](LICENSE). AttoFab is a community-driven educational
project and is not affiliated with or representing any semiconductor
foundry.
