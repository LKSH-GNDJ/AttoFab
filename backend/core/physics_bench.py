"""
core/physics_bench.py - validates the engine's output against known
ground truth on startup.

This plays the role AttoSense's calibration.py plays (establishing how
much to trust a model's output) but adapted to what actually needs
calibrating here: not a statistical confidence score, but a correctness
check against textbook physics. If the compiled engine ever drifts from
its regression-tested values (a bad rebuild, a stale binary, a future
change that silently breaks Deal-Grove), this catches it at startup
instead of silently serving wrong physics to users.
"""
from __future__ import annotations

import logging

from core.engine_bridge import EngineError, run_recipe

logger = logging.getLogger(__name__)

# A minimal, fast recipe with a known, previously-verified result: 1000C
# dry oxidation for 45 minutes on a bare column should yield ~29.8nm
# (verified against the Rust wafer1d::tests regression benchmark - see
# docs/math_references.md for the corrected Deal-Grove constants this is
# pinned to. This value legitimately changed once already, when the
# constants table was corrected against a better-sourced reference - this
# check caught the drift immediately on the next startup, which is exactly
# what it's for.)
_BENCHMARK_RECIPE = {
    "nx": 2,
    "ny": 50,
    "dx_um": 0.01,
    "dy_um": 0.005,
    "substrate": {"dopant": "Boron", "concentration_cm3": 1e15},
    "steps": [
        {"op": "oxidize", "temperature_c": 1000.0, "time_hours": 0.75, "ambient": "Dry", "mask": None},
    ],
}
_EXPECTED_OXIDE_NM = 29.8
_TOLERANCE_NM = 1.0


def run_startup_physics_check() -> bool:
    """Returns True if the engine's output matches known ground truth
    within tolerance. Logs a warning (does not raise) on mismatch, since a
    physics regression shouldn't take the whole API down - it should be
    loud in the logs and visible via /api/health."""
    try:
        wafer = run_recipe(_BENCHMARK_RECIPE)
    except EngineError as e:
        logger.warning("physics_bench_engine_error", extra={"error": str(e)})
        return False

    oxide_nm = wafer["oxide_thickness_um"][0] * 1000.0
    delta = abs(oxide_nm - _EXPECTED_OXIDE_NM)
    ok = delta <= _TOLERANCE_NM

    log = logger.info if ok else logger.warning
    log(
        "physics_bench_result",
        extra={"expected_nm": _EXPECTED_OXIDE_NM, "actual_nm": oxide_nm, "delta_nm": delta, "ok": ok},
    )
    return ok
