"""
core/recipe_pipeline.py - the 3-stage pipeline every simulation request
goes through, mirroring the shape of AttoSense's 3-stage NLU pipeline but
for physics recipes instead of natural-language intents:

  1. validate  - shape/sanity-check the recipe (Pydantic already handled
                 field-level validation at the API boundary; this stage
                 catches cross-field issues Pydantic can't, e.g. an empty
                 mask range)
  2. execute   - run it through the real Rust engine (core/engine_bridge)
  3. summarize - reduce the (potentially large) wafer grid into a compact,
                 human-readable summary: material composition, oxide
                 stats, dopant peak concentrations - the numbers a caller
                 actually wants at a glance, without shipping the full
                 nx*ny grid every time.
"""
from __future__ import annotations

import logging
from typing import Any

from core.engine_bridge import EngineError, run_recipe

logger = logging.getLogger(__name__)


class RecipeValidationError(ValueError):
    pass


def validate(recipe: dict[str, Any]) -> dict[str, Any]:
    nx = recipe["nx"]
    steps = recipe["steps"]

    for i, step in enumerate(steps):
        mask = step.get("mask")
        if isinstance(mask, dict) and "range" in mask:
            lo, hi = mask["range"]
            if not (0 <= lo < hi <= nx):
                raise RecipeValidationError(
                    f"step {i} ({step.get('op')}): mask range [{lo}, {hi}) is out of bounds for nx={nx}"
                )
        elif isinstance(mask, list) and len(mask) not in (0, nx):
            raise RecipeValidationError(
                f"step {i} ({step.get('op')}): explicit mask has length {len(mask)}, expected {nx}"
            )

    return recipe


def execute(recipe: dict[str, Any]) -> dict[str, Any]:
    try:
        return run_recipe(recipe)
    except EngineError as e:
        logger.error("recipe_execution_failed", extra={"error": str(e)})
        raise


def summarize(wafer: dict[str, Any]) -> dict[str, Any]:
    material = wafer["material"]
    n = len(material)
    counts: dict[str, int] = {}
    for m in material:
        counts[m] = counts.get(m, 0) + 1
    material_pct = {m: round(100 * c / n, 2) for m, c in counts.items()}

    oxide = wafer.get("oxide_thickness_um", [])
    oxide_stats = {
        "avg_nm": round(1000 * sum(oxide) / len(oxide), 2) if oxide else 0.0,
        "max_nm": round(1000 * max(oxide), 2) if oxide else 0.0,
        "min_nm": round(1000 * min(oxide), 2) if oxide else 0.0,
    }

    species_peaks = {}
    for dopant, values in wafer.get("species", {}).items():
        if values:
            species_peaks[dopant] = f"{max(values):.3e}"

    return {
        "material_pct": material_pct,
        "oxide_nm": oxide_stats,
        "species_peak_cm3": species_peaks,
        "process_steps": wafer.get("process_log", []),
    }


def run_pipeline(recipe: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    """Runs all three stages. Returns (wafer, summary)."""
    validated = validate(recipe)
    wafer = execute(validated)
    summary = summarize(wafer)
    return wafer, summary
