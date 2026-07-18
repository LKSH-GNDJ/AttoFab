#!/usr/bin/env python3
"""
bot.py - AttoFab CLI: run a process recipe directly through the Rust
engine, no backend server required.

Usage:
    python bot.py recipe.json
    python bot.py recipe.json --out wafer.json
    cat recipe.json | python bot.py -

Prints a human-readable summary to stdout; optionally writes the full
wafer JSON with --out.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "backend"))

from core.engine_bridge import EngineError, run_recipe  # noqa: E402
from core.recipe_pipeline import RecipeValidationError, summarize, validate  # noqa: E402


def load_recipe(source: str) -> dict:
    if source == "-":
        return json.loads(sys.stdin.read())
    return json.loads(Path(source).read_text())


def print_summary(summary: dict, run_label: str) -> None:
    print(f"\n=== AttoFab run: {run_label} ===")
    print("\nMaterial composition:")
    for mat, pct in summary["material_pct"].items():
        print(f"  {mat:<20} {pct:>6.2f}%")

    ox = summary["oxide_nm"]
    print(f"\nOxide thickness: avg={ox['avg_nm']}nm  max={ox['max_nm']}nm  min={ox['min_nm']}nm")

    if summary["species_peak_cm3"]:
        print("\nPeak dopant concentration:")
        for dopant, peak in summary["species_peak_cm3"].items():
            print(f"  {dopant:<20} {peak} /cm^3")

    print("\nProcess steps:")
    for step in summary["process_steps"]:
        print(f"  - {step}")
    print()


def main() -> int:
    parser = argparse.ArgumentParser(description="Run an AttoFab process recipe from the command line.")
    parser.add_argument("recipe", help='Path to a recipe JSON file, or "-" to read from stdin')
    parser.add_argument("--out", help="Write the full wafer JSON to this path")
    args = parser.parse_args()

    try:
        recipe = load_recipe(args.recipe)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        print(f"error: could not read recipe: {e}", file=sys.stderr)
        return 1

    try:
        validated = validate(recipe)
        wafer = run_recipe(validated)
        summary = summarize(wafer)
    except RecipeValidationError as e:
        print(f"error: invalid recipe: {e}", file=sys.stderr)
        return 1
    except EngineError as e:
        print(f"error: engine failed: {e}", file=sys.stderr)
        return 1

    print_summary(summary, args.recipe)

    if args.out:
        Path(args.out).write_text(json.dumps(wafer))
        print(f"Full wafer state written to {args.out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
