#!/usr/bin/env python3
"""Structural validation for calibration JSON files; not scientific validation."""
import json
import sys
from pathlib import Path

def load(path):
    with path.open(encoding="utf-8") as f:
        return json.load(f)

def finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and value == value and abs(value) != float("inf")

def validate_file(path):
    data = load(path)
    name = path.as_posix().lower()
    if "oxidation" in name:
        # Accept either the compact dry/wet schema or the PRD's parameter-pack schema.
        if "dry" in data or "wet" in data:
            for ambient in ("dry", "wet"):
                if ambient not in data: raise ValueError(f"missing {ambient} model")
                for key in ("a_um", "b_um2_hr", "ea_eV", "eb_eV"):
                    if key not in data[ambient] or not finite_number(data[ambient][key]): raise ValueError(f"{ambient}.{key} must be numeric")
        elif not isinstance(data, dict) or not data:
            raise ValueError("oxidation calibration must be a non-empty object")
    elif "implant" in name:
        if not isinstance(data, dict): raise ValueError("implant calibration must be an object")
        for key in ("species", "energies_keV", "range_nm", "straggle_nm"):
            if key not in data: raise ValueError(f"missing {key}")
        for key in ("energies_keV", "range_nm", "straggle_nm"):
            if not isinstance(data[key], list) or not all(finite_number(x) for x in data[key]): raise ValueError(f"{key} must be a numeric list")
        if len({len(data[k]) for k in ("energies_keV", "range_nm", "straggle_nm")}) != 1: raise ValueError("energy/range/straggle lengths differ")
    elif not isinstance(data, (dict, list)):
        raise ValueError("calibration root must be an object or list")
    print(f"OK (schema only): {path}")

def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/calibration")
    if not root.exists():
        print(f"ERROR: calibration path does not exist: {root}", file=sys.stderr); return 2
    files = sorted(root.rglob("*.json")) if root.is_dir() else [root]
    if not files:
        print(f"ERROR: no JSON calibration files found under {root}", file=sys.stderr); return 2
    failed = 0
    for path in files:
        try: validate_file(path)
        except Exception as exc:
            failed += 1; print(f"FAIL: {path}: {exc}", file=sys.stderr)
    print(f"Structural validation: {len(files)-failed}/{len(files)} passed; this does not validate physical accuracy.")
    return 1 if failed else 0

if __name__ == "__main__": raise SystemExit(main())
