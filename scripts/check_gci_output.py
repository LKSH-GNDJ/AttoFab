#!/usr/bin/env python3
"""Validate the gci_runner output contract on a synthetic unit-test sequence.

This checks runner wiring/schema and a known manufactured convergence series;
it is not a physical-model validation or a release accuracy claim.
"""
import json
import math
import sys
from pathlib import Path

if len(sys.argv) != 2:
    raise SystemExit("usage: check_gci_output.py <output.json>")
path = Path(sys.argv[1])
data = json.loads(path.read_text(encoding="utf-8"))
report = data.get("report")
required = {"observed_order", "extrapolated_value", "fine_gci_percent", "medium_gci_percent", "monotonic_convergence"}
if not isinstance(report, dict) or not required.issubset(report):
    raise SystemExit(f"invalid GCI output schema; missing {sorted(required - set(report or {}))}")
for key in required - {"monotonic_convergence"}:
    if not isinstance(report[key], (int, float)) or not math.isfinite(report[key]):
        raise SystemExit(f"{key} must be a finite number")
if not isinstance(report["monotonic_convergence"], bool):
    raise SystemExit("monotonic_convergence must be boolean")
if abs(report["observed_order"] - 2.0) > 1e-6:
    raise SystemExit(f"manufactured second-order series failed: p={report['observed_order']}")
if abs(report["extrapolated_value"] - 1.0) > 1e-6:
    raise SystemExit(f"manufactured extrapolation failed: value={report['extrapolated_value']}")
if report["fine_gci_percent"] < 0 or report["medium_gci_percent"] < 0:
    raise SystemExit("GCI percentages must be non-negative")
print("GCI runner schema and manufactured second-order sequence: PASS (synthetic only)")
