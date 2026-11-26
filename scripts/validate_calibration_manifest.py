#!/usr/bin/env python3
"""Verify calibration-file identity against manifest; does not validate physics."""
import hashlib
import json
import sys
from pathlib import Path

root = Path(__file__).resolve().parents[1]
manifest_path = root / "data/calibration/manifest.json"
manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
if manifest.get("schema_version") != 1 or not isinstance(manifest.get("files"), list):
    raise SystemExit("unsupported calibration manifest schema")
errors = []
for entry in manifest["files"]:
    rel = entry.get("path")
    expected = entry.get("sha256")
    if not rel or not expected:
        errors.append(f"entry missing path/hash: {entry!r}")
        continue
    path = root / rel
    if not path.is_file():
        errors.append(f"missing calibration file: {rel}")
        continue
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual != expected:
        errors.append(f"hash mismatch: {rel}")
    if entry.get("evidence_status") not in {"unverified", "traceable_reference", "experimental"}:
        errors.append(f"invalid evidence_status for {rel}")
    if entry.get("evidence_status") != "unverified" and not entry.get("source"):
        errors.append(f"traceable evidence requires source: {rel}")
if errors:
    print("\n".join(errors), file=sys.stderr)
    raise SystemExit(1)
print(f"Calibration manifest identity checks: PASS ({len(manifest['files'])} files); physical provenance/accuracy not established")
