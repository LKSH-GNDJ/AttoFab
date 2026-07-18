"""
core/engine_bridge.py - calls the real Rust physics engine.

Deliberately a subprocess boundary (JSON in via stdin, JSON out via
stdout) rather than a compiled PyO3 extension: any Python environment
that can spawn a process is compatible, with zero coupling to the Rust
toolchain version used to build the backend's Python dependencies. This
mirrors AttoSense's multimodal.py role - a bridge to a heavier external
computation - just swapping a hosted API call for a local compiled binary.
"""
from __future__ import annotations

import json
import logging
import os
import subprocess
from pathlib import Path

logger = logging.getLogger(__name__)

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


class EngineError(RuntimeError):
    pass


def _find_binary() -> Path:
    override = os.environ.get("ATTOFAB_ENGINE_BINARY")
    if override:
        p = Path(override)
        if p.exists():
            return p
        raise EngineError(f"ATTOFAB_ENGINE_BINARY set to {p}, but no such file exists")

    candidates = [
        REPO_ROOT / "target" / "release" / "recipe_runner",
        REPO_ROOT / "target" / "debug" / "recipe_runner",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise EngineError(
        "recipe_runner binary not found. Build it with "
        "`cargo build -p attofab-core --bin recipe_runner` "
        f"(looked in: {', '.join(str(c) for c in candidates)})"
    )


def run_recipe(recipe: dict, timeout_seconds: float = 30.0) -> dict:
    """Run a process recipe through the real Rust engine and return the
    resulting wafer state as a dict (parsed from the engine's JSON output).

    Raises EngineError if the binary is missing, times out, exits non-zero,
    or produces output that isn't valid JSON.
    """
    binary = _find_binary()
    payload = json.dumps(recipe)

    try:
        result = subprocess.run(
            [str(binary)],
            input=payload,
            capture_output=True,
            text=True,
            timeout=timeout_seconds,
        )
    except subprocess.TimeoutExpired as e:
        raise EngineError(f"engine timed out after {timeout_seconds}s") from e

    if result.returncode != 0:
        logger.error("engine_run_failed", extra={"stderr": result.stderr, "returncode": result.returncode})
        raise EngineError(f"engine exited with code {result.returncode}: {result.stderr.strip()}")

    for line in result.stderr.splitlines():
        logger.info("engine_step", extra={"engine_log": line})

    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError as e:
        raise EngineError(f"engine produced invalid JSON: {e}") from e
