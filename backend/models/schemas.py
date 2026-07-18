"""
models/schemas.py - Pydantic v2 request/response models.

Steps are intentionally permissive (`extra="allow"`, no per-op field
validation) rather than a fully-typed union of every step shape: the Rust
engine's serde-tagged enum is already the strict validator (it will reject
an unknown `op` or missing required field with a clear error), so
duplicating that validation in Python would just be two sources of truth
drifting apart. Python's job here is shape-level sanity (nx/ny/steps
exist and are reasonable), not step-level physics validation.
"""
from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class SubstrateSpec(BaseModel):
    dopant: str = Field(description='e.g. "Boron", "Phosphorus", "Arsenic", "Antimony"')
    concentration_cm3: float = Field(gt=0)


class StepSpec(BaseModel):
    model_config = ConfigDict(extra="allow")
    op: str


class RecipeRequest(BaseModel):
    nx: int = Field(gt=0, le=2000, description="Grid width in voxels")
    ny: int = Field(gt=0, le=2000, description="Grid depth in voxels")
    dx_um: float = Field(gt=0)
    dy_um: float = Field(gt=0)
    substrate: SubstrateSpec
    steps: list[StepSpec] = Field(min_length=1, max_length=200)

    def to_engine_payload(self) -> dict[str, Any]:
        return self.model_dump(mode="json")


class RunSummary(BaseModel):
    id: int
    created_at: str
    nx: int
    ny: int
    summary: dict[str, Any]


class SimulateResponse(BaseModel):
    run_id: int
    summary: dict[str, Any]
    wafer: dict[str, Any]


class HealthResponse(BaseModel):
    status: str
    engine_ok: bool
    physics_bench_ok: bool
