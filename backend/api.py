"""
api.py - AttoFab backend: all routes, middleware, startup lifecycle.

Run with: uvicorn api:app --reload --port 8000  (from the backend/ dir)
or:       ./start_backend.sh  (from the repo root)
"""
from __future__ import annotations

import logging
from contextlib import asynccontextmanager

from fastapi import Depends, FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from sqlalchemy.ext.asyncio import AsyncSession

from core.auth import ApiKeyMiddleware
from core.database import get_run, get_session, init_db, list_runs, save_run
from core.logging_config import configure_logging
from core.physics_bench import run_startup_physics_check
from core.recipe_pipeline import RecipeValidationError, run_pipeline
from models.schemas import HealthResponse, RecipeRequest, RunSummary, SimulateResponse

configure_logging()
logger = logging.getLogger(__name__)

_startup_state = {"engine_ok": False, "physics_bench_ok": False}


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("startup_begin")
    await init_db()
    try:
        _startup_state["physics_bench_ok"] = run_startup_physics_check()
        _startup_state["engine_ok"] = True
    except Exception as e:  # engine genuinely unreachable (binary missing, etc.)
        logger.warning("startup_engine_check_failed", extra={"error": str(e)})
        _startup_state["engine_ok"] = False
    logger.info("startup_complete", extra=dict(_startup_state))
    yield
    logger.info("shutdown")


app = FastAPI(
    title="AttoFab API",
    description="Backend for AttoFab, an open-source, education-first electronics fabrication simulator.",
    version="0.1.0",
    lifespan=lifespan,
)

app.add_middleware(ApiKeyMiddleware)
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # local/education tool - tighten if deployed publicly
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/api/health", response_model=HealthResponse)
async def health() -> HealthResponse:
    return HealthResponse(
        status="ok",
        engine_ok=_startup_state["engine_ok"],
        physics_bench_ok=_startup_state["physics_bench_ok"],
    )


@app.post("/api/simulate", response_model=SimulateResponse)
async def simulate(request: RecipeRequest, session: AsyncSession = Depends(get_session)) -> SimulateResponse:
    recipe = request.to_engine_payload()
    try:
        wafer, summary = run_pipeline(recipe)
    except RecipeValidationError as e:
        raise HTTPException(status_code=422, detail=str(e)) from e
    except Exception as e:
        logger.error("simulate_failed", extra={"error": str(e)})
        raise HTTPException(status_code=502, detail=f"engine error: {e}") from e

    run = await save_run(session, nx=request.nx, ny=request.ny, recipe=recipe, summary=summary, wafer=wafer)
    return SimulateResponse(run_id=run.id, summary=summary, wafer=wafer)


@app.get("/api/runs", response_model=list[RunSummary])
async def get_runs(limit: int = 50, session: AsyncSession = Depends(get_session)) -> list[RunSummary]:
    runs = await list_runs(session, limit=limit)
    return [RunSummary(**r.summary()) for r in runs]


@app.get("/api/runs/{run_id}")
async def get_run_detail(run_id: int, session: AsyncSession = Depends(get_session)) -> dict:
    run = await get_run(session, run_id)
    if run is None:
        raise HTTPException(status_code=404, detail=f"run {run_id} not found")
    import json

    return {
        "id": run.id,
        "created_at": run.created_at.isoformat(),
        "recipe": json.loads(run.recipe_json),
        "summary": json.loads(run.summary_json),
        "wafer": json.loads(run.wafer_json),
    }
