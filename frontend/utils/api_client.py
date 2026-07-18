"""
utils/api_client.py - thin wrapper over the AttoFab backend HTTP API.
"""
from __future__ import annotations

import os

import requests

BASE_URL = os.environ.get("ATTOFAB_API_BASE_URL", "http://127.0.0.1:8000")


def _headers() -> dict:
    key = os.environ.get("ATTOFAB_API_KEY")
    return {"X-API-Key": key} if key else {}


def health() -> dict:
    r = requests.get(f"{BASE_URL}/api/health", headers=_headers(), timeout=5)
    r.raise_for_status()
    return r.json()


def simulate(recipe: dict) -> dict:
    r = requests.post(f"{BASE_URL}/api/simulate", json=recipe, headers=_headers(), timeout=60)
    r.raise_for_status()
    return r.json()


def list_runs(limit: int = 50) -> list[dict]:
    r = requests.get(f"{BASE_URL}/api/runs", params={"limit": limit}, headers=_headers(), timeout=10)
    r.raise_for_status()
    return r.json()


def get_run(run_id: int) -> dict:
    r = requests.get(f"{BASE_URL}/api/runs/{run_id}", headers=_headers(), timeout=10)
    r.raise_for_status()
    return r.json()
