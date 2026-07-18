"""
components/sidebar.py - shared sidebar: logo, backend health, nav hint.
"""
from __future__ import annotations

from pathlib import Path

import streamlit as st

from utils.api_client import health

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


def render_sidebar() -> None:
    logo_path = REPO_ROOT / "attofab_logo.svg"
    if logo_path.exists():
        st.sidebar.image(str(logo_path), use_container_width=True)

    st.sidebar.markdown("---")

    try:
        h = health()
        if h.get("engine_ok") and h.get("physics_bench_ok"):
            st.sidebar.success("Backend online \u00b7 physics check OK")
        elif h.get("engine_ok"):
            st.sidebar.warning("Backend online \u00b7 physics bench FAILED")
        else:
            st.sidebar.error("Backend online \u00b7 engine unreachable")
    except Exception:
        st.sidebar.error("Backend unreachable \u2014 start it with ./start_backend.sh")

    st.sidebar.markdown("---")
    st.sidebar.caption("AttoFab \u00b7 open-source fabrication simulator")
