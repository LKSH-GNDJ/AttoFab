"""
pages/1_Run_History.py - browse past simulation runs.
Streamlit auto-discovers files in pages/ and adds them to the sidebar nav.
"""
from __future__ import annotations

import sys
from pathlib import Path

import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from components.sidebar import render_sidebar  # noqa: E402
from utils.api_client import get_run, list_runs  # noqa: E402
from utils.visualizer import render_cross_section  # noqa: E402

st.set_page_config(page_title="AttoFab \u2014 Run History", page_icon="\U0001f4dc", layout="wide")
render_sidebar()

st.title("Run History")

try:
    runs = list_runs()
except Exception as e:
    st.error(f"Could not reach backend: {e}")
    st.stop()

if not runs:
    st.info("No runs yet \u2014 go to the Simulate page and run a recipe.")
    st.stop()

for run in runs:
    with st.expander(f"Run #{run['id']} \u2014 {run['created_at']} \u2014 {run['nx']}x{run['ny']} voxels"):
        st.json(run["summary"])
        if st.button("Load cross-section", key=f"load-{run['id']}"):
            detail = get_run(run["id"])
            fig = render_cross_section(detail["wafer"], title=f"Run #{run['id']}")
            st.pyplot(fig)
