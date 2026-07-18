"""
app.py - AttoFab Streamlit frontend (alternative to frontend_react/).

Run with: streamlit run app.py  (from the frontend/ dir)
or:       ./start_frontend.sh   (from the repo root)
"""
from __future__ import annotations

import json

import streamlit as st

from components.sidebar import render_sidebar
from utils.api_client import simulate
from utils.visualizer import render_cross_section

st.set_page_config(page_title="AttoFab \u2014 Simulate", page_icon="\U0001f9ea", layout="wide")
render_sidebar()

st.title("AttoFab Simulator")
st.caption("Open-source, education-first electronics fabrication simulator.")

DEFAULT_RECIPE = {
    "nx": 60,
    "ny": 80,
    "dx_um": 0.01,
    "dy_um": 0.01,
    "substrate": {"dopant": "Boron", "concentration_cm3": 1e15},
    "steps": [
        {"op": "oxidize", "temperature_c": 1000, "time_hours": 0.75, "ambient": "Dry"},
        {"op": "implant", "dopant": "Phosphorus", "dose_cm2": 1e15, "energy_kev": 80},
        {"op": "anneal", "temperature_c": 1000, "time_minutes": 20},
    ],
}

col_form, col_result = st.columns([1, 1.4])

with col_form:
    st.subheader("Recipe")
    recipe_text = st.text_area(
        "Process recipe (JSON)",
        value=json.dumps(DEFAULT_RECIPE, indent=2),
        height=420,
    )
    run_clicked = st.button("Run simulation", type="primary")

with col_result:
    st.subheader("Result")
    if run_clicked:
        try:
            recipe = json.loads(recipe_text)
        except json.JSONDecodeError as e:
            st.error(f"Invalid JSON: {e}")
        else:
            with st.spinner("Running through the Rust physics engine\u2026"):
                try:
                    result = simulate(recipe)
                except Exception as e:
                    st.error(f"Simulation failed: {e}")
                else:
                    st.success(f"Run #{result['run_id']} complete")
                    fig = render_cross_section(result["wafer"], title=f"Run #{result['run_id']}")
                    st.pyplot(fig)

                    st.subheader("Summary")
                    st.json(result["summary"])
    else:
        st.info("Edit the recipe and click **Run simulation**.")
