#!/usr/bin/env bash
# Start the AttoFab Streamlit frontend (alternative to frontend_react/).
set -e
cd "$(dirname "$0")/frontend"
exec python3 -m streamlit run app.py
