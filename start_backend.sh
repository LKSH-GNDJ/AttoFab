#!/usr/bin/env bash
# Start the AttoFab FastAPI backend.
set -e
cd "$(dirname "$0")"

if [ ! -f target/debug/recipe_runner ] && [ ! -f target/release/recipe_runner ]; then
  echo "Building recipe_runner (first run only)..."
  cargo build -p attofab-core --bin recipe_runner
fi

cd backend
exec python3 -m uvicorn api:app --host "${ATTOFAB_BACKEND_HOST:-127.0.0.1}" --port "${ATTOFAB_BACKEND_PORT:-8000}" --reload
