#!/usr/bin/env bash
# Start the AttoFab React frontend (recommended UI).
set -e
cd "$(dirname "$0")"
if [ ! -d node_modules ]; then
  echo "Installing dependencies (first run only)..."
  npm install
fi
exec npm run dev
