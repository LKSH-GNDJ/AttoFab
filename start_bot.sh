#!/usr/bin/env bash
# Run the AttoFab CLI. Usage: ./start_bot.sh recipe.json
set -e
cd "$(dirname "$0")"
exec python3 bot.py "$@"
