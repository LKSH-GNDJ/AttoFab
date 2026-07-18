@echo off
cd /d "%~dp0"

if not exist target\debug\recipe_runner.exe if not exist target\release\recipe_runner.exe (
  echo Building recipe_runner ^(first run only^)...
  cargo build -p attofab-core --bin recipe_runner
)

cd backend
python -m uvicorn api:app --host 127.0.0.1 --port 8000 --reload
