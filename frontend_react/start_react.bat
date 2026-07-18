@echo off
cd /d "%~dp0"
if not exist node_modules (
  echo Installing dependencies ^(first run only^)...
  npm install
)
npm run dev
