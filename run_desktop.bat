@echo off
setlocal
cd /d "%~dp0"
cargo run -p app_desktop
if errorlevel 1 pause
