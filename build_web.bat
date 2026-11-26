@echo off
setlocal
cd /d "%~dp0"
rustup target add wasm32-unknown-unknown
where wasm-pack >nul 2>nul
if errorlevel 1 cargo install wasm-pack
wasm-pack build app_web --target web --out-dir pkg
if errorlevel 1 pause
