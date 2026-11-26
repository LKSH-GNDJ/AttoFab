$ErrorActionPreference = 'Stop'
cargo run --release -p sim_core --bin blind_validator -- data/blind_test/demo_synthetic.json validation_results/demo_synthetic
if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
