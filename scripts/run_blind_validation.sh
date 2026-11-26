#!/usr/bin/env bash
set -euo pipefail
cargo run --release -p sim_core --bin blind_validator -- \
  data/blind_test/demo_synthetic.json validation_results/demo_synthetic
