#!/usr/bin/env bash
set -euo pipefail
cargo test -p sim_core recovers_second_order_sequence
cargo bench -p sim_core --bench benchmarks
