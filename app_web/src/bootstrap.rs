// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Platform-agnostic session bootstrap logic shared by the WASM entry
//! point in `lib.rs`. Kept free of `wasm-bindgen`/`web-sys` so it can be
//! unit-tested on a native host target instead of requiring a browser or a
//! `wasm32-unknown-unknown` toolchain just to validate the session-state
//! wiring.

use sim_core::engine::Engine;
use sim_core::physics::Stage;
use sim_core::wafer::WaferState;
use sim_ui::process_builder::ProcessFlow;
use sim_ui::table_view::ParameterTable;

/// Everything a running AttoFab session needs, independent of which
/// platform shell (desktop/web/android) is hosting it.
pub struct Session {
    pub engine: Engine,
    pub flow: ProcessFlow,
    pub table: ParameterTable,
    pub wafer: WaferState,
}

impl Session {
    /// Construct a fresh session with the default recipe, default
    /// parameter table, and a blank 300mm silicon substrate. The web shell
    /// calls this once from its `#[wasm_bindgen(start)]` entry point;
    /// native/android call it from their own bootstraps.
    pub fn new_default() -> Self {
        Self {
            engine: Engine::new_default(),
            flow: ProcessFlow::default_recipe(),
            table: ParameterTable::default_table(),
            wafer: WaferState::new_substrate(300.0, 500.0, 10.0),
        }
    }

    /// Mark a stage dirty in both the engine's dependency graph and the
    /// process-flow sidebar's per-node status dot, keeping the two in sync
    /// the way a UI event handler must every time a table cell changes.
    pub fn mark_dirty(&mut self, stage: Stage) {
        self.engine.mark_dirty(stage);
        for s in self.engine.dirty_stages() {
            self.flow.mark_dirty(s);
        }
    }
}

impl Default for Session {
    fn default() -> Self {
        Self::new_default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_session_has_no_dirty_stages() {
        let session = Session::new_default();
        assert!(session.engine.dirty_stages().is_empty());
    }

    #[test]
    fn marking_lithography_dirty_flags_its_sidebar_node_too() {
        let mut session = Session::new_default();
        session.mark_dirty(Stage::Lithography);
        let node = session
            .flow
            .nodes
            .iter()
            .find(|n| n.stage == Stage::Lithography)
            .unwrap();
        assert!(node.dirty);
    }
}
