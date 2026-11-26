// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Platform-agnostic session bootstrap for the Android shell. Mirrors
//! `app_web::bootstrap` deliberately: both platform shells wrap the same
//! `sim_core`/`sim_ui` session state, and keeping the wiring logic free of
//! `android_activity`/NDK types means it can be unit-tested on the host
//! target without an Android emulator or NDK toolchain.

use sim_core::engine::Engine;
use sim_core::physics::Stage;
use sim_core::wafer::WaferState;
use sim_ui::process_builder::ProcessFlow;
use sim_ui::table_view::ParameterTable;

/// Everything a running AttoFab session needs, independent of platform.
pub struct Session {
    pub engine: Engine,
    pub flow: ProcessFlow,
    pub table: ParameterTable,
    pub wafer: WaferState,
}

impl Session {
    pub fn new_default() -> Self {
        Self {
            engine: Engine::new_default(),
            flow: ProcessFlow::default_recipe(),
            table: ParameterTable::default_table(),
            wafer: WaferState::new_substrate(300.0, 500.0, 10.0),
        }
    }

    pub fn mark_dirty(&mut self, stage: Stage) {
        self.engine.mark_dirty(stage);
        for s in self.engine.dirty_stages() {
            self.flow.mark_dirty(s);
        }
    }

    /// A short human-readable status line, used both for the Android
    /// shell's on-screen debug overlay and for `android_logger` output.
    pub fn status_line(&self) -> String {
        format!(
            "AttoFab session: {} stages, {} SVO node(s)",
            self.flow.nodes.len(),
            self.wafer.node_count()
        )
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
    fn status_line_reports_stage_and_node_counts() {
        let session = Session::new_default();
        let line = session.status_line();
        assert!(line.contains("stages"));
        assert!(line.contains("SVO node"));
    }
}
