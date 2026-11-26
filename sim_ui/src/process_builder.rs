//! Process Flow Sidebar — sequential recipe view of the six fabrication
//! stages as a reorderable, toggleable node chain.

use sim_core::physics::Stage;

/// One entry in the process flow sidebar.
#[derive(Debug, Clone)]
pub struct FlowNode {
    pub stage: Stage,
    pub enabled: bool,
    pub expanded: bool,
    /// True if this stage's parameters have changed since it last ran
    /// (mirrors `sim_core::engine::Engine::is_dirty`), shown as a status dot.
    pub dirty: bool,
}

/// The full orderable process flow. Order matters: it defines the recipe
/// sequence, independent of the engine's dependency graph (which tracks
/// data dependencies, not click-order).
#[derive(Debug, Clone)]
pub struct ProcessFlow {
    pub nodes: Vec<FlowNode>,
}

impl ProcessFlow {
    pub fn default_recipe() -> Self {
        Self {
            nodes: Stage::DEFAULT_RECIPE
                .iter()
                .map(|&stage| FlowNode { stage, enabled: true, expanded: false, dirty: false })
                .collect(),
        }
    }

    /// Move a node from `from` to `to` index (drag-and-drop reorder).
    pub fn reorder(&mut self, from: usize, to: usize) {
        if from >= self.nodes.len() || to >= self.nodes.len() {
            return;
        }
        let node = self.nodes.remove(from);
        self.nodes.insert(to, node);
    }

    pub fn toggle_enabled(&mut self, index: usize) {
        if let Some(node) = self.nodes.get_mut(index) {
            node.enabled = !node.enabled;
        }
    }

    pub fn set_expanded(&mut self, index: usize, expanded: bool) {
        if let Some(node) = self.nodes.get_mut(index) {
            node.expanded = expanded;
        }
    }

    pub fn mark_dirty(&mut self, stage: Stage) {
        if let Some(node) = self.nodes.iter_mut().find(|n| n.stage == stage) {
            node.dirty = true;
        }
    }

    pub fn clear_dirty(&mut self, stage: Stage) {
        if let Some(node) = self.nodes.iter_mut().find(|n| n.stage == stage) {
            node.dirty = false;
        }
    }

    /// The active (enabled) stages, in recipe order — what the engine
    /// should actually execute.
    pub fn active_sequence(&self) -> Vec<Stage> {
        self.nodes.iter().filter(|n| n.enabled).map(|n| n.stage).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_recipe_has_all_seven_stages_enabled() {
        let flow = ProcessFlow::default_recipe();
        assert_eq!(flow.nodes.len(), 7);
        assert_eq!(flow.active_sequence().len(), 7);
    }

    #[test]
    fn disabling_a_stage_removes_it_from_active_sequence() {
        let mut flow = ProcessFlow::default_recipe();
        flow.toggle_enabled(2); // Lithography
        assert_eq!(flow.active_sequence().len(), 6);
    }

    #[test]
    fn reorder_moves_node_to_new_position() {
        let mut flow = ProcessFlow::default_recipe();
        let first_stage = flow.nodes[0].stage;
        flow.reorder(0, 3);
        assert_eq!(flow.nodes[3].stage, first_stage);
    }
}
