// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Which stages a move invalidates (PRD G-10).
//!
//! Dragging a stage is a gesture, and a gesture cannot be tested without a pointer.
//! What *can* be tested is the rule the drag exists to serve: moving a stage
//! invalidates every stage downstream of the new position, because each of them was
//! run against a stack the moved stage no longer leaves behind.
//!
//! The rule is deliberately conservative. A stage that is *upstream* of the move
//! keeps its result; a stage at or after the earliest changed position does not.

use crate::StageId;

/// What a move would do, before it is made.
#[derive(Debug, Clone, PartialEq)]
pub struct MovePlan {
    /// The order after the move.
    pub order: Vec<StageId>,
    /// The stages whose results would no longer be valid.
    pub invalidated: Vec<StageId>,
    /// Whether the move changes anything at all.
    pub changed: bool,
}

/// Plan a move of the stage at `from` to `to`.
///
/// Returns `None` when either index is out of range: a drag that lands nowhere is a
/// cancelled gesture, not a move to the end.
pub fn plan(order: &[StageId], from: usize, to: usize) -> Option<MovePlan> {
    if from >= order.len() || to >= order.len() {
        return None;
    }
    let mut next = order.to_vec();
    let moved = next.remove(from);
    next.insert(to, moved);
    let changed = from != to;
    // everything from the earlier of the two positions onwards was run against a
    // stack that this move changes
    let first = from.min(to);
    let invalidated = if changed { next[first..].to_vec() } else { Vec::new() };
    Some(MovePlan { order: next, invalidated, changed })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn order() -> Vec<StageId> {
        StageId::ALL.to_vec()
    }

    #[test]
    fn moving_a_stage_down_invalidates_from_the_new_position_on() {
        let o = order();
        let p = plan(&o, 1, 3).unwrap();
        assert!(p.changed);
        assert_eq!(p.order[3], o[1], "the stage landed where it was dragged");
        assert_eq!(p.invalidated, p.order[1..].to_vec(), "everything from the earlier slot");
        assert!(!p.invalidated.contains(&o[0]), "the stage before the move is untouched");
    }

    #[test]
    fn moving_a_stage_up_invalidates_from_where_it_landed() {
        let o = order();
        let p = plan(&o, 4, 1).unwrap();
        assert_eq!(p.order[1], o[4]);
        assert_eq!(p.invalidated.first(), Some(&o[4]), "the moved stage itself is stale");
        assert_eq!(p.invalidated.len(), o.len() - 1);
    }

    #[test]
    fn a_move_to_the_same_place_invalidates_nothing() {
        let o = order();
        let p = plan(&o, 2, 2).unwrap();
        assert!(!p.changed);
        assert!(p.invalidated.is_empty(), "nothing moved, nothing is stale");
        assert_eq!(p.order, o);
    }

    #[test]
    fn a_move_off_the_end_is_refused() {
        let o = order();
        assert!(plan(&o, 0, o.len()).is_none());
        assert!(plan(&o, o.len(), 0).is_none());
        assert!(plan(&[], 0, 0).is_none());
    }

    #[test]
    fn the_plan_never_loses_or_duplicates_a_stage() {
        let o = order();
        for from in 0..o.len() {
            for to in 0..o.len() {
                let p = plan(&o, from, to).unwrap();
                let mut sorted = p.order.clone();
                sorted.sort_by_key(|s| *s as u8);
                let mut orig = o.clone();
                orig.sort_by_key(|s| *s as u8);
                assert_eq!(sorted, orig, "a move changed the set of stages");
            }
        }
    }
}
