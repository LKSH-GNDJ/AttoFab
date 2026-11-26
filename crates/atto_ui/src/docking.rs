// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Undocking a panel into a floating window (PRD G-28).
//!
//! A panel is either docked or floating, never both and never neither — that is the
//! invariant worth enforcing, because the failure mode is a panel that vanishes when
//! the window it floated over is closed. Docking again restores the slot it came
//! from, so an undock is always reversible.

/// Where a panel is.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Placement {
    /// In its slot in the main window.
    Docked,
    /// In a window of its own, at a position and size.
    Floating { x: f64, y: f64, w: f64, h: f64 },
}

/// A panel's placement and the slot it belongs to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Panel {
    placement: Placement,
    /// The docked size, kept while floating so redocking restores it.
    docked_size: (f64, f64),
    /// The screen it may float within.
    screen: (f64, f64),
}

impl Panel {
    /// A docked panel of a given size on a screen of a given size.
    pub fn new(w: f64, h: f64, screen: (f64, f64)) -> Self {
        Panel { placement: Placement::Docked, docked_size: (w, h), screen }
    }

    /// Where the panel is.
    pub fn placement(&self) -> Placement {
        self.placement
    }

    /// Whether it is floating.
    pub fn is_floating(&self) -> bool {
        matches!(self.placement, Placement::Floating { .. })
    }

    /// Float the panel at a position, clamped onto the screen.
    ///
    /// A panel dragged past the edge is pulled back rather than allowed to hang off
    /// it: a window whose title bar is off-screen cannot be recovered by dragging.
    pub fn undock(&mut self, x: f64, y: f64) {
        let (w, h) = self.docked_size;
        let (sw, sh) = self.screen;
        let cx = x.max(0.0).min((sw - w).max(0.0));
        let cy = y.max(0.0).min((sh - h).max(0.0));
        self.placement = Placement::Floating { x: cx, y: cy, w, h };
    }

    /// Put the panel back in its slot.
    pub fn dock(&mut self) {
        self.placement = Placement::Docked;
    }

    /// The panel's current rectangle.
    pub fn rect(&self) -> (f64, f64, f64, f64) {
        match self.placement {
            Placement::Docked => (0.0, 0.0, self.docked_size.0, self.docked_size.1),
            Placement::Floating { x, y, w, h } => (x, y, w, h),
        }
    }

    /// Check the panel's invariants, returning every violation.
    pub fn check(&self) -> Vec<String> {
        let mut out = Vec::new();
        let (_, _, w, h) = self.rect();
        if !(w > 0.0) || !(h > 0.0) {
            out.push(format!("the panel has no area: {w:.0}x{h:.0}"));
        }
        if let Placement::Floating { x, y, .. } = self.placement {
            let (sw, sh) = self.screen;
            if x < 0.0 || y < 0.0 {
                out.push(format!("the floating panel is off the screen at ({x:.0}, {y:.0})"));
            }
            if x + w > sw + 1e-9 || y + h > sh + 1e-9 {
                out.push(format!("the floating panel hangs off a {sw:.0}x{sh:.0} screen"));
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn panel() -> Panel {
        Panel::new(320.0, 600.0, (1280.0, 720.0))
    }

    #[test]
    fn a_panel_starts_docked_and_round_trips() {
        let mut p = panel();
        assert_eq!(p.placement(), Placement::Docked);
        p.undock(100.0, 80.0);
        assert!(p.is_floating());
        p.dock();
        assert_eq!(p.placement(), Placement::Docked);
        assert_eq!(p.rect(), (0.0, 0.0, 320.0, 600.0), "the slot is restored");
    }

    #[test]
    fn a_panel_dragged_off_the_screen_is_pulled_back() {
        // a window whose title bar is off-screen cannot be dragged back
        let mut p = panel();
        p.undock(5000.0, 5000.0);
        assert!(p.check().is_empty(), "{:?}", p.check());
        let (x, y, _, _) = p.rect();
        assert!(x <= 1280.0 - 320.0, "x stayed on screen: {x}");
        assert!(y <= 720.0 - 600.0, "y stayed on screen: {y}");
    }

    #[test]
    fn a_negative_position_is_pulled_back_too() {
        let mut p = panel();
        p.undock(-50.0, -20.0);
        assert!(p.check().is_empty(), "{:?}", p.check());
        assert_eq!(p.rect().0, 0.0);
    }

    #[test]
    fn a_panel_is_never_both_docked_and_floating() {
        let mut p = panel();
        p.undock(10.0, 10.0);
        let floating = p.is_floating();
        p.dock();
        assert!(floating && !p.is_floating(), "the states are exclusive");
    }

    #[test]
    fn a_panel_with_no_area_is_reported() {
        let p = Panel::new(0.0, 600.0, (1280.0, 720.0));
        assert!(!p.check().is_empty());
    }
}
