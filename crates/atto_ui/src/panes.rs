// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! One, two and four viewport panes over the same state (PRD G-08).
//!
//! The panes show one session from different angles, so the *only* thing that
//! differs between them is the rectangle each occupies. That makes the layout a
//! pure function of the window rectangle and the pane count — and it makes the two
//! things that can go wrong checkable: a gap between panes, and a pane with no area.

/// A rectangle in window coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rect {
    pub x: f64,
    pub y: f64,
    pub w: f64,
    pub h: f64,
}

impl Rect {
    /// The area, which must be positive for a pane to be drawable.
    pub fn area(&self) -> f64 {
        self.w * self.h
    }

    /// Whether two rectangles overlap.
    pub fn overlaps(&self, o: &Rect) -> bool {
        self.x < o.x + o.w && o.x < self.x + self.w && self.y < o.y + o.h && o.y < self.y + self.h
    }
}

/// How many panes the viewport is showing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PaneCount {
    One,
    Two,
    Four,
}

impl PaneCount {
    /// The number of panes.
    pub fn n(self) -> usize {
        match self {
            PaneCount::One => 1,
            PaneCount::Two => 2,
            PaneCount::Four => 4,
        }
    }
}

/// Divide a window into panes.
///
/// Two panes split side by side; four split into quadrants. The division is exact:
/// the panes tile the window with no gap and no overlap, which is what stops a
/// one-pixel seam appearing between viewports at odd window sizes.
pub fn layout(window: Rect, count: PaneCount) -> Vec<Rect> {
    let (x, y, w, h) = (window.x, window.y, window.w, window.h);
    match count {
        PaneCount::One => vec![window],
        PaneCount::Two => vec![
            Rect { x, y, w: w * 0.5, h },
            Rect { x: x + w * 0.5, y, w: w - w * 0.5, h },
        ],
        PaneCount::Four => {
            let hw = w * 0.5;
            let hh = h * 0.5;
            vec![
                Rect { x, y, w: hw, h: hh },
                Rect { x: x + hw, y, w: w - hw, h: hh },
                Rect { x, y: y + hh, w: hw, h: h - hh },
                Rect { x: x + hw, y: y + hh, w: w - hw, h: h - hh },
            ]
        }
    }
}

/// Check a layout: every pane has area, nothing overlaps, and the panes cover the
/// window. Returns every violation.
pub fn check(window: Rect, panes: &[Rect]) -> Vec<String> {
    let mut out = Vec::new();
    for (i, p) in panes.iter().enumerate() {
        if !(p.area() > 0.0) {
            out.push(format!("pane {i} has no area"));
        }
    }
    for i in 0..panes.len() {
        for j in i + 1..panes.len() {
            if panes[i].overlaps(&panes[j]) {
                out.push(format!("panes {i} and {j} overlap"));
            }
        }
    }
    let covered: f64 = panes.iter().map(|p| p.area()).sum();
    let total = window.area();
    if total > 0.0 && (covered - total).abs() / total > 1e-9 {
        out.push(format!("the panes cover {covered:.3} of {total:.3}"));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn window() -> Rect {
        Rect { x: 0.0, y: 0.0, w: 1280.0, h: 720.0 }
    }

    #[test]
    fn every_count_tiles_the_window_exactly() {
        for c in [PaneCount::One, PaneCount::Two, PaneCount::Four] {
            let panes = layout(window(), c);
            assert_eq!(panes.len(), c.n());
            assert!(check(window(), &panes).is_empty(), "{c:?}: {:?}", check(window(), &panes));
        }
    }

    #[test]
    fn an_odd_window_size_leaves_no_seam() {
        // halves of an odd width must still sum to the whole, or a seam appears
        let w = Rect { x: 0.0, y: 0.0, w: 1281.0, h: 721.0 };
        for c in [PaneCount::Two, PaneCount::Four] {
            let panes = layout(w, c);
            assert!(check(w, &panes).is_empty(), "{c:?}: {:?}", check(w, &panes));
        }
    }

    #[test]
    fn four_panes_are_quadrants() {
        let panes = layout(window(), PaneCount::Four);
        assert_eq!(panes[0], Rect { x: 0.0, y: 0.0, w: 640.0, h: 360.0 });
        assert_eq!(panes[3], Rect { x: 640.0, y: 360.0, w: 640.0, h: 360.0 });
    }

    #[test]
    fn a_window_with_no_area_is_reported_rather_than_tiled() {
        let tiny = Rect { x: 0.0, y: 0.0, w: 0.0, h: 100.0 };
        let panes = layout(tiny, PaneCount::Two);
        assert!(!check(tiny, &panes).is_empty(), "zero width cannot be tiled");
    }

    #[test]
    fn an_overlap_is_detected() {
        let a = Rect { x: 0.0, y: 0.0, w: 10.0, h: 10.0 };
        let b = Rect { x: 5.0, y: 0.0, w: 10.0, h: 10.0 };
        assert!(a.overlaps(&b));
        assert!(check(Rect { x: 0.0, y: 0.0, w: 20.0, h: 10.0 }, &[a, b])
            .iter()
            .any(|m| m.contains("overlap")));
    }
}
