// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The dockable layout: recipe, viewport, inspector and bottom dock (PRD G-01).
//!
//! The layout is a tree of splits rather than a table of rectangles, because that is
//! what makes it dockable: moving a panel means changing one split, and the rest of
//! the window follows. The invariant that has to hold at every step is that the four
//! areas tile the window exactly — no gap, no overlap, and each one present once.

use crate::panes::Rect;

/// A region of the window.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Area {
    /// The recipe list.
    Recipe,
    /// The 3-D viewport.
    Viewport,
    /// The inspector.
    Inspector,
    /// The timeline and console.
    BottomDock,
}

impl Area {
    /// All four, in reading order.
    pub const ALL: [Area; 4] = [Area::Recipe, Area::Viewport, Area::Inspector, Area::BottomDock];

    /// The name a panel shows.
    pub fn label(self) -> &'static str {
        match self {
            Area::Recipe => "Recipe",
            Area::Viewport => "Viewport",
            Area::Inspector => "Inspector",
            Area::BottomDock => "Bottom dock",
        }
    }
}

/// Where each area sits, given the window.
#[derive(Debug, Clone, PartialEq)]
pub struct Layout {
    /// `(area, rect)` in draw order.
    pub areas: Vec<(Area, Rect)>,
}

impl Layout {
    /// The default arrangement: recipe left, inspector right, dock along the bottom,
    /// viewport in what is left.
    ///
    /// The fractions are of the window, so the layout scales with it; the bottom dock
    /// takes a fixed share rather than a fixed height, for the same reason.
    pub fn default_for(window: Rect) -> Self {
        let (x, y, w, h) = (window.x, window.y, window.w, window.h);
        let left = w * 0.20;
        let right = w * 0.22;
        let bottom = h * 0.25;
        let mid_w = w - left - right;
        let top_h = h - bottom;
        Layout {
            areas: vec![
                (Area::Recipe, Rect { x, y, w: left, h: top_h }),
                (Area::Viewport, Rect { x: x + left, y, w: mid_w, h: top_h }),
                (Area::Inspector, Rect { x: x + left + mid_w, y, w: right, h: top_h }),
                (Area::BottomDock, Rect { x, y: y + top_h, w, h: bottom }),
            ],
        }
    }

    /// The rectangle an area occupies.
    pub fn rect(&self, area: Area) -> Option<Rect> {
        self.areas.iter().find(|(a, _)| *a == area).map(|(_, r)| *r)
    }

    /// Check the layout, returning every violation.
    pub fn check(&self, window: Rect) -> Vec<String> {
        let mut out = Vec::new();
        for a in Area::ALL {
            let n = self.areas.iter().filter(|(b, _)| *b == a).count();
            if n == 0 {
                out.push(format!("{} is missing", a.label()));
            } else if n > 1 {
                out.push(format!("{} appears {n} times", a.label()));
            }
        }
        for i in 0..self.areas.len() {
            for j in i + 1..self.areas.len() {
                if self.areas[i].1.overlaps(&self.areas[j].1) {
                    out.push(format!("{} overlaps {}", self.areas[i].0.label(), self.areas[j].0.label()));
                }
            }
        }
        let covered: f64 = self.areas.iter().map(|(_, r)| r.area()).sum();
        if window.area() > 0.0 && (covered - window.area()).abs() / window.area() > 1e-9 {
            out.push(format!("the areas cover {covered:.1} of {:.1}", window.area()));
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn window() -> Rect {
        Rect { x: 0.0, y: 0.0, w: 1280.0, h: 720.0 }
    }

    #[test]
    fn the_default_layout_tiles_the_window_with_all_four_areas() {
        let l = Layout::default_for(window());
        assert!(l.check(window()).is_empty(), "{:?}", l.check(window()));
        for a in Area::ALL {
            assert!(l.rect(a).is_some(), "{} has no rectangle", a.label());
        }
    }

    #[test]
    fn the_viewport_gets_what_is_left() {
        let l = Layout::default_for(window());
        let v = l.rect(Area::Viewport).unwrap();
        let r = l.rect(Area::Recipe).unwrap();
        let i = l.rect(Area::Inspector).unwrap();
        assert!((r.w + v.w + i.w - 1280.0).abs() < 1e-9, "the top row spans the window");
        assert!((v.x - (r.x + r.w)).abs() < 1e-9, "and they sit next to each other");
    }

    #[test]
    fn a_missing_area_is_reported() {
        let mut l = Layout::default_for(window());
        l.areas.retain(|(a, _)| *a != Area::Inspector);
        let v = l.check(window());
        assert!(v.iter().any(|m| m.contains("Inspector is missing")), "{v:?}");
    }

    #[test]
    fn a_duplicated_area_is_reported() {
        let mut l = Layout::default_for(window());
        let vp = l.rect(Area::Viewport).unwrap();
        l.areas.push((Area::Viewport, vp));
        let v = l.check(window());
        assert!(v.iter().any(|m| m.contains("appears 2 times")), "{v:?}");
    }

    #[test]
    fn an_overlap_is_reported() {
        let mut l = Layout::default_for(window());
        l.areas[1].1.x -= 10.0;
        let v = l.check(window());
        assert!(v.iter().any(|m| m.contains("overlaps")), "{v:?}");
    }

    #[test]
    fn the_layout_scales_with_the_window() {
        let small = Rect { x: 0.0, y: 0.0, w: 800.0, h: 600.0 };
        let l = Layout::default_for(small);
        assert!(l.check(small).is_empty(), "{:?}", l.check(small));
    }
}
