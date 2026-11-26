// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Scaling from 80 % to 200 % without breaking the layout (PRD G-26).
//!
//! "Without layout breakage" has to mean something checkable. Here it means three
//! things at every scale in the range: no dimension goes negative, nothing that must
//! be visible collapses below a pixel, and the fixed rows plus the flexible row
//! still add up to the window.

/// The range the requirement names.
pub const MIN_SCALE: f64 = 0.8;
/// The top of the range.
pub const MAX_SCALE: f64 = 2.0;

/// A row of the layout: a fixed pixel height that scales, or the flexible remainder.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Row {
    /// A row with a design height in pixels at 100 %.
    Fixed(f64),
    /// The row that absorbs what is left.
    Flexible,
}

/// Heights in pixels at a given scale, for a window of `window_px`.
///
/// The flexible row takes the remainder, floored at zero: a window too small for the
/// fixed rows must produce an empty flexible row, never a negative one.
pub fn heights(rows: &[Row], window_px: f64, scale: f64) -> Vec<f64> {
    let fixed: f64 = rows
        .iter()
        .filter_map(|r| match r {
            Row::Fixed(h) => Some(h * scale),
            Row::Flexible => None,
        })
        .sum();
    let flex = (window_px - fixed).max(0.0);
    rows.iter()
        .map(|r| match r {
            Row::Fixed(h) => h * scale,
            Row::Flexible => flex,
        })
        .collect()
}

/// Check the layout at one scale, returning every violation.
pub fn check(rows: &[Row], window_px: f64, scale: f64) -> Vec<String> {
    let mut out = Vec::new();
    // a 5 % sweep accumulates error, so the ends are compared with a tolerance
    if !(scale >= MIN_SCALE - 1e-9 && scale <= MAX_SCALE + 1e-9) {
        out.push(format!("scale {scale} is outside the supported range"));
    }
    let hs = heights(rows, window_px, scale);
    for (i, h) in hs.iter().enumerate() {
        if !h.is_finite() {
            out.push(format!("row {i} is not a number at scale {scale}"));
        } else if *h < 0.0 {
            out.push(format!("row {i} is negative ({h:.1}) at scale {scale}"));
        } else if *h > 0.0 && *h < 1.0 {
            out.push(format!("row {i} collapses below a pixel ({h:.2}) at scale {scale}"));
        }
    }
    let total: f64 = hs.iter().sum();
    if total > window_px + 1e-9 {
        out.push(format!("the rows take {total:.1} px of a {window_px:.1} px window at {scale}"));
    }
    out
}

/// Check every scale in the supported range, in 5 % steps.
///
/// The whole range rather than its ends, because a layout can break in the middle —
/// at 80 % a fixed row can vanish, and at 200 % the fixed rows can exceed the window.
pub fn sweep(rows: &[Row], window_px: f64) -> Vec<String> {
    let mut out = Vec::new();
    let mut s = MIN_SCALE;
    while s <= MAX_SCALE + 1e-9 {
        out.extend(check(rows, window_px, s));
        s += 0.05;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn layout() -> Vec<Row> {
        vec![Row::Fixed(40.0), Row::Fixed(120.0), Row::Flexible, Row::Fixed(24.0)]
    }

    #[test]
    fn a_sane_layout_survives_the_whole_range() {
        let v = sweep(&layout(), 900.0);
        assert!(v.is_empty(), "{v:?}");
    }

    #[test]
    fn the_flexible_row_absorbs_the_change() {
        let a = heights(&layout(), 900.0, 1.0);
        let b = heights(&layout(), 900.0, 2.0);
        assert!((b[0] - a[0] * 2.0).abs() < 1e-9, "a fixed row doubles");
        assert!(b[2] < a[2], "and the flexible row gives up the difference");
        assert!((b.iter().sum::<f64>() - 900.0).abs() < 1e-9, "the total is unchanged");
    }

    #[test]
    fn a_window_too_small_gives_an_empty_row_not_a_negative_one() {
        // this is the breakage the requirement is about
        let hs = heights(&layout(), 100.0, 2.0);
        assert!(hs[2] >= 0.0, "the flexible row went negative: {}", hs[2]);
        assert!(check(&layout(), 100.0, 2.0).iter().any(|m| m.contains("take")));
    }

    #[test]
    fn a_row_that_vanishes_at_80_percent_is_reported() {
        // a 1 px design row becomes 0.8 px at 80 %, which is not drawable
        let rows = vec![Row::Fixed(1.0), Row::Flexible];
        let v = check(&rows, 900.0, 0.8);
        assert!(v.iter().any(|m| m.contains("below a pixel")), "{v:?}");
    }

    #[test]
    fn a_scale_outside_the_range_is_refused() {
        assert!(!check(&layout(), 900.0, 0.5).is_empty());
        assert!(!check(&layout(), 900.0, 3.0).is_empty());
        assert!(check(&layout(), 900.0, MIN_SCALE).is_empty());
        assert!(check(&layout(), 900.0, MAX_SCALE).is_empty());
    }
}
