// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Covering a whole layout with tiles, and turning many layers into one mask.
//!
//! # Why this exists
//!
//! The operator grid is capped at [`crate::MAX_GRID_AXIS`] cells per axis. At a 16 nm
//! cell that is 3,072 nm of process window, so a layout spanning 100 mm cannot be
//! simulated in one go at that resolution — it would need 6.25 million cells per axis.
//! Simulating a whole file therefore means **sweeping tiles** and reporting each one,
//! rather than enlarging a single simulation.
//!
//! This module does the two things that sweep needs and nothing else:
//!
//! * [`TileGrid::cover`] lays a grid of tiles over a bounding box;
//! * [`rasterise_union`] turns **every layer at once** into one boolean mask, which is
//!   what a combined etch sees.
//!
//! It works on plain polygons, so it does not depend on the GDS reader and can be tested
//! without one.

/// An axis-aligned box in nanometres.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bbox {
    pub min: (f64, f64),
    pub max: (f64, f64),
}

impl Bbox {
    /// A box from two corners, in any order.
    pub fn new(a: (f64, f64), b: (f64, f64)) -> Option<Self> {
        if !a.0.is_finite() || !a.1.is_finite() || !b.0.is_finite() || !b.1.is_finite() {
            return None;
        }
        let bbox = Bbox {
            min: (a.0.min(b.0), a.1.min(b.1)),
            max: (a.0.max(b.0), a.1.max(b.1)),
        };
        if bbox.width() <= 0.0 || bbox.height() <= 0.0 {
            return None;
        }
        Some(bbox)
    }

    /// Extent along x, nm.
    pub fn width(&self) -> f64 { self.max.0 - self.min.0 }
    /// Extent along y, nm.
    pub fn height(&self) -> f64 { self.max.1 - self.min.1 }
    /// The middle of the box.
    pub fn centre(&self) -> (f64, f64) {
        ((self.min.0 + self.max.0) * 0.5, (self.min.1 + self.max.1) * 0.5)
    }

    /// The box containing both.
    pub fn union(&self, other: &Bbox) -> Bbox {
        Bbox {
            min: (self.min.0.min(other.min.0), self.min.1.min(other.min.1)),
            max: (self.max.0.max(other.max.0), self.max.1.max(other.max.1)),
        }
    }
}

/// A closed polygon in nanometres. The last point joins the first.
#[derive(Debug, Clone, PartialEq)]
pub struct Polygon {
    pub layer: i32,
    pub datatype: i32,
    pub points: Vec<(f64, f64)>,
}

impl Polygon {
    /// The polygon's own bounding box, or `None` if it has fewer than three points.
    pub fn bbox(&self) -> Option<Bbox> {
        if self.points.len() < 3 {
            return None;
        }
        let mut min = (f64::INFINITY, f64::INFINITY);
        let mut max = (f64::NEG_INFINITY, f64::NEG_INFINITY);
        for p in &self.points {
            min = (min.0.min(p.0), min.1.min(p.1));
            max = (max.0.max(p.0), max.1.max(p.1));
        }
        Bbox::new(min, max)
    }

    /// Even-odd point-in-polygon test.
    pub fn contains(&self, x: f64, y: f64) -> bool {
        let n = self.points.len();
        if n < 3 {
            return false;
        }
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = self.points[i];
            let (xj, yj) = self.points[j];
            // the standard crossing test: count edges crossed by the ray to +x
            if (yi > y) != (yj > y) {
                let x_at = (xj - xi) * (y - yi) / (yj - yi) + xi;
                if x < x_at {
                    inside = !inside;
                }
            }
            j = i;
        }
        inside
    }
}

/// A grid of square tiles covering a bounding box.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TileGrid {
    /// Tile edge, nm.
    pub tile_nm: f64,
    /// Tiles along x and y.
    pub nx: usize,
    pub ny: usize,
    /// Lower-left corner of tile (0, 0).
    pub origin: (f64, f64),
}

impl TileGrid {
    /// Cover `bbox` with square tiles of edge `tile_nm`.
    ///
    /// Returns `None` for a tile edge that is not a positive finite number. A box smaller
    /// than one tile still gets one tile — a sweep of a small layout is one simulation,
    /// not zero.
    pub fn cover(bbox: Bbox, tile_nm: f64) -> Option<Self> {
        if !(tile_nm > 0.0) || !tile_nm.is_finite() {
            return None;
        }
        let nx = ((bbox.width() / tile_nm).ceil() as usize).max(1);
        let ny = ((bbox.height() / tile_nm).ceil() as usize).max(1);
        Some(TileGrid { tile_nm, nx, ny, origin: bbox.min })
    }

    /// How many tiles the sweep runs.
    pub fn count(&self) -> usize { self.nx * self.ny }

    /// The centre of tile `(ix, iy)`.
    pub fn centre(&self, ix: usize, iy: usize) -> (f64, f64) {
        (
            self.origin.0 + (ix as f64 + 0.5) * self.tile_nm,
            self.origin.1 + (iy as f64 + 0.5) * self.tile_nm,
        )
    }

    /// The box a tile covers.
    pub fn tile_bbox(&self, ix: usize, iy: usize) -> Bbox {
        let c = self.centre(ix, iy);
        let h = self.tile_nm * 0.5;
        Bbox { min: (c.0 - h, c.1 - h), max: (c.0 + h, c.1 + h) }
    }

    /// Every tile as `(ix, iy, centre)`.
    pub fn tiles(&self) -> impl Iterator<Item = (usize, usize, (f64, f64))> + '_ {
        (0..self.ny).flat_map(move |iy| (0..self.nx).map(move |ix| (ix, iy, self.centre(ix, iy))))
    }
}

/// Rasterise the **union of every polygon** over `region` into a `cells x cells` mask.
///
/// `true` means at least one polygon covers the cell centre — that is the mask a combined
/// etch sees when all layers are used together. Layers are not distinguished here on
/// purpose: the caller decides which polygons to pass.
///
/// Returns `None` for a cell count of zero, a non-finite region, or a region with no area.
pub fn rasterise_union(polys: &[Polygon], region: Bbox, cells: usize) -> Option<Vec<bool>> {
    if cells == 0 || !region.width().is_finite() || !region.height().is_finite() {
        return None;
    }
    let dx = region.width() / cells as f64;
    let dy = region.height() / cells as f64;
    let mut mask = vec![false; cells * cells];
    for iy in 0..cells {
        let y = region.min.1 + (iy as f64 + 0.5) * dy;
        for ix in 0..cells {
            let x = region.min.0 + (ix as f64 + 0.5) * dx;
            mask[iy * cells + ix] = polys.iter().any(|p| p.contains(x, y));
        }
    }
    Some(mask)
}

/// The bounding box of every polygon given, or `None` if none has area.
pub fn extent(polys: &[Polygon]) -> Option<Bbox> {
    let mut acc: Option<Bbox> = None;
    for p in polys {
        if let Some(b) = p.bbox() {
            acc = Some(match acc { Some(a) => a.union(&b), None => b });
        }
    }
    acc
}

/// Which GDS layers are present, as `(layer, datatype)`, sorted.
pub fn layers(polys: &[Polygon]) -> Vec<(i32, i32)> {
    let mut v: Vec<(i32, i32)> = polys.iter().map(|p| (p.layer, p.datatype)).collect();
    v.sort_unstable();
    v.dedup();
    v
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(layer: i32, x0: f64, y0: f64, side: f64) -> Polygon {
        Polygon {
            layer,
            datatype: 0,
            points: vec![(x0, y0), (x0 + side, y0), (x0 + side, y0 + side), (x0, y0 + side)],
        }
    }

    #[test]
    fn a_tile_grid_covers_the_whole_extent() {
        let b = Bbox::new((0.0, 0.0), (10_000.0, 6_000.0)).unwrap();
        let g = TileGrid::cover(b, 2_000.0).unwrap();
        assert_eq!((g.nx, g.ny), (5, 3));
        assert_eq!(g.count(), 15);
        // the last tile still reaches the far corner
        let last = g.tile_bbox(g.nx - 1, g.ny - 1);
        assert!(last.max.0 >= 10_000.0 && last.max.1 >= 6_000.0);
    }

    #[test]
    fn a_layout_smaller_than_one_tile_is_one_tile() {
        let b = Bbox::new((0.0, 0.0), (100.0, 100.0)).unwrap();
        let g = TileGrid::cover(b, 2_000.0).unwrap();
        assert_eq!((g.nx, g.ny), (1, 1));
        assert_eq!(g.count(), 1);
    }

    #[test]
    fn tile_centres_are_evenly_spaced_and_inside_their_tile() {
        let b = Bbox::new((0.0, 0.0), (4_000.0, 4_000.0)).unwrap();
        let g = TileGrid::cover(b, 1_000.0).unwrap();
        let mut seen = 0;
        for (ix, iy, c) in g.tiles() {
            let tb = g.tile_bbox(ix, iy);
            assert!(c.0 > tb.min.0 && c.0 < tb.max.0);
            assert!(c.1 > tb.min.1 && c.1 < tb.max.1);
            seen += 1;
        }
        assert_eq!(seen, g.count());
    }

    #[test]
    fn a_tile_edge_with_no_meaning_is_refused() {
        let b = Bbox::new((0.0, 0.0), (10.0, 10.0)).unwrap();
        assert!(TileGrid::cover(b, 0.0).is_none());
        assert!(TileGrid::cover(b, -1.0).is_none());
        assert!(TileGrid::cover(b, f64::NAN).is_none());
        assert!(Bbox::new((0.0, 0.0), (0.0, 10.0)).is_none(), "no width");
    }

    #[test]
    fn the_union_of_two_layers_keeps_both() {
        // the point of the union mask: a feature on either layer is still a feature
        let polys = vec![square(2, 0.0, 0.0, 100.0), square(95, 200.0, 0.0, 100.0)];
        let region = Bbox::new((0.0, 0.0), (300.0, 100.0)).unwrap();
        let mask = rasterise_union(&polys, region, 30).unwrap();
        let covered = mask.iter().filter(|b| **b).count();
        // two 100x100 squares side by side in a 300x100 region cover two thirds of it:
        // 30 columns by 30 rows, ten columns of each square, so about 600 of 900
        assert!(covered > 560 && covered < 640, "covered {covered} of 900");
        // a cell in the gap between them is empty
        assert!(!mask[5 * 30 + 15], "the gap must not be covered");
    }

    #[test]
    fn overlapping_polygons_do_not_double_count() {
        let polys = vec![square(1, 0.0, 0.0, 100.0), square(2, 50.0, 0.0, 100.0)];
        let region = Bbox::new((0.0, 0.0), (150.0, 100.0)).unwrap();
        let mask = rasterise_union(&polys, region, 15).unwrap();
        assert!(mask.iter().all(|b| *b), "the two overlap into full coverage");
    }

    #[test]
    fn the_extent_and_layers_describe_the_layout() {
        let polys = vec![square(2, 0.0, 0.0, 100.0), square(95, 500.0, 500.0, 100.0)];
        let e = extent(&polys).unwrap();
        assert_eq!(e.min, (0.0, 0.0));
        assert_eq!(e.max, (600.0, 600.0));
        assert_eq!(layers(&polys), vec![(2, 0), (95, 0)]);
    }

    #[test]
    fn a_point_test_respects_the_polygon_edges() {
        let p = square(1, 0.0, 0.0, 100.0);
        assert!(p.contains(50.0, 50.0));
        assert!(!p.contains(-1.0, 50.0));
        assert!(!p.contains(101.0, 50.0));
        assert!(!p.contains(50.0, 101.0));
    }

    #[test]
    fn an_empty_polygon_set_has_no_extent() {
        assert!(extent(&[]).is_none());
        let degenerate = vec![Polygon { layer: 1, datatype: 0, points: vec![(0.0, 0.0), (1.0, 1.0)] }];
        assert!(extent(&degenerate).is_none(), "two points are not a polygon");
    }
}
