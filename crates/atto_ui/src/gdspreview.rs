// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! What a dropped GDS file would show (PRD G-21).
//!
//! The drop itself is a gesture. What the file *is* — which layers it has, how many
//! polygons on each, and how big the drawing is — is a pure function of the parsed
//! layout, and it is the part that decides whether loading it is a good idea.
//!
//! The requirement adds a constraint worth enforcing here: the file loads **as the
//! mask only**, so the preview reports what would be taken and nothing else.

use atto_geom::gds::Layout;

/// One layer, summarised.
#[derive(Debug, Clone, PartialEq)]
pub struct LayerPreview {
    pub layer: i16,
    pub datatype: i16,
    pub polygons: usize,
    /// The layer's extent in microns, `(x0, y0, x1, y1)`.
    pub extent_um: (f64, f64, f64, f64),
}

/// The whole preview.
#[derive(Debug, Clone, PartialEq)]
pub struct Preview {
    pub layers: Vec<LayerPreview>,
    /// The file's extent across every layer, microns.
    pub extent_um: (f64, f64, f64, f64),
    /// Database units per micron, which decides the coordinate resolution.
    pub units_per_micron: f64,
}

impl Preview {
    /// The layers, in ascending order.
    pub fn layer_numbers(&self) -> Vec<i16> {
        self.layers.iter().map(|l| l.layer).collect()
    }

    /// Total polygons across every layer.
    pub fn polygons(&self) -> usize {
        self.layers.iter().map(|l| l.polygons).sum()
    }

    /// Whether the file is worth loading: something on a layer, and a unit.
    pub fn usable(&self) -> bool {
        self.polygons() > 0 && self.units_per_micron > 0.0
    }
}

/// Summarise a layout for the drop preview.
///
/// Returns `None` for a layout with no structure to draw — an empty file is a
/// failed drop, not a preview of nothing.
pub fn preview(layout: &Layout) -> Option<Preview> {
    if layout.structures.is_empty() || !(layout.units_per_micron > 0.0) {
        return None;
    }
    let mut layers: Vec<LayerPreview> = Vec::new();
    for s in &layout.structures {
        for p in &s.polygons {
            if p.points.is_empty() {
                continue;
            }
            let (mut x0, mut y0) = (f64::MAX, f64::MAX);
            let (mut x1, mut y1) = (f64::MIN, f64::MIN);
            for (x, y) in &p.points {
                x0 = x0.min(*x);
                y0 = y0.min(*y);
                x1 = x1.max(*x);
                y1 = y1.max(*y);
            }
            match layers.iter_mut().find(|l| l.layer == p.layer && l.datatype == p.datatype) {
                Some(l) => {
                    l.polygons += 1;
                    l.extent_um = (
                        l.extent_um.0.min(x0),
                        l.extent_um.1.min(y0),
                        l.extent_um.2.max(x1),
                        l.extent_um.3.max(y1),
                    );
                }
                None => layers.push(LayerPreview {
                    layer: p.layer,
                    datatype: p.datatype,
                    polygons: 1,
                    extent_um: (x0, y0, x1, y1),
                }),
            }
        }
    }
    if layers.is_empty() {
        return None;
    }
    layers.sort_by_key(|l| (l.layer, l.datatype));
    let mut extent = (f64::MAX, f64::MAX, f64::MIN, f64::MIN);
    for l in &layers {
        extent = (
            extent.0.min(l.extent_um.0),
            extent.1.min(l.extent_um.1),
            extent.2.max(l.extent_um.2),
            extent.3.max(l.extent_um.3),
        );
    }
    Some(Preview { layers, extent_um: extent, units_per_micron: layout.units_per_micron })
}

#[cfg(test)]
mod tests {
    use super::*;
    use atto_geom::gds::{Polygon, Structure};

    fn layout(layers: &[(i16, (f64, f64, f64, f64))]) -> Layout {
        Layout {
            units_per_micron: 1000.0,
            structures: vec![Structure {
                name: "TOP".into(),
                polygons: layers
                    .iter()
                    .map(|(l, (x0, y0, x1, y1))| Polygon {
                        layer: *l,
                        datatype: 0,
                        points: vec![(*x0, *y0), (*x1, *y0), (*x1, *y1), (*x0, *y1)],
                    })
                    .collect(),
            }],
        }
    }

    #[test]
    fn the_preview_lists_each_layer_once_in_order() {
        let p = preview(&layout(&[(5, (0.0, 0.0, 1.0, 1.0)), (1, (0.0, 0.0, 2.0, 2.0)), (5, (1.0, 1.0, 3.0, 3.0))])).unwrap();
        assert_eq!(p.layer_numbers(), vec![1, 5], "ascending, one entry per layer");
        assert_eq!(p.polygons(), 3, "two on layer 5, one on layer 1");
        let five = p.layers.iter().find(|l| l.layer == 5).unwrap();
        assert_eq!(five.polygons, 2);
        assert_eq!(five.extent_um, (0.0, 0.0, 3.0, 3.0), "the layer extent spans both");
    }

    #[test]
    fn the_file_extent_covers_every_layer() {
        let p = preview(&layout(&[(1, (-5.0, -5.0, 0.0, 0.0)), (2, (0.0, 0.0, 10.0, 2.0))])).unwrap();
        assert_eq!(p.extent_um, (-5.0, -5.0, 10.0, 2.0));
        assert!(p.usable());
    }

    #[test]
    fn an_empty_file_is_a_failed_drop_not_a_preview_of_nothing() {
        let empty = Layout { units_per_micron: 1000.0, structures: vec![] };
        assert!(preview(&empty).is_none());
        let no_polys = Layout {
            units_per_micron: 1000.0,
            structures: vec![Structure { name: "TOP".into(), polygons: vec![] }],
        };
        assert!(preview(&no_polys).is_none());
    }

    #[test]
    fn a_file_without_a_unit_cannot_be_loaded_as_a_mask() {
        let mut l = layout(&[(1, (0.0, 0.0, 1.0, 1.0))]);
        l.units_per_micron = 0.0;
        assert!(preview(&l).is_none(), "no unit means no real size");
    }
}
