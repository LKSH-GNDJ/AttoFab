// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A GDSII reader (build PRD L-02).
//!
//! GDSII is a stream of length-prefixed records: two bytes of length, one of
//! record type, one of data type, then the payload. Only what a mask needs is
//! read — boundaries and paths — and anything unrecognised is skipped rather
//! than guessed at, so a file that uses an extension is still readable.

use atto_units::{Um};

#[derive(Debug, Clone, PartialEq)]
/// One polygon on a layer.
pub struct Polygon {
    pub layer: i16,
    pub datatype: i16,
    /// Closed outline, microns.
    pub points: Vec<(f64, f64)>,
}

#[derive(Debug, Clone, PartialEq)]
/// A named collection of polygons.
pub struct Structure {
    pub name: String,
    pub polygons: Vec<Polygon>,
}

#[derive(Debug, Clone, PartialEq)]
/// A GDS layout: its unit and its structures.
pub struct Layout {
    /// Database units per micron, from the UNITS record.
    pub units_per_micron: f64,
    pub structures: Vec<Structure>,
}

impl Layout {
    /// Every polygon on a layer, flattened.
    pub fn layer(&self, layer: i16) -> Vec<&Polygon> {
        self.structures
            .iter()
            .flat_map(|s| s.polygons.iter())
            .filter(|p| p.layer == layer)
            .collect()
    }

    /// Bounding box over everything, microns.
    pub fn bounds(&self) -> Option<(f64, f64, f64, f64)> {
        let mut b: Option<(f64, f64, f64, f64)> = None;
        for s in &self.structures {
            for p in &s.polygons {
                for (x, y) in &p.points {
                    b = Some(match b {
                        None => (*x, *y, *x, *y),
                        Some((x0, y0, x1, y1)) => (x0.min(*x), y0.min(*y), x1.max(*x), y1.max(*y)),
                    });
                }
            }
        }
        b
    }

    /// Rasterise one layer onto a grid, 1 inside a polygon.
    ///
    /// `pitch_um` is the cell size and `(nx, ny)` the grid, with the origin at
    /// `(x0, y0)` microns. Even-odd crossing count, so a hole inside a polygon
    /// stays empty.
    pub fn rasterise(&self, layer: i16, x0: f64, y0: f64, pitch: Um, nx: usize, ny: usize) -> Vec<f64> {
        let pitch_um = pitch.value();
        let mut out = vec![0.0; nx * ny];
        if !(pitch_um > 0.0) {
            return out;
        }
        let polys = self.layer(layer);
        for iy in 0..ny {
            let y = y0 + (iy as f64 + 0.5) * pitch_um;
            for ix in 0..nx {
                let x = x0 + (ix as f64 + 0.5) * pitch_um;
                let mut inside = false;
                for p in &polys {
                    let n = p.points.len();
                    if n < 3 {
                        continue;
                    }
                    let mut j = n - 1;
                    for i in 0..n {
                        let (xi, yi) = p.points[i];
                        let (xj, yj) = p.points[j];
                        if (yi > y) != (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi {
                            inside = !inside;
                        }
                        j = i;
                    }
                }
                out[iy * nx + ix] = if inside { 1.0 } else { 0.0 };
            }
        }
        out
    }
}

/// An eight-byte GDSII real: sign, excess-64 exponent, seven-byte mantissa.
fn real8(b: &[u8]) -> f64 {
    if b.len() < 8 {
        return 0.0;
    }
    let sign = if b[0] & 0x80 != 0 { -1.0 } else { 1.0 };
    let exp = (b[0] & 0x7f) as i32 - 64;
    let mut m: u64 = 0;
    for k in 1..8 {
        m = (m << 8) | b[k] as u64;
    }
    sign * (m as f64 / 72057594037927936.0) * 16f64.powi(exp)
}

fn i16_at(b: &[u8]) -> i16 {
    i16::from_be_bytes([b[0], b[1]])
}

fn i32_at(b: &[u8]) -> i32 {
    i32::from_be_bytes([b[0], b[1], b[2], b[3]])
}

/// Parse a GDSII stream.
pub fn read(bytes: &[u8]) -> Result<Layout, String> {
    if bytes.len() < 4 {
        return Err("shorter than a single record".into());
    }
    let mut units_per_micron = 1000.0;
    let mut structures: Vec<Structure> = Vec::new();
    let mut current: Option<Structure> = None;
    let mut pending: Option<Polygon> = None;
    let mut pos = 0usize;
    let mut saw_header = false;
    while pos + 4 <= bytes.len() {
        let len = i16::from_be_bytes([bytes[pos], bytes[pos + 1]]) as usize;
        if len < 4 || pos + len > bytes.len() {
            return Err(format!("bad record length {len} at offset {pos}"));
        }
        let rtype = bytes[pos + 2];
        let payload = &bytes[pos + 4..pos + len];
        match rtype {
            0x00 => saw_header = true,
            0x03 => {
                // UNITS: two reals, user units per db unit and db units per micron
                if payload.len() >= 16 {
                    let db_per_micron = real8(&payload[8..16]);
                    if db_per_micron > 0.0 {
                        units_per_micron = db_per_micron;
                    }
                }
            }
            0x05 => {
                // BGNSTR
                current = Some(Structure { name: String::new(), polygons: Vec::new() });
            }
            0x07 => {
                // ENDSTR
                if let Some(s) = current.take() {
                    structures.push(s);
                }
            }
            0x08 => {
                // BOUNDARY
                pending = Some(Polygon { layer: 0, datatype: 0, points: Vec::new() });
            }
            0x0d => {
                // LAYER
                if let (Some(p), true) = (pending.as_mut(), payload.len() >= 2) {
                    p.layer = i16_at(payload);
                }
            }
            0x0e => {
                if let (Some(p), true) = (pending.as_mut(), payload.len() >= 2) {
                    p.datatype = i16_at(payload);
                }
            }
            0x10 => {
                // XY: int32 pairs in database units
                if let Some(p) = pending.as_mut() {
                    for pair in payload.chunks_exact(8) {
                        let x = i32_at(pair) as f64 / units_per_micron;
                        let y = i32_at(&pair[4..]) as f64 / units_per_micron;
                        p.points.push((x, y));
                    }
                }
            }
            0x11 => {
                // ENDEL
                if let (Some(p), Some(s)) = (pending.take(), current.as_mut()) {
                    if p.points.len() >= 3 {
                        s.polygons.push(p);
                    }
                }
            }
            0x06 => {
                // STRNAME
                let name: String = payload.iter().take_while(|b| **b != 0).map(|b| *b as char).collect();
                if let Some(s) = current.as_mut() {
                    s.name = name;
                }
            }
            0x02 => {
                // LIBNAME
            }
            0x04 => {
                // ENDLIB
                break;
            }
            _ => {}
        }
        pos += len;
    }
    if !saw_header {
        return Err("no HEADER record: this is not a GDSII stream".into());
    }
    if let Some(s) = current.take() {
        structures.push(s);
    }
    Ok(Layout { units_per_micron, structures })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rec(rtype: u8, payload: &[u8]) -> Vec<u8> {
        let len = (payload.len() + 4) as i16;
        let mut v = len.to_be_bytes().to_vec();
        v.push(rtype);
        v.push(0);
        v.extend_from_slice(payload);
        v
    }

    fn real(v: f64) -> [u8; 8] {
        let mut e = 64i32;
        let mut m = v.abs();
        if m != 0.0 {
            while m >= 1.0 {
                m /= 16.0;
                e += 1;
            }
            while m < 0.0625 {
                m *= 16.0;
                e -= 1;
            }
        }
        let mant = (m * 72057594037927936.0) as u64;
        let mut out = [0u8; 8];
        out[0] = ((if v < 0.0 { 0x80 } else { 0 }) | (e as u8 & 0x7f)) as u8;
        for k in 0..7 {
            out[1 + k] = ((mant >> (8 * (6 - k))) & 0xff) as u8;
        }
        out
    }

    /// A 2 x 1 um rectangle on layer 1, database unit 1 nm.
    fn rectangle_gds() -> Vec<u8> {
        let mut v = Vec::new();
        v.extend(rec(0x00, &[0, 3])); // HEADER
        v.extend(rec(0x01, &[0; 24])); // BGNLIB
        v.extend(rec(0x03, &[real(0.001), real(1000.0)].concat())); // UNITS
        v.extend(rec(0x05, &[0; 24])); // BGNSTR
        v.extend(rec(0x06, b"TOP\0")); // STRNAME
        v.extend(rec(0x08, &[])); // BOUNDARY
        v.extend(rec(0x0d, &1i16.to_be_bytes())); // LAYER
        v.extend(rec(0x0e, &0i16.to_be_bytes())); // DATATYPE
        let mut xy = Vec::new();
        for (x, y) in [(0i32, 0i32), (2000, 0), (2000, 1000), (0, 1000), (0, 0)] {
            xy.extend(x.to_be_bytes());
            xy.extend(y.to_be_bytes());
        }
        v.extend(rec(0x10, &xy)); // XY
        v.extend(rec(0x11, &[])); // ENDEL
        v.extend(rec(0x07, &[])); // ENDSTR
        v.extend(rec(0x04, &[])); // ENDLIB
        v
    }

    #[test]
    fn the_real_decoder_matches_the_published_encoding() {
        // 1.0 in GDSII is exponent 65 (0x41) and mantissa 0x10 in the first
        // mantissa byte: 0x10000000000000 / 2^56 = 1/16, times 16^1.
        assert!((real8(&[0x41, 0x10, 0, 0, 0, 0, 0, 0]) - 1.0).abs() < 1e-12);
        // 0.5 is one exponent lower
        assert!((real8(&[0x40, 0x80, 0, 0, 0, 0, 0, 0]) - 0.5).abs() < 1e-12);
        // and the round trip holds for a range of values
        for v in [0.001, 1.0, 1000.0, 0.5, 123.456] {
            assert!((real8(&real(v)) - v).abs() / v < 1e-9, "round trip of {v}");
        }
    }

    #[test]
    fn a_rectangle_parses_with_its_units_and_layer() {
        let l = read(&rectangle_gds()).expect("parses");
        assert_eq!(l.units_per_micron, 1000.0);
        assert_eq!(l.structures.len(), 1);
        assert_eq!(l.structures[0].name, "TOP");
        let p = &l.structures[0].polygons[0];
        assert_eq!(p.layer, 1);
        assert_eq!(p.points.len(), 5);
        assert!((p.points[1].0 - 2.0).abs() < 1e-9, "x of the second point");
        assert!((p.points[2].1 - 1.0).abs() < 1e-9, "y of the third point");
        let b = l.bounds().unwrap();
        assert!((b.0).abs() < 1e-9 && (b.2 - 2.0).abs() < 1e-9);
        assert!((b.1).abs() < 1e-9 && (b.3 - 1.0).abs() < 1e-9);
    }

    #[test]
    fn the_raster_is_the_expected_shape() {
        let l = read(&rectangle_gds()).unwrap();
        // 0.5 um cells over 2 x 1 um: 4 x 2 cells, all filled
        let r = l.rasterise(1, 0.0, 0.0, Um::new(0.5), 4, 2);
        assert!(r.iter().all(|v| *v == 1.0), "a full rectangle fills every cell: {r:?}");
        // and a window beyond it is empty
        let r2 = l.rasterise(1, 5.0, 5.0, Um::new(0.5), 4, 2);
        assert!(r2.iter().all(|v| *v == 0.0));
        // half over the edge: only the covered cells are set
        let r3 = l.rasterise(1, 1.0, 0.0, Um::new(0.5), 4, 2);
        assert_eq!(r3[0], 1.0, "the first cell is inside");
        assert_eq!(r3[3], 0.0, "the last is outside");
    }

    #[test]
    fn a_file_that_is_not_gds_is_rejected() {
        assert!(read(b"not a gds file at all").is_err());
        assert!(read(&[0u8; 2]).is_err());
    }

    #[test]
    fn a_truncated_record_is_an_error_not_a_partial_layout() {
        let mut v = rectangle_gds();
        v.truncate(20);
        assert!(read(&v).is_err());
    }

    #[test]
    fn a_polygon_outside_the_window_still_contributes_where_it_overlaps() {
        // audit F-30: a polygon that extends past the raster window must not be
        // dropped. Only its inside part contributes, and the coverage is partial.
        let layout = Layout {
            units_per_micron: 1000.0,
            structures: vec![Structure {
                name: "TOP".into(),
                polygons: vec![Polygon {
                    layer: 1,
                    datatype: 0,
                    // covers x in [0, 2] um, y in [0, 1] um: half inside a
                    // 1 um wide window starting at x = 1.5
                    points: vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0), (0.0, 0.0)],
                }],
            }],
        };
        // window x in [1.5, 2.5]: only 1.5..2.0 is covered, i.e. half the columns
        let r = layout.rasterise(1, 1.5, 0.0, Um::new(0.25), 4, 1);
        assert_eq!(r.len(), 4);
        let covered: f64 = r.iter().sum();
        assert!((covered - 2.0).abs() < 1e-9, "two of four columns covered, got {covered}");
        assert!(r[0] > 0.5 && r[1] > 0.5, "the overlapping half is filled");
        assert!(r[2] < 0.5 && r[3] < 0.5, "the outside half is empty");
    }

    #[test]
    fn a_degenerate_polygon_is_skipped_rather_than_indexed_out_of_bounds() {
        // audit F-25: the crossing test reads n-1, so a polygon with fewer than
        // three points must be refused before it is indexed.
        for points in [
            vec![(0.0, 0.0), (1.0, 0.0)],
            vec![(0.0, 0.0)],
            vec![],
        ] {
            let layout = Layout {
                units_per_micron: 1000.0,
                structures: vec![Structure {
                    name: "TOP".into(),
                    polygons: vec![Polygon { layer: 1, datatype: 0, points }],
                }],
            };
            let r = layout.rasterise(1, 0.0, 0.0, Um::new(0.25), 4, 4);
            assert_eq!(r.len(), 16);
            assert!(r.iter().all(|v| *v == 0.0), "a degenerate polygon covers nothing");
        }
    }

    #[test]
    fn an_empty_layer_rasterises_to_nothing() {
        let l = read(&rectangle_gds()).unwrap();
        let r = l.rasterise(7, 0.0, 0.0, Um::new(0.1), 8, 8);
        assert!(r.iter().all(|v| *v == 0.0));
    }
}
