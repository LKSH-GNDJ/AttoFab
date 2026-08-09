//! GDSII / OASIS Layout Parser (Phase 3).
//!
//! A minimal, self-contained GDSII stream reader - deliberately hand-rolled
//! rather than depending on an external crate (e.g. gds21), since the
//! GDSII binary format is a small, fully-specified, deterministic record
//! format (not a statistical fitting problem like Pearson-IV), so a
//! direct implementation is both tractable and independently auditable.
//! OASIS is not implemented (GDSII is the format actually used by the
//! open-source flows AttoFab targets, e.g. OpenLane output).
//!
//! GDSII stream structure: a sequence of records, each:
//!   [2 bytes: record length, big-endian, includes this 4-byte header]
//!   [1 byte:  record type]
//!   [1 byte:  data type]
//!   [(length-4) bytes: data, interpreted per data type]
//!
//! This parser extracts BOUNDARY elements (polygon outlines) per layer,
//! which is what's needed to drive a lithography mask - it does not
//! attempt to resolve SREF/AREF (cell instance references) or PATH
//! elements, which are out of scope for the mask-rasterization use case
//! this module serves.

use std::fmt;

// GDSII record type codes (high byte of the 2-byte record type field, per
// the standard GDSII stream spec).
mod rtype {
    pub const HEADER: u8 = 0x00;
    #[allow(dead_code)]
    pub const BGNLIB: u8 = 0x01;
    #[allow(dead_code)]
    pub const LIBNAME: u8 = 0x02;
    pub const UNITS: u8 = 0x03;
    pub const ENDLIB: u8 = 0x04;
    pub const BGNSTR: u8 = 0x05;
    pub const STRNAME: u8 = 0x06;
    pub const ENDSTR: u8 = 0x07;
    pub const BOUNDARY: u8 = 0x08;
    pub const LAYER: u8 = 0x0d;
    pub const XY: u8 = 0x10;
    pub const ENDEL: u8 = 0x11;
}

#[derive(Debug)]
pub enum GdsError {
    UnexpectedEof,
    InvalidRecordLength,
    NotAGdsFile,
}

impl fmt::Display for GdsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GdsError::UnexpectedEof => write!(f, "unexpected end of GDSII stream"),
            GdsError::InvalidRecordLength => write!(f, "invalid GDSII record length"),
            GdsError::NotAGdsFile => write!(f, "file does not start with a GDSII HEADER record"),
        }
    }
}

impl std::error::Error for GdsError {}

/// Decode an 8-byte GDSII "real" number: 1 sign bit, 7-bit excess-64
/// base-16 exponent, 56-bit mantissa. This is NOT IEEE 754.
///
/// value = sign * (mantissa / 2^56) * 16^(exponent - 64)
fn decode_real8(bytes: &[u8; 8]) -> f64 {
    let sign = if bytes[0] & 0x80 != 0 { -1.0 } else { 1.0 };
    let exponent = (bytes[0] & 0x7f) as i32 - 64;
    let mut mantissa: u64 = 0;
    for &b in &bytes[1..8] {
        mantissa = (mantissa << 8) | b as u64;
    }
    let mantissa_frac = mantissa as f64 / (1u64 << 56) as f64;
    sign * mantissa_frac * 16f64.powi(exponent)
}

/// Encode a positive f64 into GDSII 8-byte real format. Used by the test
/// suite to build synthetic GDSII files with known values, and available
/// for anything that needs to *write* GDSII in the future.
pub fn encode_real8(value: f64) -> [u8; 8] {
    if value == 0.0 {
        return [0u8; 8];
    }
    let sign = value < 0.0;
    let mut v = value.abs();
    let mut exponent = 64i32;

    // Normalize so 1/16 <= v < 1 (i.e. the mantissa's top nibble is
    // nonzero once scaled to 56 bits), adjusting the base-16 exponent.
    while v >= 1.0 {
        v /= 16.0;
        exponent += 1;
    }
    while v < 1.0 / 16.0 {
        v *= 16.0;
        exponent -= 1;
    }

    let mantissa = (v * (1u64 << 56) as f64).round() as u64;
    let mut bytes = [0u8; 8];
    bytes[0] = (exponent as u8 & 0x7f) | if sign { 0x80 } else { 0 };
    for i in 0..7 {
        bytes[7 - i] = ((mantissa >> (i * 8)) & 0xff) as u8;
    }
    bytes
}

struct RawRecord {
    rtype: u8,
    #[allow(dead_code)]
    dtype: u8,
    data: Vec<u8>,
}

fn read_records(bytes: &[u8]) -> Result<Vec<RawRecord>, GdsError> {
    let mut records = Vec::new();
    let mut i = 0usize;
    while i < bytes.len() {
        if i + 4 > bytes.len() {
            return Err(GdsError::UnexpectedEof);
        }
        let len = u16::from_be_bytes([bytes[i], bytes[i + 1]]) as usize;
        if len < 4 {
            return Err(GdsError::InvalidRecordLength);
        }
        if i + len > bytes.len() {
            return Err(GdsError::UnexpectedEof);
        }
        let rtype = bytes[i + 2];
        let dtype = bytes[i + 3];
        let data = bytes[i + 4..i + len].to_vec();
        records.push(RawRecord { rtype, dtype, data });
        i += len;
    }
    Ok(records)
}

fn parse_i32_array(data: &[u8]) -> Vec<i32> {
    data.chunks_exact(4).map(|c| i32::from_be_bytes([c[0], c[1], c[2], c[3]])).collect()
}

#[derive(Debug, Clone)]
pub struct GdsPolygon {
    pub layer: i16,
    /// Polygon vertices in database units (integer GDSII coordinates).
    pub points_dbu: Vec<(i32, i32)>,
}

#[derive(Debug, Clone)]
pub struct GdsStructure {
    pub name: String,
    pub polygons: Vec<GdsPolygon>,
}

#[derive(Debug, Clone)]
pub struct GdsLibrary {
    /// Micrometers per database unit - use this to convert `points_dbu`
    /// into physical coordinates for mask rasterization.
    pub um_per_dbu: f64,
    pub structures: Vec<GdsStructure>,
}

/// Parse a GDSII byte stream, extracting BOUNDARY polygons per structure.
/// SREF/AREF (cell references) and PATH elements are not resolved.
pub fn parse_gds(bytes: &[u8]) -> Result<GdsLibrary, GdsError> {
    let records = read_records(bytes)?;
    if records.is_empty() || records[0].rtype != rtype::HEADER {
        return Err(GdsError::NotAGdsFile);
    }

    let mut um_per_dbu = 0.001; // GDSII default-ish fallback (1000 dbu/um) if UNITS is missing
    let mut structures = Vec::new();

    let mut cur_struct_name: Option<String> = None;
    let mut cur_polygons: Vec<GdsPolygon> = Vec::new();

    let mut in_boundary = false;
    let mut cur_layer: i16 = 0;
    let mut cur_points: Vec<(i32, i32)> = Vec::new();

    for rec in &records {
        match rec.rtype {
            t if t == rtype::UNITS => {
                // UNITS data is two REAL8 values: [user_units_per_dbu, meters_per_dbu].
                if rec.data.len() >= 16 {
                    let meters_per_dbu_bytes: [u8; 8] = rec.data[8..16].try_into().unwrap();
                    let meters_per_dbu = decode_real8(&meters_per_dbu_bytes);
                    um_per_dbu = meters_per_dbu * 1e6;
                }
            }
            t if t == rtype::BGNSTR => {
                cur_struct_name = None;
                cur_polygons = Vec::new();
            }
            t if t == rtype::STRNAME => {
                let name = String::from_utf8_lossy(&rec.data).trim_end_matches('\0').to_string();
                cur_struct_name = Some(name);
            }
            t if t == rtype::BOUNDARY => {
                in_boundary = true;
                cur_layer = 0;
                cur_points.clear();
            }
            t if t == rtype::LAYER && in_boundary => {
                if rec.data.len() >= 2 {
                    cur_layer = i16::from_be_bytes([rec.data[0], rec.data[1]]);
                }
            }
            t if t == rtype::XY && in_boundary => {
                let coords = parse_i32_array(&rec.data);
                cur_points = coords.chunks_exact(2).map(|c| (c[0], c[1])).collect();
            }
            t if t == rtype::ENDEL && in_boundary => {
                cur_polygons.push(GdsPolygon { layer: cur_layer, points_dbu: cur_points.clone() });
                in_boundary = false;
            }
            t if t == rtype::ENDSTR => {
                if let Some(name) = cur_struct_name.take() {
                    structures.push(GdsStructure { name, polygons: std::mem::take(&mut cur_polygons) });
                }
            }
            t if t == rtype::ENDLIB => break,
            _ => {}
        }
    }

    Ok(GdsLibrary { um_per_dbu, structures })
}

/// Rasterize a single polygon (vertices in um, already unit-converted)
/// onto a 2D boolean mask (length nx*ny, row-major y*nx+x) using the
/// standard even-odd point-in-polygon test, evaluated at each grid cell's
/// center - correctly handles arbitrary simple polygons, not just
/// axis-aligned rectangles (standard-cell layouts are mostly rectangles,
/// but this doesn't assume that).
pub fn rasterize_polygon_um(
    points_um: &[(f64, f64)],
    nx: usize,
    ny: usize,
    dx_um: f64,
    dy_um: f64,
    mask: &mut [bool],
) {
    if points_um.len() < 3 {
        return;
    }
    let n = points_um.len();
    for yi in 0..ny {
        let py = (yi as f64 + 0.5) * dy_um;
        for xi in 0..nx {
            let px = (xi as f64 + 0.5) * dx_um;
            let mut inside = false;
            let mut j = n - 1;
            for i in 0..n {
                let (xi_pt, yi_pt) = points_um[i];
                let (xj_pt, yj_pt) = points_um[j];
                if (yi_pt > py) != (yj_pt > py) {
                    let x_cross = (xj_pt - xi_pt) * (py - yi_pt) / (yj_pt - yi_pt) + xi_pt;
                    if px < x_cross {
                        inside = !inside;
                    }
                }
                j = i;
            }
            if inside {
                mask[yi * nx + xi] = true;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(rtype: u8, dtype: u8, data: &[u8]) -> Vec<u8> {
        let len = 4 + data.len();
        let mut out = Vec::with_capacity(len);
        out.extend_from_slice(&(len as u16).to_be_bytes());
        out.push(rtype);
        out.push(dtype);
        out.extend_from_slice(data);
        out
    }

    /// Hand-assembles a minimal, valid GDSII stream in memory (no external
    /// file needed): HEADER, BGNLIB, LIBNAME, UNITS, BGNSTR, STRNAME,
    /// BOUNDARY(layer=1, a 2x1 um rectangle), ENDEL, ENDSTR, ENDLIB.
    fn synthetic_gds_with_rectangle() -> Vec<u8> {
        let mut buf = Vec::new();
        buf.extend(record(rtype::HEADER, 0x02, &600i16.to_be_bytes()));
        buf.extend(record(rtype::BGNLIB, 0x02, &[0u8; 24])); // 12x int16 timestamp fields, zeroed
        buf.extend(record(rtype::LIBNAME, 0x06, b"TESTLIB\0"));

        // UNITS: user-units-per-dbu, meters-per-dbu. 1 dbu = 1nm = 1e-9 m.
        let mut units_data = Vec::new();
        units_data.extend_from_slice(&encode_real8(0.001)); // 1000 dbu per user-unit
        units_data.extend_from_slice(&encode_real8(1e-9)); // meters per dbu
        buf.extend(record(rtype::UNITS, 0x05, &units_data));

        buf.extend(record(rtype::BGNSTR, 0x02, &[0u8; 24]));
        buf.extend(record(rtype::STRNAME, 0x06, b"TOP\0"));

        buf.extend(record(rtype::BOUNDARY, 0x00, &[]));
        buf.extend(record(rtype::LAYER, 0x02, &1i16.to_be_bytes()));
        // Rectangle from (0,0) to (2000, 1000) dbu = (2um, 1um) at 1nm/dbu,
        // closed polygon (first point repeated at the end, per GDSII spec).
        let xy: Vec<i32> = vec![0, 0, 2000, 0, 2000, 1000, 0, 1000, 0, 0];
        let mut xy_data = Vec::new();
        for v in xy {
            xy_data.extend_from_slice(&v.to_be_bytes());
        }
        buf.extend(record(rtype::XY, 0x03, &xy_data));
        buf.extend(record(rtype::ENDEL, 0x00, &[]));

        buf.extend(record(rtype::ENDSTR, 0x00, &[]));
        buf.extend(record(rtype::ENDLIB, 0x00, &[]));
        buf
    }

    #[test]
    fn real8_decodes_known_literal_byte_patterns() {
        // 1.0 in GDSII real8: sign=0, exponent=65 (0x41, excess-64 => 1),
        // mantissa top byte 0x10 (=16 decimal => 16/256 top-byte-fraction,
        // i.e. mantissa/2^56 = 1/16), giving 1/16 * 16^1 = 1.0 exactly.
        // This is the standard textbook example for the format, not just
        // a self-consistency check against our own encoder.
        let one = [0x41, 0x10, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        assert!((decode_real8(&one) - 1.0).abs() < 1e-12);

        // -2.0: same exponent, mantissa doubled (top byte 0x20), sign bit set.
        let neg_two = [0xC1, 0x20, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00];
        assert!((decode_real8(&neg_two) - (-2.0)).abs() < 1e-12);
    }

    #[test]
    fn encode_decode_real8_round_trips() {
        for v in [1.0, 0.001, 1e-9, 123.456, 0.0625, 16.0] {
            let encoded = encode_real8(v);
            let decoded = decode_real8(&encoded);
            let rel_err = (decoded - v).abs() / v.abs();
            assert!(rel_err < 1e-9, "round-trip failed for {v}: got {decoded}");
        }
    }

    #[test]
    fn parses_synthetic_gds_and_extracts_rectangle() {
        let bytes = synthetic_gds_with_rectangle();
        let lib = parse_gds(&bytes).expect("should parse");

        assert!((lib.um_per_dbu - 0.001).abs() < 1e-12, "expected 1nm/dbu, got {}", lib.um_per_dbu);
        assert_eq!(lib.structures.len(), 1);
        let s = &lib.structures[0];
        assert_eq!(s.name, "TOP");
        assert_eq!(s.polygons.len(), 1);
        assert_eq!(s.polygons[0].layer, 1);
        assert_eq!(s.polygons[0].points_dbu, vec![(0, 0), (2000, 0), (2000, 1000), (0, 1000), (0, 0)]);
    }

    #[test]
    fn rejects_non_gds_data() {
        let garbage = b"not a gds file at all, just some text".to_vec();
        assert!(parse_gds(&garbage).is_err());
    }

    #[test]
    fn rasterize_rectangle_produces_correct_mask() {
        // 2um x 1um rectangle on a 10x4 grid at 1um/cell - should mark
        // (x=0,y=0) and (x=1,y=0) as inside (cell centers 0.5,0.5 and
        // 1.5,0.5 both fall within [0,2]x[0,1]), and nothing else.
        let points_um = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 1.0), (0.0, 1.0), (0.0, 0.0)];
        let nx = 10;
        let ny = 4;
        let mut mask = vec![false; nx * ny];
        rasterize_polygon_um(&points_um, nx, ny, 1.0, 1.0, &mut mask);

        let expected_true: Vec<usize> = vec![0, 1]; // (x=0,y=0) and (x=1,y=0), row-major indices 0 and 1
        for i in 0..mask.len() {
            let expected = expected_true.contains(&i);
            assert_eq!(mask[i], expected, "mask[{i}] (x={}, y={}) mismatch", i % nx, i / nx);
        }
    }

    #[test]
    fn rasterize_rectangle_covers_multiple_rows() {
        // Same rectangle but taller (2um in y) so it should cover two full
        // rows (y=0 and y=1) at x=0,1.
        let points_um = vec![(0.0, 0.0), (2.0, 0.0), (2.0, 2.0), (0.0, 2.0), (0.0, 0.0)];
        let nx = 5;
        let ny = 5;
        let mut mask = vec![false; nx * ny];
        rasterize_polygon_um(&points_um, nx, ny, 1.0, 1.0, &mut mask);

        for y in 0..5 {
            for x in 0..5 {
                let expected = (x == 0 || x == 1) && (y == 0 || y == 1);
                assert_eq!(mask[y * nx + x], expected, "({x},{y}) mismatch");
            }
        }
    }
}
