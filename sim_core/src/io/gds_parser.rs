// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Binary GDSII stream-format parser and 2D polygon boolean cleaning.
//!
//! GDSII files are a sequence of fixed-format records: a 2-byte big-endian
//! length, a 1-byte record type, a 1-byte data type, followed by
//! `length - 4` bytes of payload. This module implements just enough of the
//! spec to walk STRUCTUREs and extract BOUNDARY (polygon) elements with
//! their LAYER and XY vertex data — sufficient to feed the lithography
//! fracturing stage. It intentionally does not implement the full GDSII
//! element set (PATH, SREF/AREF, TEXT) in this first pass.

use serde::{Deserialize, Serialize};
use thiserror::Error;

/// A single polygon extracted from a GDSII BOUNDARY element.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Polygon {
    /// GDS layer number.
    pub layer: i32,
    /// Vertices in database units (nm, after applying the file's UNITS
    /// record), in order, forming a closed ring (first point may or may not
    /// repeat as last point — callers should not assume either).
    pub vertices: Vec<(f64, f64)>,
}

/// A parsed GDSII library: flattened polygons from every structure, plus the
/// database-unit-to-meter scale factor.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GdsLibrary {
    pub polygons: Vec<Polygon>,
    /// Meters per database unit.
    pub db_unit_m: f64,
}

#[derive(Debug, Error)]
pub enum GdsError {
    #[error("unexpected end of file while reading record header")]
    UnexpectedEof,
    #[error("record length {0} is shorter than the 4-byte header")]
    RecordTooShort(u16),
    #[error("truncated record payload: expected {expected} bytes, found {found}")]
    TruncatedPayload { expected: usize, found: usize },
    #[error("{0}")]
    Format(String),
}

// GDSII record type codes relevant to boundary/layer/xy/units extraction.

/// Parse a GDSII byte stream into polygons in nanometres (hierarchy flattened, UNITS applied).
/// Thin wrapper over [`crate::io::gds`], kept for existing callers.
pub fn parse(bytes: &[u8]) -> Result<GdsLibrary, GdsError> {
    let lay = crate::io::gds::parse(bytes).map_err(|e| match e {
        crate::io::gds::GdsError::Eof => GdsError::UnexpectedEof,
        other => GdsError::Format(other.to_string()),
    })?;
    let flat = lay.flatten(None).map_err(|e| GdsError::Format(e.to_string()))?;
    Ok(GdsLibrary { polygons: flat.into_iter().map(|p| Polygon { layer: p.layer, vertices: p.pts }).collect(), db_unit_m: 1e-9 })
}

#[cfg(test)]
mod tests {
    const REC_BOUNDARY: u8 = 0x08;
    const REC_LAYER: u8 = 0x0D;
    const REC_XY: u8 = 0x10;
    const REC_ENDEL: u8 = 0x11;
    use super::*;

    fn write_record(out: &mut Vec<u8>, record_type: u8, data_type: u8, payload: &[u8]) {
        let length = (payload.len() + 4) as u16;
        out.extend_from_slice(&length.to_be_bytes());
        out.push(record_type);
        out.push(data_type);
        out.extend_from_slice(payload);
    }

    #[test]
    fn parses_a_single_square_boundary() {
        let mut bytes = Vec::new();
        // BOUNDARY (no data)
        write_record(&mut bytes, REC_BOUNDARY, 0, &[]);
        // LAYER = 1 (int16)
        write_record(&mut bytes, REC_LAYER, 2, &1i16.to_be_bytes());
        // XY: a closed square, four int32 pairs
        let mut xy_payload = Vec::new();
        for &(x, y) in &[(0, 0), (10, 0), (10, 10), (0, 10), (0, 0)] {
            xy_payload.extend_from_slice(&(x as i32).to_be_bytes());
            xy_payload.extend_from_slice(&(y as i32).to_be_bytes());
        }
        write_record(&mut bytes, REC_XY, 3, &xy_payload);
        // ENDEL
        write_record(&mut bytes, REC_ENDEL, 0, &[]);

        let lib = parse(&bytes).unwrap();
        assert_eq!(lib.polygons.len(), 1);
        assert_eq!(lib.polygons[0].layer, 1);
        assert_eq!(lib.polygons[0].vertices.len(), 4, "the repeated closing vertex is implicit in GDSII");
    }

    #[test]
    fn rejects_truncated_stream() {
        let bytes = vec![0x00]; // too short even for a header
        assert!(matches!(parse(&bytes), Err(GdsError::UnexpectedEof)));
    }
}
