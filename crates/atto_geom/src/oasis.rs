// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! OASIS container parsing (L-12), as far as it can go without zlib.
//!
//! OASIS (SEMI P39) is a record stream behind a magic string. Each record is a varint id, a varint
//! length and the payload. CBLOCK records may carry **zlib-compressed** bytes, and nearly every real
//! file uses them.
//!
//! This module parses the container: the magic string, the record framing, the version and the
//! table offsets. It **refuses a compressed CBLOCK with a message saying so** rather than returning
//! garbage — a decompressor is an external dependency this crate does not have, and inventing one
//! would be worse than not having the format.
//!
//! So L-12 is partial: a file can be opened and its structure read, and a real-world file stops at
//! the first compressed block with an honest reason.

/// The OASIS magic string, including its line endings.
pub const MAGIC: &[u8] = b"%SEMI-OASIS\r\n";

/// Why an OASIS file could not be read.
#[derive(Debug, Clone, PartialEq)]
pub enum OasisError {
    /// The file does not start with the magic string.
    NotOasis,
    /// A record ran past the end of the data.
    Truncated { at: usize },
    /// The file uses zlib compression, which this crate cannot decompress.
    Compressed { at: usize },
    /// A varint was longer than the format allows.
    BadVarint { at: usize },
}

impl std::fmt::Display for OasisError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            OasisError::NotOasis => write!(f, "oasis: the file does not start with %SEMI-OASIS"),
            OasisError::Truncated { at } => write!(f, "oasis: a record at byte {at} runs past the end"),
            OasisError::Compressed { at } => {
                write!(f, "oasis: a CBLOCK at byte {at} is zlib-compressed; this build has no decompressor")
            }
            OasisError::BadVarint { at } => write!(f, "oasis: a varint at byte {at} is malformed"),
        }
    }
}

/// One record, as framed in the file.
#[derive(Debug, Clone, PartialEq)]
pub struct Record {
    /// The record id.
    pub id: u64,
    /// Byte offset of the record's payload.
    pub offset: usize,
    /// The payload, uncompressed.
    pub payload: Vec<u8>,
}

/// The container: the records of a file that this build can read.
#[derive(Debug, Clone, PartialEq)]
pub struct OasisFile {
    /// Records in file order, up to the first one that could not be read.
    pub records: Vec<Record>,
    /// True when parsing stopped at a compressed block rather than at the end.
    pub stopped_at_compression: bool,
}

impl OasisFile {
    /// Whether the whole file was read.
    pub fn is_complete(&self) -> bool {
        !self.stopped_at_compression
    }
}

/// Read an unsigned varint (7 bits per byte, high bit continues).
fn varint(data: &[u8], at: usize) -> Result<(u64, usize), OasisError> {
    let mut value: u64 = 0;
    let mut shift = 0;
    let mut i = at;
    loop {
        if i >= data.len() {
            return Err(OasisError::Truncated { at });
        }
        let b = data[i];
        if shift >= 64 {
            return Err(OasisError::BadVarint { at });
        }
        value |= ((b & 0x7f) as u64) << shift;
        shift += 7;
        i += 1;
        if b & 0x80 == 0 {
            return Ok((value, i));
        }
    }
}

/// Parse the container of an OASIS file.
///
/// Stops at the first compressed CBLOCK and reports it: everything before is returned, so a caller
/// can say how far it got.
pub fn parse(data: &[u8]) -> Result<OasisFile, OasisError> {
    if data.len() < MAGIC.len() || &data[..MAGIC.len()] != MAGIC {
        return Err(OasisError::NotOasis);
    }
    let mut records = Vec::new();
    let mut i = MAGIC.len();
    let mut stopped = false;
    while i < data.len() {
        let start = i;
        let (id, after_id) = varint(data, i)?;
        let (len, after_len) = varint(data, after_id)?;
        let payload_start = after_len;
        let payload_end = payload_start + len as usize;
        if payload_end > data.len() {
            return Err(OasisError::Truncated { at: start });
        }
        // record id 34 is CBLOCK in P39; its payload carries a compression byte
        if id == 34 {
            if let Some(&flag) = data.get(payload_start) {
                if flag == 0 {
                    // uncompressed: the rest of the payload is the block
                    records.push(Record {
                        id,
                        offset: start,
                        payload: data[payload_start + 1..payload_end].to_vec(),
                    });
                    i = payload_end;
                    continue;
                }
                stopped = true;
                return Ok(OasisFile { records, stopped_at_compression: stopped });
            }
        }
        records.push(Record { id, offset: start, payload: data[payload_start..payload_end].to_vec() });
        i = payload_end;
    }
    Ok(OasisFile { records, stopped_at_compression: stopped })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enc(mut v: u64) -> Vec<u8> {
        let mut out = Vec::new();
        loop {
            let mut b = (v & 0x7f) as u8;
            v >>= 7;
            if v != 0 {
                b |= 0x80;
            }
            out.push(b);
            if v == 0 {
                return out;
            }
        }
    }

    fn file_with(records: &[(u64, &[u8])]) -> Vec<u8> {
        let mut out = MAGIC.to_vec();
        for (id, payload) in records {
            out.extend(enc(*id));
            out.extend(enc(payload.len() as u64));
            out.extend(*payload);
        }
        out
    }

    #[test]
    fn a_file_without_the_magic_string_is_refused() {
        assert_eq!(parse(b"not an oasis file"), Err(OasisError::NotOasis));
        assert_eq!(parse(b""), Err(OasisError::NotOasis));
    }

    #[test]
    fn the_records_of_an_uncompressed_file_are_read() {
        let data = file_with(&[(1, b"start"), (2, b"cell"), (3, b"end")]);
        let f = parse(&data).expect("a file");
        assert!(f.is_complete(), "nothing was skipped");
        assert_eq!(f.records.len(), 3, "three records");
        assert_eq!(f.records[0].id, 1);
        assert_eq!(f.records[0].payload, b"start".to_vec());
        assert_eq!(f.records[2].payload, b"end".to_vec());
    }

    #[test]
    fn an_uncompressed_cblock_is_read_through() {
        // compression flag 0 means the block is stored, not deflated
        let mut payload = vec![0u8];
        payload.extend(b"raw block contents");
        let data = file_with(&[(1, b"start"), (34, &payload), (3, b"end")]);
        let f = parse(&data).expect("a file");
        assert!(f.is_complete(), "an uncompressed block does not stop the read");
        assert_eq!(f.records.len(), 3);
        assert_eq!(f.records[1].payload, b"raw block contents".to_vec(), "the flag byte is stripped");
    }

    #[test]
    fn a_compressed_cblock_stops_the_read_and_says_why() {
        // the case every real file hits: flag 1 is zlib
        let mut payload = vec![1u8];
        payload.extend(b"deflated bytes");
        let data = file_with(&[(1, b"start"), (34, &payload), (3, b"end")]);
        let f = parse(&data).expect("a file");
        assert!(!f.is_complete(), "the read stopped");
        assert_eq!(f.records.len(), 1, "with what came before the block");
        assert!(f.stopped_at_compression, "and says where");
    }

    #[test]
    fn a_truncated_record_is_refused() {
        let mut data = MAGIC.to_vec();
        data.extend(enc(1));
        data.extend(enc(100)); // claims 100 bytes
        data.extend(b"short");
        assert!(matches!(parse(&data), Err(OasisError::Truncated { .. })));
    }

    #[test]
    fn a_multi_byte_varint_round_trips() {
        for v in [0u64, 1, 127, 128, 300, 16_384, 1_000_000] {
            let bytes = enc(v);
            let (back, _) = varint(&bytes, 0).expect("a varint");
            assert_eq!(back, v, "varint {v} round trips through {bytes:?}");
        }
    }
}
