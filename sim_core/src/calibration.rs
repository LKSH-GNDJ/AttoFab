// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Calibration-data audit (provenance, not physics).
//!
//! The bundled `data/calibration` pack is *unverified*: its manifest records file
//! hashes (identity) but no source, conditions, uncertainty or licence. This module
//! reads that manifest, re-hashes the files, cross-checks the implant table against the
//! built-in range table, and classifies the pack. A stage may only be called
//! `Calibrated` when the pack is classified [`Evidence::Experimental`]; no bundled pack is.

use crate::params::Species;
use serde::{Deserialize, Serialize};
use std::path::Path;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Evidence {
    /// Hash/identity only; source and conditions undocumented.
    Unverified,
    /// Published reference values with documented source, not fitted to a tool.
    Reference,
    /// Measured on a tool with documented conditions and uncertainty.
    Experimental,
}

#[derive(Debug, Clone, Serialize)]
pub struct FileAudit { pub path: String, pub declared_status: String, pub hash_ok: Option<bool>, pub has_source: bool, pub has_conditions: bool, pub has_uncertainty: bool }

#[derive(Debug, Clone, Serialize)]
pub struct TableDeviation { pub species: String, pub energy_kev: f64, pub bundled_rp_nm: f64, pub builtin_rp_nm: f64, pub deviation_pct: f64 }

#[derive(Debug, Clone, Serialize)]
pub struct CalibrationStatus {
    pub evidence: Evidence,
    pub files: Vec<FileAudit>,
    pub implant_deviations: Vec<TableDeviation>,
    pub notes: Vec<String>,
}

impl CalibrationStatus {
    pub fn max_implant_deviation_pct(&self) -> f64 { self.implant_deviations.iter().map(|d| d.deviation_pct.abs()).fold(0.0, f64::max) }
    pub fn allows_calibrated_label(&self) -> bool { self.evidence == Evidence::Experimental }
}

fn sha256_hex(bytes: &[u8]) -> String {
    // Minimal SHA-256 (no extra dependency); verified against a known vector in the tests.
    const K: [u32; 64] = [0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2];
    let mut h: [u32; 8] = [0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19];
    let mut msg = bytes.to_vec();
    let bitlen = (bytes.len() as u64) * 8;
    msg.push(0x80);
    while msg.len() % 64 != 56 { msg.push(0); }
    msg.extend_from_slice(&bitlen.to_be_bytes());
    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 { w[i] = u32::from_be_bytes([chunk[4*i], chunk[4*i+1], chunk[4*i+2], chunk[4*i+3]]); }
        for i in 16..64 {
            let s0 = w[i-15].rotate_right(7) ^ w[i-15].rotate_right(18) ^ (w[i-15] >> 3);
            let s1 = w[i-2].rotate_right(17) ^ w[i-2].rotate_right(19) ^ (w[i-2] >> 10);
            w[i] = w[i-16].wrapping_add(s0).wrapping_add(w[i-7]).wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
        }
        for i in 0..8 { h[i] = h[i].wrapping_add(v[i]); }
    }
    h.iter().map(|x| format!("{x:08x}")).collect()
}

/// Audit a calibration directory containing `manifest.json` (paths in the manifest are relative to `root`).
pub fn audit(root: &Path, dir: &Path) -> Result<CalibrationStatus, String> {
    let manifest_path = dir.join("manifest.json");
    let text = std::fs::read_to_string(&manifest_path).map_err(|e| format!("{}: {e}", manifest_path.display()))?;
    let v: serde_json::Value = serde_json::from_str(&text).map_err(|e| e.to_string())?;
    let mut files = Vec::new();
    let mut notes = Vec::new();
    for f in v["files"].as_array().ok_or("manifest has no `files` array")? {
        let path = f["path"].as_str().unwrap_or("").to_string();
        let hash_ok = std::fs::read(root.join(&path)).ok().map(|b| Some(sha256_hex(&b)) == f["sha256"].as_str().map(|s| s.to_lowercase()));
        files.push(FileAudit {
            declared_status: f["evidence_status"].as_str().unwrap_or("unknown").to_string(), hash_ok,
            has_source: !f["source"].is_null(), has_conditions: !f["conditions"].is_null(), has_uncertainty: !f["uncertainty"].is_null(), path,
        });
    }
    // Cross-check the bundled implant range table against the built-in one.
    let mut dev = Vec::new();
    let imp = dir.join("implant");
    if let Ok(rd) = std::fs::read_dir(&imp) {
        for e in rd.flatten() {
            let Ok(t) = std::fs::read_to_string(e.path()) else { continue };
            let Ok(j) = serde_json::from_str::<serde_json::Value>(&t) else { continue };
            let sp = match j["species"].as_str() { Some("Boron") => Species::Boron, Some("Phosphorus") => Species::Phosphorus, Some("Arsenic") => Species::Arsenic, _ => continue };
            let (es, rs) = (j["energies_keV"].as_array(), j["range_nm"].as_array());
            if let (Some(es), Some(rs)) = (es, rs) {
                for (e, r) in es.iter().zip(rs) {
                    let (e, r) = (e.as_f64().unwrap_or(0.0), r.as_f64().unwrap_or(0.0));
                    if e <= 0.0 || r <= 0.0 { continue; }
                    let b = crate::pipeline::implant_range_nm(sp, e);
                    dev.push(TableDeviation { species: sp.label().into(), energy_kev: e, bundled_rp_nm: r, builtin_rp_nm: b, deviation_pct: (r - b) / b * 100.0 });
                }
            }
        }
    }
    let all_documented = !files.is_empty() && files.iter().all(|f| f.has_source && f.has_conditions && f.has_uncertainty);
    let any_hash_bad = files.iter().any(|f| f.hash_ok == Some(false));
    if any_hash_bad { notes.push("at least one bundled file no longer matches its recorded SHA-256".into()); }
    if !all_documented { notes.push("source, conditions and uncertainty are not documented for every file: the pack cannot be called calibrated".into()); }
    let evidence = if files.iter().all(|f| f.declared_status == "experimental") && all_documented && !any_hash_bad { Evidence::Experimental }
        else if all_documented && !any_hash_bad { Evidence::Reference } else { Evidence::Unverified };
    Ok(CalibrationStatus { evidence, files, implant_deviations: dev, notes })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sha256_matches_known_vectors() {
        assert_eq!(sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
        assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    }

    #[test]
    fn bundled_pack_is_unverified_and_its_implant_table_is_flagged() {
        let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("..");
        let st = audit(&root, &root.join("data/calibration")).expect("audit");
        assert_eq!(st.evidence, Evidence::Unverified);
        assert!(!st.allows_calibrated_label());
        assert!(st.files.iter().all(|f| f.hash_ok == Some(true)), "hashes identify the bytes: {:?}", st.files);
        assert!(st.max_implant_deviation_pct() > 15.0, "bundled boron range table disagrees with the built-in one: {:.1}%", st.max_implant_deviation_pct());
    }

    #[test]
    fn missing_manifest_is_an_error_not_a_pass() {
        assert!(audit(Path::new("."), Path::new("/definitely/not/here")).is_err());
    }
}
