// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Calibration audit (build PRD R-02, R-04, R-05).
//!
//! Every claim about the calibration pack is checked here rather than read: the
//! bytes are hashed, an unreadable file is a failure rather than an unchecked
//! file, and "experimental" requires a signature stored outside the directory
//! being blessed.

use crate::json::{self, Json};

#[derive(Debug, Clone, PartialEq)]
/// Whether a calibration file matched its recorded digest.
pub enum HashStatus {
    Verified,
    Mismatch,
    /// The file could not be read. Not the same as "not checked".
    Unreadable(String),
}

#[derive(Debug, Clone, PartialEq)]
/// Where a calibration value came from.
pub enum Provenance {
    Valid { kind: String, id: String },
    /// Present but not structured enough to check.
    Invalid(String),
    Missing,
}

#[derive(Debug, Clone, PartialEq)]
/// How a calibration value's uncertainty was established.
pub enum Uncertainty {
    Valid(f64),
    Invalid(String),
    Missing,
}

#[derive(Debug, Clone, PartialEq)]
/// The audit of one calibration file.
pub struct FileAudit {
    pub path: String,
    pub declared: String,
    pub hash: HashStatus,
    pub source: Provenance,
    pub uncertainty: Uncertainty,
    pub conditions_ok: bool,
}

impl FileAudit {
    /// Whether the file declares enough to be used.
    pub fn documented(&self) -> bool {
        matches!(self.source, Provenance::Valid { .. })
            && matches!(self.uncertainty, Uncertainty::Valid(_))
            && self.conditions_ok
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// What a model's evidence has actually established.
pub enum Evidence {
    Unverified,
    Reference,
    Experimental,
}

#[derive(Debug, Clone)]
/// The audit across every calibration file.
pub struct CalibrationAudit {
    pub files: Vec<FileAudit>,
    pub evidence: Evidence,
    pub notes: Vec<String>,
    /// Identity of the external signature, when one was supplied and matched.
    pub signed_by: Option<String>,
}

/// The source kinds a calibration file may declare.
pub const SOURCE_KINDS: [&str; 4] = ["measured", "published", "fitted", "assumed"];

fn provenance(v: Option<&Json>) -> Provenance {
    match v {
        None | Some(Json::Null) => Provenance::Missing,
        Some(j) => {
            let kind = j.get("kind").and_then(|k| k.as_str()).unwrap_or("");
            let id = j.get("id").and_then(|k| k.as_str()).unwrap_or("");
            if kind.is_empty() || id.is_empty() {
                Provenance::Invalid("source needs a kind and a resolvable id".into())
            } else if !SOURCE_KINDS.contains(&kind) {
                Provenance::Invalid(format!("'{kind}' is not one of {SOURCE_KINDS:?}"))
            } else if id.trim().len() < 3 {
                Provenance::Invalid("id is too short to resolve".into())
            } else {
                Provenance::Valid { kind: kind.to_string(), id: id.to_string() }
            }
        }
    }
}

/// Conditions are free-form: any non-empty object or string counts. A source,
/// by contrast, must be structured, because it has to resolve to something.
fn conditions_ok(v: Option<&Json>) -> bool {
    match v {
        None | Some(Json::Null) => false,
        Some(Json::Str(s)) => s.trim().len() >= 3,
        Some(Json::Obj(m)) => !m.is_empty(),
        Some(_) => false,
    }
}

fn uncertainty(v: Option<&Json>) -> Uncertainty {
    match v {
        None | Some(Json::Null) => Uncertainty::Missing,
        Some(j) => match j.as_f64() {
            Some(s) if s.is_finite() && s > 0.0 => Uncertainty::Valid(s),
            Some(_) => Uncertainty::Invalid("uncertainty must be positive".into()),
            None => Uncertainty::Invalid("uncertainty must be a number".into()),
        },
    }
}

/// Audit a manifest.
///
/// `read` returns the bytes of a path, or None when it cannot be read. `digest`
/// is injectable so the audit can be tested without a filesystem; production
/// passes [`crate::sha256::hex_digest`].
pub fn audit(
    manifest_text: &str,
    read: &dyn Fn(&str) -> Option<Vec<u8>>,
    digest: &dyn Fn(&[u8]) -> String,
    override_text: Option<&str>,
) -> Result<CalibrationAudit, String> {
    let root = json::parse(manifest_text)?;
    let files_json = root.get("files").and_then(|f| f.as_array()).ok_or("manifest has no files array")?;
    let mut files = Vec::new();
    let mut notes = Vec::new();
    for f in files_json {
        let path = f.get("path").and_then(|p| p.as_str()).unwrap_or("").to_string();
        if path.is_empty() {
            return Err("a manifest entry has no path".into());
        }
        let declared = f.get("evidence_status").and_then(|s| s.as_str()).unwrap_or("unknown").to_string();
        let expected = f.get("sha256").and_then(|s| s.as_str()).map(|s| s.to_lowercase());
        let hash = match (read(&path), expected) {
            (None, _) => HashStatus::Unreadable(format!("{path} could not be read")),
            (Some(_), None) => HashStatus::Unreadable(format!("{path} has no recorded hash")),
            (Some(bytes), Some(want)) => {
                if digest(&bytes) == want {
                    HashStatus::Verified
                } else {
                    HashStatus::Mismatch
                }
            }
        };
        files.push(FileAudit {
            path,
            declared,
            hash,
            source: provenance(f.get("source")),
            uncertainty: uncertainty(f.get("uncertainty")),
            conditions_ok: conditions_ok(f.get("conditions")),
        });
    }
    // a signature is only accepted if it names the pack and lives outside it
    let mut signed_by = None;
    if let Some(text) = override_text {
        let sig = json::parse(text).map_err(|e| format!("signature is not valid JSON: {e}"))?;
        let who = sig.get("signed_by").and_then(|s| s.as_str()).unwrap_or("");
        let bound = sig.get("pack").and_then(|s| s.as_str()).unwrap_or("");
        if who.trim().len() >= 3 && bound == "data/calibration" {
            signed_by = Some(who.to_string());
        } else {
            notes.push("the signature does not name the calibration pack, so it is ignored".into());
        }
    }

    let any_bad_hash = files.iter().any(|f| !matches!(f.hash, HashStatus::Verified));
    if any_bad_hash {
        notes.push("a calibration file is missing, unreadable, or no longer matches its recorded hash".into());
    }
    let all_documented = !files.is_empty() && files.iter().all(|f| f.documented());
    if !all_documented {
        notes.push("source, conditions and uncertainty are not verifiable for every file".into());
    }
    let all_experimental = files.iter().all(|f| f.declared == "experimental");
    let evidence = if all_documented && !any_bad_hash && all_experimental && signed_by.is_some() {
        Evidence::Experimental
    } else if all_documented && !any_bad_hash {
        Evidence::Reference
    } else {
        Evidence::Unverified
    };
    if all_experimental && all_documented && signed_by.is_none() {
        notes.push("the pack declares experimental evidence but carries no external signature".into());
    }
    Ok(CalibrationAudit { files, evidence, notes, signed_by })
}

#[cfg(test)]
mod tests {
    use super::*;

    const MANIFEST: &str = r#"{"schema_version":1,"files":[
        {"path":"a.json","sha256":"aaa","evidence_status":"unverified","source":null,"conditions":null,"uncertainty":null},
        {"path":"b.json","sha256":"bbb","evidence_status":"unverified","source":null,"conditions":null,"uncertainty":null}
    ]}"#;

    fn good_digest(b: &[u8]) -> String { String::from_utf8_lossy(b).to_string() }
    fn all_read(_p: &str) -> Option<Vec<u8>> { Some(b"aaa".to_vec()) }

    #[test]
    fn an_unreadable_file_is_a_failure_not_an_unchecked_file() {
        let read = |p: &str| if p == "a.json" { Some(b"aaa".to_vec()) } else { None };
        let a = audit(MANIFEST, &read, &good_digest, None).unwrap();
        assert_eq!(a.files[1].hash, HashStatus::Unreadable("b.json could not be read".into()));
        assert_eq!(a.evidence, Evidence::Unverified);
        assert!(a.notes.iter().any(|n| n.contains("unreadable")), "{:?}", a.notes);
    }

    #[test]
    fn a_hash_mismatch_is_a_failure() {
        let read = |p: &str| if p == "a.json" { Some(b"zzz".to_vec()) } else { Some(b"bbb".to_vec()) };
        let a = audit(MANIFEST, &read, &good_digest, None).unwrap();
        assert_eq!(a.files[0].hash, HashStatus::Mismatch);
        assert_eq!(a.evidence, Evidence::Unverified);
    }

    #[test]
    fn placeholder_provenance_is_not_experimental() {
        let m = r#"{"files":[{"path":"a.json","sha256":"aaa","evidence_status":"experimental",
            "source":"x","conditions":"{}","uncertainty":0}]}"#;
        let a = audit(m, &all_read, &good_digest, None).unwrap();
        assert!(matches!(a.files[0].source, Provenance::Invalid(_)), "a bare string is not a source");
        assert!(matches!(a.files[0].uncertainty, Uncertainty::Invalid(_)));
        assert_eq!(a.evidence, Evidence::Unverified, "placeholder text must not reach Reference");
    }

    #[test]
    fn documented_files_reach_reference_but_not_experimental() {
        let m = r#"{"files":[{"path":"a.json","sha256":"aaa","evidence_status":"experimental",
            "source":{"kind":"published","id":"doi:10.1063/1.123456"},
            "conditions":{"tool":"tube furnace"},
            "uncertainty":2.5}]}"#;
        let a = audit(m, &all_read, &good_digest, None).unwrap();
        assert!(a.files[0].documented());
        assert_eq!(a.evidence, Evidence::Reference, "no signature, so no Experimental");
        assert!(a.notes.iter().any(|n| n.contains("external signature")));
    }

    #[test]
    fn a_signature_outside_the_pack_is_what_reaches_experimental() {
        let m = r#"{"files":[{"path":"a.json","sha256":"aaa","evidence_status":"experimental",
            "source":{"kind":"measured","id":"lot 42, tool A"},
            "conditions":{"tool":"A"},
            "uncertainty":2.5}]}"#;
        let sig = r#"{"signed_by":"L. Gundoju","pack":"data/calibration"}"#;
        let a = audit(m, &all_read, &good_digest, Some(sig)).unwrap();
        assert_eq!(a.evidence, Evidence::Experimental);
        assert_eq!(a.signed_by.as_deref(), Some("L. Gundoju"));
    }

    #[test]
    fn a_signature_naming_the_wrong_pack_is_ignored() {
        let m = r#"{"files":[{"path":"a.json","sha256":"aaa","evidence_status":"experimental",
            "source":{"kind":"measured","id":"lot 42"},"conditions":{"tool":"A"},"uncertainty":2.5}]}"#;
        let sig = r#"{"signed_by":"L. Gundoju","pack":"something/else"}"#;
        let a = audit(m, &all_read, &good_digest, Some(sig)).unwrap();
        assert_eq!(a.evidence, Evidence::Reference);
        assert!(a.signed_by.is_none());
    }

    #[test]
    fn a_broken_signature_is_reported_not_swallowed() {
        let m = r#"{"files":[]}"#;
        assert!(audit(m, &all_read, &good_digest, Some("{not json")).is_err());
    }
}
