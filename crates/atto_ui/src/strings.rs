// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! User-facing strings, in one table (PRD G-27).
//!
//! A word that appears in the code and again in the interface is a word that will
//! eventually disagree with itself. Everything the operator reads comes from here,
//! and `missing()` reports anything a panel would have to invent.

use crate::console::Severity;
use crate::{StageId, StageStatus};

/// One entry: a stable key and the text shown.
pub type Entry = (&'static str, &'static str);

/// The status words, as constants so the enum can name them without a lookup.
pub const STATUS_NOT_RUN: &str = "Not run";
/// The blocked status.
pub const STATUS_BLOCKED: &str = "Blocked";
/// The committed status.
pub const STATUS_COMMITTED: &str = "Committed";
/// The stale status.
pub const STATUS_STALE: &str = "Stale";
/// The failed status.
pub const STATUS_FAILED: &str = "Failed";

/// Every user-facing string.
pub const TABLE: &[Entry] = &[
    // stage names
    ("stage.prep", "Substrate prep"),
    ("stage.deposition", "Deposition"),
    ("stage.lithography", "Lithography"),
    ("stage.etching", "Etching"),
    ("stage.doping", "Doping"),
    ("stage.anneal", "Anneal"),
    ("stage.metrology", "Metrology"),
    ("stage.oxidation", "Oxidation"),
    // status words
    ("status.notrun", STATUS_NOT_RUN),
    ("status.blocked", STATUS_BLOCKED),
    ("status.committed", STATUS_COMMITTED),
    ("status.stale", STATUS_STALE),
    ("status.failed", STATUS_FAILED),
    // severities
    ("severity.info", "info"),
    ("severity.warning", "warning"),
    ("severity.error", "error"),
    // the few verbs the interface uses
    ("action.commit", "Commit"),
    ("action.undo", "Undo"),
    ("action.redo", "Redo"),
    ("action.compare", "Compare"),
    ("action.run", "Run the whole recipe"),
];

/// The text for a key, if the table has one.
pub fn get(key: &str) -> Option<&'static str> {
    TABLE.iter().find(|(k, _)| *k == key).map(|(_, v)| *v)
}

/// The key for a stage.
pub fn stage_key(id: StageId) -> &'static str {
    match id {
        StageId::Prep => "stage.prep",
        StageId::Deposition => "stage.deposition",
        StageId::Lithography => "stage.lithography",
        StageId::Etching => "stage.etching",
        StageId::Doping => "stage.doping",
        StageId::Anneal => "stage.anneal",
        StageId::Metrology => "stage.metrology",
        StageId::Oxidation => "stage.oxidation",
    }
}

/// The key for a status.
pub fn status_key(s: StageStatus) -> &'static str {
    match s {
        StageStatus::NotRun => "status.notrun",
        StageStatus::Blocked => "status.blocked",
        StageStatus::Committed => "status.committed",
        StageStatus::Stale => "status.stale",
        StageStatus::Failed => "status.failed",
    }
}

/// The key for a severity.
pub fn severity_key(s: Severity) -> &'static str {
    match s {
        Severity::Info => "severity.info",
        Severity::Warning => "severity.warning",
        Severity::Error => "severity.error",
    }
}

/// Keys the interface would ask for that the table does not have.
///
/// Empty means every word the panels need is in the table, which is what G-27
/// asks for; anything else is a word the code would have to invent.
pub fn missing() -> Vec<&'static str> {
    let mut out = Vec::new();
    for id in StageId::ALL {
        let k = stage_key(id);
        if get(k).is_none() {
            out.push(k);
        }
    }
    for s in [
        StageStatus::NotRun,
        StageStatus::Blocked,
        StageStatus::Committed,
        StageStatus::Stale,
        StageStatus::Failed,
    ] {
        let k = status_key(s);
        if get(k).is_none() {
            out.push(k);
        }
    }
    for s in [Severity::Info, Severity::Warning, Severity::Error] {
        let k = severity_key(s);
        if get(k).is_none() {
            out.push(k);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_word_a_panel_needs_is_in_the_table() {
        // this is G-27 itself: no panel invents its own text
        assert!(missing().is_empty(), "the table lacks: {:?}", missing());
    }

    #[test]
    fn the_table_has_no_duplicate_keys_or_empty_text() {
        let mut keys: Vec<&str> = TABLE.iter().map(|(k, _)| *k).collect();
        let n = keys.len();
        keys.sort();
        keys.dedup();
        assert_eq!(keys.len(), n, "two entries share a key");
        for (k, v) in TABLE {
            assert!(!k.is_empty(), "an entry has no key");
            assert!(!v.trim().is_empty(), "entry {k} has no text");
        }
    }

    #[test]
    fn the_table_agrees_with_the_enums_it_serves() {
        // Stage names are still owned by StageId, so this comparison is real.
        for id in StageId::ALL {
            assert_eq!(get(stage_key(id)), Some(id.label()), "stage {id:?}");
        }
        // Severity and status now *read* the table, so comparing them to it would
        // be circular. What can still be checked is that the words are the ones a
        // reader expects, spelled out here once.
        assert_eq!(get("severity.error"), Some("error"));
        assert_eq!(get("status.failed"), Some("Failed"));
        assert_eq!(crate::StageStatus::Failed.word(), "Failed");
    }

    #[test]
    fn an_unknown_key_is_none_not_an_empty_string() {
        assert_eq!(get("no.such.key"), None);
    }
}
