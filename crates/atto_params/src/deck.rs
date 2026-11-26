// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A recipe deck: a text a process engineer can write, diff and review (W-09).
//!
//! A JSON recipe is machine-readable and human-hostile: every review is a diff of punctuation, and
//! a one-value change moves three lines. The format here is the one a fab actually uses — a stage
//! header, then `key value` lines under it:
//!
//! ```text
//! # a ten-line recipe
//! stage substrate_prep
//!   rca_dilution 1.0
//!   hf_dip true
//! stage deposition
//!   material sio2
//!   thickness_nm 96
//! stage lithography
//!   na 0.85
//!   dose_ratio 2.0
//! ```
//!
//! The parser is deliberately strict and the printer is deliberately plain: a deck that goes in
//! and comes out unchanged is the property that makes a diff trustworthy, and it is the test.

/// One stage and its parameters, in the order they were written.
#[derive(Debug, Clone, PartialEq)]
pub struct DeckStage {
    pub name: String,
    /// Key/value pairs, kept in source order so a round trip preserves them.
    pub params: Vec<(String, String)>,
}

/// A parsed deck.
#[derive(Debug, Clone, PartialEq)]
pub struct Deck {
    /// Free-text comments, preserved so a round trip is exact.
    pub header: Vec<String>,
    pub stages: Vec<DeckStage>,
}

/// Why a deck could not be parsed.
#[derive(Debug, Clone, PartialEq)]
pub struct DeckError {
    /// One-based line number, so an editor can jump to it.
    pub line: usize,
    pub what: String,
}

impl std::fmt::Display for DeckError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "deck line {}: {}", self.line, self.what)
    }
}

/// Parse a deck.
pub fn parse(text: &str) -> Result<Deck, DeckError> {
    let mut header = Vec::new();
    let mut stages: Vec<DeckStage> = Vec::new();
    for (i, raw) in text.lines().enumerate() {
        let n = i + 1;
        let line = raw.trim_end();
        let trimmed = line.trim();
        if trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with('#') {
            if stages.is_empty() {
                header.push(trimmed.to_string());
            }
            continue;
        }
        let indented = line.starts_with(' ') || line.starts_with('\t');
        let mut parts = trimmed.split_whitespace();
        let first = parts.next().unwrap_or("");
        if !indented {
            if first != "stage" {
                return Err(DeckError { line: n, what: format!("expected 'stage <name>', found '{first}'") });
            }
            let name = parts.next().ok_or(DeckError { line: n, what: "a stage needs a name".into() })?;
            if parts.next().is_some() {
                return Err(DeckError { line: n, what: "a stage header takes one name".into() });
            }
            stages.push(DeckStage { name: name.to_string(), params: Vec::new() });
        } else {
            let Some(stage) = stages.last_mut() else {
                return Err(DeckError { line: n, what: "a parameter before any stage".into() });
            };
            let key = first.to_string();
            let value = parts.collect::<Vec<_>>().join(" ");
            if value.is_empty() {
                return Err(DeckError { line: n, what: format!("'{key}' has no value") });
            }
            if stage.params.iter().any(|(k, _)| *k == key) {
                return Err(DeckError { line: n, what: format!("'{key}' is set twice in stage '{}'", stage.name) });
            }
            stage.params.push((key, value));
        }
    }
    if stages.is_empty() {
        return Err(DeckError { line: 1, what: "the deck has no stages".into() });
    }
    Ok(Deck { header, stages })
}

/// Print a deck back to text.
pub fn to_text(deck: &Deck) -> String {
    let mut out = String::new();
    for h in &deck.header {
        out.push_str(h);
        out.push('\n');
    }
    for s in &deck.stages {
        out.push_str(&format!("stage {}\n", s.name));
        for (k, v) in &s.params {
            out.push_str(&format!("  {k} {v}\n"));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    const TEN_LINES: &str = "# a ten-line recipe\nstage substrate_prep\n  rca_dilution 1.0\n  hf_dip true\nstage deposition\n  material sio2\n  thickness_nm 96\nstage lithography\n  na 0.85\n  dose_ratio 2.0\n";

    #[test]
    fn a_ten_line_deck_parses_to_three_stages() {
        let d = parse(TEN_LINES).expect("a deck");
        assert_eq!(d.stages.len(), 3, "three stages");
        assert_eq!(d.stages[0].name, "substrate_prep");
        assert_eq!(d.stages[0].params.len(), 2);
        assert_eq!(d.stages[1].params[0], ("material".to_string(), "sio2".to_string()));
        assert_eq!(d.stages[2].params[1], ("dose_ratio".to_string(), "2.0".to_string()));
    }

    #[test]
    fn a_deck_round_trips_unchanged() {
        // W-09's property: what goes in comes out, so a diff means something
        let d = parse(TEN_LINES).expect("a deck");
        let again = to_text(&d);
        assert_eq!(again, TEN_LINES, "round trip is exact");
        let d2 = parse(&again).expect("a deck");
        assert_eq!(d, d2, "and parsing it again gives the same structure");
    }

    #[test]
    fn comments_and_blank_lines_are_preserved_or_ignored() {
        let with_blanks = "# header\n\nstage a\n  x 1\n\n# a note\nstage b\n  y 2\n";
        let d = parse(with_blanks).expect("a deck");
        assert_eq!(d.stages.len(), 2, "blank lines and mid-deck comments do not make stages");
        assert_eq!(d.header, vec!["# header".to_string()], "the leading comment is kept");
    }

    #[test]
    fn a_parameter_value_may_contain_spaces() {
        let d = parse("stage a\n  note two words here\n").expect("a deck");
        assert_eq!(d.stages[0].params[0].1, "two words here");
    }

    #[test]
    fn a_malformed_deck_names_the_line() {
        assert_eq!(parse("stage a\n  x\n").unwrap_err().line, 2, "a value is missing");
        assert_eq!(parse("  x 1\n").unwrap_err().line, 1, "a parameter before any stage");
        assert_eq!(parse("stages a\n").unwrap_err().line, 1, "a misspelt header");
        assert_eq!(parse("stage a\n  x 1\n  x 2\n").unwrap_err().line, 3, "a key set twice");
        assert_eq!(parse("stage\n").unwrap_err().line, 1, "a header with no name");
        assert!(parse("# only a comment\n").is_err(), "a deck with no stages is not a recipe");
    }
}
