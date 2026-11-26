// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The command log and the project file (build PRD A-09, R-14).
//!
//! A session is not saved as a blob of state; it is saved as the list of
//! commands that produced it. Replaying that list rebuilds the session, which is
//! what makes a saved project reproducible: if a replay lands somewhere else,
//! something in the model changed.

use crate::{Evaluator, Session, StageSpec, UiError};

#[derive(Debug, Clone, PartialEq)]
/// A recorded change to the session.
pub enum Command {
    SetParam { stage: usize, key: String, value: f64 },
    SetEnabled { stage: usize, on: bool },
    Commit { stage: usize },
    Undo,
    Redo,
}

impl Command {
    /// Apply to a session. Commit needs the evaluator; the rest do not.
    pub fn apply(&self, s: &mut Session, eval: Option<Evaluator<'_>>) -> Result<(), UiError> {
        match self {
            Command::SetParam { stage, key, value } => s.set_param(*stage, key, *value),
            Command::SetEnabled { stage, on } => s.set_enabled(*stage, *on),
            Command::Commit { stage } => match eval {
                Some(e) => {
                    s.commit(*stage, e)?;
                    Ok(())
                }
                None => Err(UiError::NoSuchStage),
            },
            Command::Undo => s.undo(),
            Command::Redo => s.redo(),
        }
    }

    fn to_json(&self) -> String {
        match self {
            Command::SetParam { stage, key, value } => format!(
                r#"{{"op":"set_param","stage":{stage},"key":{},"value":{value}}}"#,
                quote(key)
            ),
            Command::SetEnabled { stage, on } => format!(r#"{{"op":"set_enabled","stage":{stage},"on":{on}}}"#),
            Command::Commit { stage } => format!(r#"{{"op":"commit","stage":{stage}}}"#),
            Command::Undo => r#"{"op":"undo"}"#.to_string(),
            Command::Redo => r#"{"op":"redo"}"#.to_string(),
        }
    }
}

fn quote(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

/// An append-only list of commands. Never edited, only extended.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Journal {
    pub commands: Vec<Command>,
}

impl Journal {
    /// An empty log.
    pub fn new() -> Self {
        Journal { commands: Vec::new() }
    }

    /// Append a command; the log is append-only.
    pub fn record(&mut self, c: Command) {
        self.commands.push(c);
    }

    /// How many commands the log holds.
    pub fn len(&self) -> usize {
        self.commands.len()
    }

    /// Whether anything has been recorded.
    pub fn is_empty(&self) -> bool {
        self.commands.is_empty()
    }

    /// Rebuild a session by replaying the log in order.
    pub fn replay(&self, eval: Evaluator<'_>) -> Result<Session, UiError> {
        let mut s = Session::new();
        for c in &self.commands {
            c.apply(&mut s, Some(eval))?;
        }
        Ok(s)
    }

    /// The project file. One JSON object, commands in order.
    pub fn to_json(&self) -> String {
        let body: Vec<String> = self.commands.iter().map(|c| c.to_json()).collect();
        format!("{{\"schema_version\":1,\"commands\":[{}]}}", body.join(","))
    }

    /// Parse a project file. A malformed file is an error, never an empty
    /// journal that would silently lose a session.
    pub fn from_json(text: &str) -> Result<Journal, String> {
        let v = crate::mini_json::parse(text)?;
        let version = v.get("schema_version").and_then(|x| x.as_f64()).ok_or("no schema_version")?;
        if version != 1.0 {
            return Err(format!("unsupported schema_version {version}"));
        }
        let arr = v.get("commands").and_then(|c| c.as_array()).ok_or("no commands array")?;
        let mut commands = Vec::with_capacity(arr.len());
        for c in arr {
            let op = c.get("op").and_then(|o| o.as_str()).ok_or("a command has no op")?;
            let stage = c.get("stage").and_then(|s| s.as_f64()).map(|s| s as usize);
            commands.push(match op {
                "set_param" => Command::SetParam {
                    stage: stage.ok_or("set_param needs a stage")?,
                    key: c.get("key").and_then(|k| k.as_str()).ok_or("set_param needs a key")?.to_string(),
                    value: c.get("value").and_then(|x| x.as_f64()).ok_or("set_param needs a value")?,
                },
                "set_enabled" => Command::SetEnabled {
                    stage: stage.ok_or("set_enabled needs a stage")?,
                    on: c.get("on").and_then(|b| b.as_bool()).ok_or("set_enabled needs on")?,
                },
                "commit" => Command::Commit { stage: stage.ok_or("commit needs a stage")? },
                "undo" => Command::Undo,
                "redo" => Command::Redo,
                other => return Err(format!("unknown op '{other}'")),
            });
        }
        Ok(Journal { commands })
    }

    /// A recipe in a form a human can read and diff.
    pub fn summary(&self) -> Vec<String> {
        self.commands
            .iter()
            .map(|c| match c {
                Command::SetParam { stage, key, value } => format!("stage {stage}: {key} = {value}"),
                Command::SetEnabled { stage, on } => format!("stage {stage}: enabled = {on}"),
                Command::Commit { stage } => format!("stage {stage}: commit"),
                Command::Undo => "undo".to_string(),
                Command::Redo => "redo".to_string(),
            })
            .collect()
    }
}

/// The recipe as it stands, so a project can be saved from a live session.
pub fn journal_of(stages: &[StageSpec]) -> Vec<(usize, String, f64)> {
    let mut out = Vec::new();
    for (i, s) in stages.iter().enumerate() {
        for (k, v) in &s.params {
            out.push((i, k.clone(), *v));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{measure_rows, StageResult};
    use atto_measure::Measurement;

    fn eval(_spec: &StageSpec, fp: u64) -> Result<StageResult, String> {
        Ok(StageResult {
            measurements: vec![Measurement::against("Film thickness", "nm", 96.0, 100.0, 10.0, 1.0)],
            model: "fake".into(),
            parameters: vec!["thickness".into()],
            inputs_fingerprint: fp,
        })
    }

    fn scripted() -> Journal {
        let mut j = Journal::new();
        j.record(Command::SetParam { stage: 1, key: "thickness".into(), value: 96.0 });
        j.record(Command::Commit { stage: 1 });
        j.record(Command::SetParam { stage: 1, key: "thickness".into(), value: 140.0 });
        j.record(Command::Undo);
        j.record(Command::SetEnabled { stage: 2, on: false });
        j
    }

    #[test]
    fn a_replay_reproduces_the_session() {
        let j = scripted();
        let a = j.replay(&eval).expect("replay");
        let b = j.replay(&eval).expect("replay");
        assert_eq!(a.fingerprint(), b.fingerprint(), "the same log gives the same session");
        assert_eq!(a.stages()[1].get("thickness"), Some(96.0), "undo was replayed");
        assert!(!a.stages()[2].enabled);
    }

    #[test]
    fn a_live_session_and_a_replayed_one_agree() {
        let mut live = Session::new();
        live.set_param(1, "thickness", 96.0).unwrap();
        live.commit(1, &eval).unwrap();
        live.set_param(1, "thickness", 140.0).unwrap();
        live.undo().unwrap();
        live.set_enabled(2, false).unwrap();
        let replayed = scripted().replay(&eval).unwrap();
        assert_eq!(live.fingerprint(), replayed.fingerprint());
    }

    #[test]
    fn the_project_file_round_trips() {
        let j = scripted();
        let text = j.to_json();
        let back = Journal::from_json(&text).expect("parses");
        assert_eq!(j, back, "commands survive the round trip");
        assert_eq!(
            j.replay(&eval).unwrap().fingerprint(),
            back.replay(&eval).unwrap().fingerprint()
        );
    }

    #[test]
    fn a_malformed_project_file_is_an_error_not_an_empty_session() {
        assert!(Journal::from_json("{").is_err());
        assert!(Journal::from_json(r#"{"commands":[]}"#).is_err(), "no schema_version");
        assert!(Journal::from_json(r#"{"schema_version":2,"commands":[]}"#).is_err(), "wrong version");
        assert!(Journal::from_json(r#"{"schema_version":1}"#).is_err(), "no commands");
        assert!(Journal::from_json(r#"{"schema_version":1,"commands":[{"op":"fly"}]}"#).is_err(), "unknown op");
        assert!(Journal::from_json(r#"{"schema_version":1,"commands":[{"op":"commit"}]}"#).is_err(), "commit without a stage");
    }

    #[test]
    fn an_empty_journal_is_a_fresh_session() {
        let j = Journal::new();
        assert!(j.is_empty());
        let s = j.replay(&eval).unwrap();
        assert_eq!(s.fingerprint(), Session::new().fingerprint());
    }

    #[test]
    fn the_summary_is_readable() {
        let lines = scripted().summary();
        assert_eq!(lines[0], "stage 1: thickness = 96");
        assert_eq!(lines[1], "stage 1: commit");
        assert_eq!(lines[3], "undo");
    }

    #[test]
    fn a_commit_without_an_evaluator_is_an_error() {
        let mut s = Session::new();
        assert_eq!(Command::Commit { stage: 1 }.apply(&mut s, None), Err(UiError::NoSuchStage));
    }

    #[test]
    fn saving_from_a_live_session_lists_its_parameters() {
        let mut s = Session::new();
        s.set_param(1, "thickness", 96.0).unwrap();
        s.set_param(1, "coverage", 0.3).unwrap();
        let rows = journal_of(s.stages());
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().any(|(i, k, v)| *i == 1 && k == "thickness" && *v == 96.0));
    }

    #[test]
    fn a_replay_that_commits_produces_measurements() {
        let s = scripted().replay(&eval).unwrap();
        let rows = measure_rows(s.result(1).expect("stage 1 committed"));
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].value, Some(96.0));
    }
}
