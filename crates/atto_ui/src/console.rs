// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The console: a log with severity, source and a filter (PRD G-17).
//!
//! The requirement has one part that is easy to get wrong: the console must be the
//! only place an error surfaces. That means an error cannot be silently dropped by
//! a filter — filtering hides a line from the view, it does not discard it. The log
//! is append-only and bounded; the *view* is what the filter narrows.

use crate::Colour;

/// How serious a line is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// Something happened; nothing to do.
    Info,
    /// Something looks wrong but the run continued.
    Warning,
    /// The run stopped.
    Error,
}

impl Severity {
    /// The word a panel draws.
    pub fn word(self) -> &'static str {
        match self {
            Severity::Info => "info",
            Severity::Warning => "warning",
            Severity::Error => "error",
        }
    }

    /// The colour a panel draws it in.
    pub fn colour(self) -> Colour {
        match self {
            Severity::Info => Colour::Grey,
            Severity::Warning => Colour::Amber,
            Severity::Error => Colour::Red,
        }
    }
}

/// One line in the log.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Line {
    pub severity: Severity,
    /// Where it came from: a stage name, or a module.
    pub source: String,
    pub message: String,
}

/// What the view is narrowed to.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct Filter {
    /// Show only this severity and above. `None` shows everything.
    pub at_least: Option<Severity>,
    /// Show only this source, if set.
    pub source: Option<String>,
}

impl Filter {
    /// Whether a line passes.
    pub fn accepts(&self, line: &Line) -> bool {
        if let Some(min) = self.at_least {
            if line.severity < min {
                return false;
            }
        }
        if let Some(src) = &self.source {
            if &line.source != src {
                return false;
            }
        }
        true
    }
}

/// The log and its view.
#[derive(Debug, Clone, PartialEq)]
pub struct Console {
    lines: Vec<Line>,
    filter: Filter,
    /// The most lines kept. Older ones are dropped, newest kept.
    limit: usize,
}

impl Default for Console {
    fn default() -> Self {
        Console::new(500)
    }
}

impl Console {
    /// A console that keeps at most `limit` lines.
    pub fn new(limit: usize) -> Self {
        Console { lines: Vec::new(), filter: Filter::default(), limit: limit.max(1) }
    }

    /// Append a line. The log is append-only; the limit drops the oldest.
    pub fn push(&mut self, severity: Severity, source: &str, message: &str) {
        self.lines.push(Line {
            severity,
            source: source.to_string(),
            message: message.to_string(),
        });
        while self.lines.len() > self.limit {
            self.lines.remove(0);
        }
    }

    /// Every line, filtered or not.
    pub fn len(&self) -> usize {
        self.lines.len()
    }

    /// Whether the log is empty.
    pub fn is_empty(&self) -> bool {
        self.lines.is_empty()
    }

    /// The filter in force.
    pub fn filter(&self) -> &Filter {
        &self.filter
    }

    /// Narrow the view.
    pub fn set_filter(&mut self, filter: Filter) {
        self.filter = filter;
    }

    /// The lines the view shows.
    pub fn visible(&self) -> Vec<&Line> {
        self.lines.iter().filter(|l| self.filter.accepts(l)).collect()
    }

    /// How many lines the filter hides.
    pub fn hidden(&self) -> usize {
        self.lines.len() - self.visible().len()
    }

    /// Every error ever logged, filtered or not.
    ///
    /// This is the method that makes "the console is the only place an error
    /// surfaces" checkable: an error cannot be lost by narrowing the view.
    pub fn errors(&self) -> Vec<&Line> {
        self.lines.iter().filter(|l| l.severity == Severity::Error).collect()
    }

    /// Whether anything has gone wrong.
    pub fn has_errors(&self) -> bool {
        self.lines.iter().any(|l| l.severity == Severity::Error)
    }

    /// The sources that have logged something, in first-seen order.
    pub fn sources(&self) -> Vec<String> {
        let mut out: Vec<String> = Vec::new();
        for l in &self.lines {
            if !out.contains(&l.source) {
                out.push(l.source.clone());
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn console() -> Console {
        let mut c = Console::new(50);
        c.push(Severity::Info, "prep", "wafer loaded");
        c.push(Severity::Warning, "etch", "rate below the table");
        c.push(Severity::Error, "etch", "the plasma model rejected the pressure");
        c.push(Severity::Info, "litho", "exposure complete");
        c
    }

    #[test]
    fn a_filter_narrows_the_view_without_discarding_a_line() {
        // the requirement: filtering hides, it does not drop
        let mut c = console();
        assert_eq!(c.len(), 4);
        c.set_filter(Filter { at_least: Some(Severity::Warning), source: None });
        assert_eq!(c.visible().len(), 2, "warning and error");
        assert_eq!(c.hidden(), 2);
        assert_eq!(c.len(), 4, "the log itself is unchanged");
    }

    #[test]
    fn an_error_cannot_be_hidden_from_the_error_list() {
        // this is what makes the console the only place an error surfaces
        let mut c = console();
        c.set_filter(Filter { at_least: Some(Severity::Error), source: Some("litho".into()) });
        assert!(c.visible().is_empty(), "the filter hides everything");
        assert_eq!(c.errors().len(), 1, "but the error is still there");
        assert!(c.has_errors());
    }

    #[test]
    fn the_filter_matches_a_source_exactly() {
        let mut c = console();
        c.set_filter(Filter { at_least: None, source: Some("etch".into()) });
        let v = c.visible();
        assert_eq!(v.len(), 2);
        assert!(v.iter().all(|l| l.source == "etch"));
    }

    #[test]
    fn the_log_is_bounded_and_keeps_the_newest() {
        let mut c = Console::new(3);
        for i in 0..5 {
            c.push(Severity::Info, "loop", &format!("line {i}"));
        }
        assert_eq!(c.len(), 3, "the limit holds");
        assert_eq!(c.visible()[2].message, "line 4", "the newest survives");
    }

    #[test]
    fn severity_orders_from_info_to_error() {
        assert!(Severity::Info < Severity::Warning);
        assert!(Severity::Warning < Severity::Error);
        assert_eq!(Severity::Error.colour(), Colour::Red);
        assert_eq!(Severity::Warning.colour(), Colour::Amber);
        assert!(!Severity::Info.word().is_empty());
    }

    #[test]
    fn the_sources_are_listed_in_first_seen_order() {
        let c = console();
        assert_eq!(c.sources(), vec!["prep".to_string(), "etch".to_string(), "litho".to_string()]);
    }
}
