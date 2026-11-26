// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A source scan for silent error swallowing (build PRD R-06).
//!
//! The audit found four places where an error became a value: `.ok().flatten()`
//! on a solver result, a `continue` where a file failed to load, a default in
//! place of a failure. A rule that lives only in a review comment does not
//! survive; this one is checked by the test suite over this repository's own
//! source.

#[derive(Debug, Clone, PartialEq)]
/// One place the lint found a silently swallowed error.
pub struct Violation {
    pub line: usize,
    pub pattern: String,
    pub text: String,
}

/// Patterns that turn an error into a non-error value.
pub const FORBIDDEN: [&str; 3] = [".ok().flatten()", "unwrap_or(0.0)", "unwrap_or_default()"];

/// Scan one source file, stopping at the test module so that test scaffolding
/// (which may legitimately unwrap) is not judged as production behaviour.
pub fn scan(text: &str) -> Vec<Violation> {
    let mut out = Vec::new();
    let lines: Vec<&str> = text.lines().collect();
    for (i, line) in lines.iter().enumerate() {
        if line.trim_start().starts_with("#[cfg(test)]") {
            break;
        }
        // A documented exemption. Silencing the lint by deleting the rule is how
        // a rule stops being a rule; an exemption on the line, with a reason, is
        // visible in review and survives.
        let exempt = line.contains("lint:allow")
            || (i > 0 && lines[i - 1].contains("lint:allow"));
        if exempt {
            continue;
        }
        let code = line.split("//").next().unwrap_or("");
        for p in FORBIDDEN {
            if code.contains(p) {
                out.push(Violation { line: i + 1, pattern: p.to_string(), text: line.trim().to_string() });
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_a_swallowed_error() {
        let src = "fn f() -> f64 {\n    let v = g().ok().flatten().unwrap_or(0.0);\n    v\n}\n";
        let v = scan(src);
        assert_eq!(v.len(), 2, "{v:?}");
        assert_eq!(v[0].line, 2);
    }

    #[test]
    fn ignores_comments_and_tests() {
        let src = "// never use .ok().flatten() here\nfn g() -> f64 { 1.0 }\n#[cfg(test)]\nmod t { fn h() { let _ = x.unwrap_or(0.0); } }\n";
        assert!(scan(src).is_empty(), "{:?}", scan(src));
    }

    #[test]
    fn a_documented_exemption_is_honoured_and_an_undocumented_one_is_not() {
        let exempt_same_line = "let v = opt.unwrap_or(0.0); // lint:allow zero is the physical default here\n";
        assert!(scan(exempt_same_line).is_empty(), "an exemption on the line is honoured");
        let exempt_above = "// lint:allow an empty cell is deliberate\nlet s = x.unwrap_or_default();\n";
        assert!(scan(exempt_above).is_empty(), "an exemption above is honoured");
        let bare = "let v = opt.unwrap_or(0.0);\n";
        assert_eq!(scan(bare).len(), 1, "no exemption, still a violation");
    }

    #[test]
    fn scans_this_workspace_clean() {
        // The repository is the subject: no production file in the domain crates
        // may turn an error into a value.
        // rustc is invoked directly here, so CARGO_MANIFEST_DIR is not set:
        // walk up from the working directory to the workspace root.
        let mut root = std::env::current_dir().unwrap();
        while !root.join("crates").is_dir() {
            assert!(root.pop(), "could not find the workspace root from {:?}", std::env::current_dir());
        }
        let mut scanned = 0;
        let mut violations = Vec::new();
        // The domain crates are the subject. `atto_cli` is a shell: it formats
        // for a human and may legitimately show a default.
        for crate_dir in std::fs::read_dir(root.join("crates")).unwrap().flatten() {
            if crate_dir.file_name() == std::ffi::OsStr::new("atto_cli") {
                continue;
            }
            let src = crate_dir.path().join("src");
            if !src.is_dir() {
                continue;
            }
            let mut stack = vec![src];
            while let Some(dir) = stack.pop() {
                for e in std::fs::read_dir(&dir).unwrap().flatten() {
                    let p = e.path();
                    if p.is_dir() {
                        stack.push(p);
                    } else if p.extension().map(|x| x == "rs").unwrap_or(false) {
                        // this module holds the patterns as literals; scanning it
                        // would always match itself
                        if p.file_name() == Some(std::ffi::OsStr::new("lint.rs")) {
                            continue;
                        }
                        // `src/bin` holds shells: they format for a human
                        if p.components().any(|c| c.as_os_str() == "bin") {
                            continue;
                        }
                        let text = std::fs::read_to_string(&p).unwrap();
                        scanned += 1;
                        for v in scan(&text) {
                            violations.push(format!("{}:{} {}", p.display(), v.line, v.pattern));
                        }
                    }
                }
            }
        }
        assert!(scanned > 8, "the scan should see the whole workspace, saw {scanned} files");
        assert!(violations.is_empty(), "silent error swallowing found:\n{}", violations.join("\n"));
    }
}
