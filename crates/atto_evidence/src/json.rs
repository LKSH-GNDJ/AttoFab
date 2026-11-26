// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A minimal JSON reader. The evidence crate has no dependencies, so that a
//! gate can be audited without auditing a dependency tree as well.

#[derive(Debug, Clone, PartialEq)]
/// A minimal JSON value: enough to read a project or a dataset, no more.
pub enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

impl Json {
    /// The value under a key, if this is an object that has one.
    pub fn get(&self, key: &str) -> Option<&Json> {
        match self {
            Json::Obj(m) => m.iter().find(|(k, _)| k == key).map(|(_, v)| v),
            _ => None,
        }
    }
    /// The number, if this is one.
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Json::Num(v) => Some(*v),
            _ => None,
        }
    }
    /// The string, if this is one.
    pub fn as_str(&self) -> Option<&str> {
        match self {
            Json::Str(s) => Some(s),
            _ => None,
        }
    }
    /// The elements, if this is an array.
    pub fn as_array(&self) -> Option<&[Json]> {
        match self {
            Json::Arr(v) => Some(v),
            _ => None,
        }
    }
    /// Whether this is the null value.
    pub fn is_null(&self) -> bool { matches!(self, Json::Null) }
    /// The boolean, if this is one.
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Json::Bool(b) => Some(*b),
            _ => None,
        }
    }
}

/// Parse a document, returning the reason rather than a partial tree.
pub fn parse(text: &str) -> Result<Json, String> {
    let b: Vec<char> = text.chars().collect();
    let mut i = 0usize;
    let v = value(&b, &mut i)?;
    skip_ws(&b, &mut i);
    if i != b.len() {
        return Err(format!("trailing characters at {i}"));
    }
    Ok(v)
}

fn skip_ws(b: &[char], i: &mut usize) {
    while *i < b.len() && b[*i].is_whitespace() {
        *i += 1;
    }
}

fn value(b: &[char], i: &mut usize) -> Result<Json, String> {
    skip_ws(b, i);
    if *i >= b.len() {
        return Err("unexpected end of input".into());
    }
    match b[*i] {
        '{' => object(b, i),
        '[' => array(b, i),
        '"' => Ok(Json::Str(string(b, i)?)),
        't' | 'f' | 'n' => literal(b, i),
        _ => number(b, i),
    }
}

fn literal(b: &[char], i: &mut usize) -> Result<Json, String> {
    let rest: String = b[*i..].iter().take(5).collect();
    for (word, v) in [("true", Json::Bool(true)), ("false", Json::Bool(false)), ("null", Json::Null)] {
        if rest.starts_with(word) {
            *i += word.len();
            return Ok(v);
        }
    }
    Err(format!("bad literal at {i}"))
}

fn number(b: &[char], i: &mut usize) -> Result<Json, String> {
    let start = *i;
    while *i < b.len() && (b[*i].is_ascii_digit() || "+-.eE".contains(b[*i])) {
        *i += 1;
    }
    let s: String = b[start..*i].iter().collect();
    s.parse::<f64>().map(Json::Num).map_err(|_| format!("bad number '{s}'"))
}

fn string(b: &[char], i: &mut usize) -> Result<String, String> {
    if *i >= b.len() {
        return Err("unexpected end of input in string".into());
    }
    if b[*i] != '"' {
        return Err("expected a string".into());
    }
    *i += 1;
    let mut out = String::new();
    while *i < b.len() {
        let c = b[*i];
        *i += 1;
        match c {
            '"' => return Ok(out),
            '\\' => {
                if *i >= b.len() {
                    return Err("unterminated escape".into());
                }
                let e = b[*i];
                *i += 1;
                out.push(match e {
                    'n' => '\n',
                    't' => '\t',
                    'r' => '\r',
                    'b' => '\u{8}',
                    'f' => '\u{c}',
                    '"' => '"',
                    '\\' => '\\',
                    '/' => '/',
                    'u' => {
                        let hex: String = b[*i..(*i + 4).min(b.len())].iter().collect();
                        *i += 4;
                        char::from_u32(u32::from_str_radix(&hex, 16).map_err(|_| "bad unicode escape".to_string())?)
                            .ok_or_else(|| "bad unicode scalar".to_string())?
                    }
                    other => return Err(format!("bad escape '\\{other}'")),
                });
            }
            _ => out.push(c),
        }
    }
    Err("unterminated string".into())
}

fn array(b: &[char], i: &mut usize) -> Result<Json, String> {
    *i += 1;
    let mut v = Vec::new();
    loop {
        skip_ws(b, i);
        if *i < b.len() && b[*i] == ']' {
            *i += 1;
            return Ok(Json::Arr(v));
        }
        v.push(value(b, i)?);
        skip_ws(b, i);
        if *i < b.len() && b[*i] == ',' {
            *i += 1;
        }
    }
}

fn object(b: &[char], i: &mut usize) -> Result<Json, String> {
    *i += 1;
    let mut m = Vec::new();
    loop {
        skip_ws(b, i);
        if *i < b.len() && b[*i] == '}' {
            *i += 1;
            return Ok(Json::Obj(m));
        }
        let k = string(b, i)?;
        skip_ws(b, i);
        if *i >= b.len() || b[*i] != ':' {
            return Err("expected ':'".into());
        }
        *i += 1;
        let v = value(b, i)?;
        m.push((k, v));
        skip_ws(b, i);
        if *i < b.len() && b[*i] == ',' {
            *i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_the_shapes_the_gate_needs() {
        let v = parse(r#"{"a": 1, "b": [1,2,3], "c": "x", "d": null, "e": true}"#).unwrap();
        assert_eq!(v.get("a").unwrap().as_f64(), Some(1.0));
        assert_eq!(v.get("b").unwrap().as_array().unwrap().len(), 3);
        assert_eq!(v.get("c").unwrap().as_str(), Some("x"));
        assert!(v.get("d").unwrap().is_null());
        assert_eq!(v.get("e").unwrap().as_bool(), Some(true));
    }

    #[test]
    fn handles_escapes_and_exponents() {
        let v = parse(r#"{"s": "a\"b\\c", "n": -1.5e-9}"#).unwrap();
        assert_eq!(v.get("s").unwrap().as_str(), Some("a\"b\\c"));
        assert!((v.get("n").unwrap().as_f64().unwrap() + 1.5e-9).abs() < 1e-20);
    }

    #[test]
    fn rejects_malformed_input() {
        assert!(parse("{").is_err());
        assert!(parse(r#"{"a": }"#).is_err());
        assert!(parse("[1,2").is_err());
        assert!(parse(r#"{"a":1} junk"#).is_err());
    }
}
