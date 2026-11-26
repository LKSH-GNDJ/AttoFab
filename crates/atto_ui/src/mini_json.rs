// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A minimal JSON reader for the project file. `atto_ui` has no dependencies, so
//! a project file can be opened by anything that can read text.

#[derive(Debug, Clone, PartialEq)]
/// A minimal JSON value, enough to read a project file.
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
    /// The boolean, if this is one.
    pub fn as_bool(&self) -> Option<bool> {
        match self {
            Json::Bool(b) => Some(*b),
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
}

/// Parse a document, returning the reason rather than a partial tree.
pub fn parse(text: &str) -> Result<Json, String> {
    let b: Vec<char> = text.chars().collect();
    let mut i = 0;
    let v = value(&b, &mut i)?;
    ws(&b, &mut i);
    if i != b.len() {
        return Err(format!("trailing characters at {i}"));
    }
    Ok(v)
}

fn ws(b: &[char], i: &mut usize) {
    while *i < b.len() && b[*i].is_whitespace() {
        *i += 1;
    }
}

fn value(b: &[char], i: &mut usize) -> Result<Json, String> {
    ws(b, i);
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
    for (w, v) in [("true", Json::Bool(true)), ("false", Json::Bool(false)), ("null", Json::Null)] {
        if rest.starts_with(w) {
            *i += w.len();
            return Ok(v);
        }
    }
    Err(format!("bad literal at {i}"))
}

fn number(b: &[char], i: &mut usize) -> Result<Json, String> {
    let s0 = *i;
    while *i < b.len() && (b[*i].is_ascii_digit() || "+-.eE".contains(b[*i])) {
        *i += 1;
    }
    let s: String = b[s0..*i].iter().collect();
    s.parse::<f64>().map(Json::Num).map_err(|_| format!("bad number '{s}'"))
}

fn string(b: &[char], i: &mut usize) -> Result<String, String> {
    if *i >= b.len() || b[*i] != '"' {
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
                    '"' => '"',
                    '\\' => '\\',
                    '/' => '/',
                    other => return Err(format!("bad escape '\\{other}'")),
                });
            }
            c => out.push(c),
        }
    }
    Err("unterminated string".into())
}

fn array(b: &[char], i: &mut usize) -> Result<Json, String> {
    *i += 1;
    let mut v = Vec::new();
    loop {
        ws(b, i);
        if *i < b.len() && b[*i] == ']' {
            *i += 1;
            return Ok(Json::Arr(v));
        }
        v.push(value(b, i)?);
        ws(b, i);
        if *i < b.len() && b[*i] == ',' {
            *i += 1;
        }
    }
}

fn object(b: &[char], i: &mut usize) -> Result<Json, String> {
    *i += 1;
    let mut m = Vec::new();
    loop {
        ws(b, i);
        if *i < b.len() && b[*i] == '}' {
            *i += 1;
            return Ok(Json::Obj(m));
        }
        let k = string(b, i)?;
        ws(b, i);
        if *i >= b.len() || b[*i] != ':' {
            return Err("expected ':'".into());
        }
        *i += 1;
        let v = value(b, i)?;
        m.push((k, v));
        ws(b, i);
        if *i < b.len() && b[*i] == ',' {
            *i += 1;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_parse_to_their_value() {
        // the module had no tests at all, and the project file depends on it
        let v = parse("{\"thickness\": 96.5}").expect("parses");
        assert_eq!(v.get("thickness").and_then(|x| x.as_f64()), Some(96.5));
        assert_eq!(parse("42").unwrap().as_f64(), Some(42.0));
        assert_eq!(parse("-7").unwrap().as_f64(), Some(-7.0));
        assert_eq!(parse("1.5e3").unwrap().as_f64(), Some(1500.0));
    }

    #[test]
    fn strings_booleans_and_nulls_parse() {
        assert_eq!(parse("\"hello\"").unwrap().as_str(), Some("hello"));
        assert_eq!(parse("true").unwrap().as_bool(), Some(true));
        assert_eq!(parse("false").unwrap().as_bool(), Some(false));
        assert!(parse("null").is_ok());
    }

    #[test]
    fn arrays_and_objects_nest() {
        let v = parse("{\"a\": [1, 2, 3]}").expect("parses");
        let arr = v.get("a").and_then(|x| x.as_array()).expect("an array");
        assert_eq!(arr.len(), 3);
        assert_eq!(arr[2].as_f64(), Some(3.0));
    }

    #[test]
    fn malformed_input_is_an_error_not_a_default() {
        // a parser that returns an empty object on failure makes a corrupt project
        // file look like a fresh one
        for bad in ["{", "{\"a\":}", "tru", "[1,", "\"unterminated"] {
            assert!(parse(bad).is_err(), "{bad:?} parsed");
        }
    }

    #[test]
    fn a_missing_key_is_none_not_zero() {
        let v = parse("{\"a\": 1}").expect("parses");
        assert!(v.get("b").is_none());
        assert_eq!(v.get("a").and_then(|x| x.as_f64()), Some(1.0));
    }
}
