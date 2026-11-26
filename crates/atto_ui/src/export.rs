// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Exporting measurements (build PRD X-07).
//!
//! Two formats, both plain text: CSV for a spreadsheet and JSON for a script.
//! A measurement that was not taken exports as an empty cell and a reason,
//! never as a zero, so an exported table cannot be mistaken for data.

use crate::{MeasureRow, Provenance};

fn esc_csv(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Header plus one row per measurement.
pub fn to_csv(rows: &[MeasureRow]) -> String {
    let mut out = String::from("name,value,unit,nominal,tolerance,sigma,status,model,reason\n");
    for r in rows {
        let value = match r.value {
            Some(v) => format!("{v:.6}"),
            None => String::new(),
        };
        // an undeclared uncertainty exports as an empty cell, deliberately
        let sigma = match r.provenance.uncertainty {
            Some(s) => format!("{s:.6}"),
            None => String::new(),
        };
        out.push_str(&format!(
            "{},{},{},{:.6},{:.6},{},{},{},{}\n",
            esc_csv(&r.label),
            value,
            esc_csv(&r.unit),
            r.nominal,
            r.tolerance,
            sigma,
            esc_csv(r.status.label()),
            esc_csv(&r.provenance.model),
            esc_csv(r.reason.as_deref().unwrap_or("")),
        ));
    }
    out
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

/// A single JSON object with a measurements array.
pub fn to_json(rows: &[MeasureRow], state_id: u64) -> String {
    let body: Vec<String> = rows
        .iter()
        .map(|r| {
            let value = match r.value {
                Some(v) => format!("{v:.6}"),
                None => "null".to_string(),
            };
            format!(
                r#"{{"name":{},"value":{},"unit":{},"nominal":{:.6},"tolerance":{:.6},"status":{},"model":{},"reason":{}}}"#,
                quote(&r.label),
                value,
                quote(&r.unit),
                r.nominal,
                r.tolerance,
                quote(r.status.label()),
                quote(&r.provenance.model),
                match &r.reason {
                    Some(x) => quote(x),
                    None => "null".to_string(),
                }
            )
        })
        .collect();
    format!(r#"{{"schema_version":1,"state_id":"{state_id:016x}","measurements":[{}]}}"#, body.join(","))
}

/// Provenance for a whole export, so a file records where its numbers came from.
pub fn provenance_line(p: &Provenance) -> String {
    format!("{} | {} | fp {:016x}", p.model, p.parameters.join(", "), p.inputs_fingerprint)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::measure_rows;
    use crate::StageResult;
    use atto_measure::Measurement;

    fn rows() -> Vec<MeasureRow> {
        let res = StageResult {
            measurements: vec![
                Measurement::against("Film thickness", "nm", 96.0, 100.0, 10.0, 1.0),
                Measurement::missing("Junction depth", "nm", "no dopant in this window"),
            ],
            model: "atto_geom::grow_film".into(),
            parameters: vec!["thickness".into()],
            inputs_fingerprint: 0x1234,
        };
        measure_rows(&res)
    }

    #[test]
    fn the_csv_has_a_header_and_one_row_per_measurement() {
        let csv = to_csv(&rows());
        let lines: Vec<&str> = csv.lines().collect();
        assert_eq!(lines.len(), 3, "header plus two measurements");
        assert!(lines[0].starts_with("name,value,unit"));
        assert!(lines[1].contains("Film thickness"));
        assert!(lines[2].contains("no dopant"));
    }

    #[test]
    fn an_unmeasured_row_exports_an_empty_cell_not_a_zero() {
        let csv = to_csv(&rows());
        let line = csv.lines().nth(2).unwrap();
        let cells: Vec<&str> = line.split(',').collect();
        assert_eq!(cells[1], "", "the value cell is empty: {line}");
        assert!(line.contains("no dopant"), "and the reason is carried");
    }

    #[test]
    fn a_name_with_a_comma_is_quoted() {
        let mut r = rows();
        r[0].label = "Thickness, film".into();
        let csv = to_csv(&r);
        assert!(csv.contains("\"Thickness, film\""), "{csv}");
    }

    #[test]
    fn the_json_parses_back_and_carries_the_state_id() {
        let json = to_json(&rows(), 0xabcdef);
        let v = crate::mini_json::parse(&json).expect("valid JSON");
        assert_eq!(v.get("schema_version").and_then(|x| x.as_f64()), Some(1.0));
        assert_eq!(v.get("state_id").and_then(|x| x.as_str()), Some("0000000000abcdef"));
        let arr = v.get("measurements").and_then(|m| m.as_array()).expect("array");
        assert_eq!(arr.len(), 2);
        assert_eq!(arr[0].get("name").and_then(|x| x.as_str()), Some("Film thickness"));
        // a missing measurement is null, not zero
        assert_eq!(arr[1].get("value"), Some(&crate::mini_json::Json::Null));
    }

    #[test]
    fn an_empty_export_is_still_valid() {
        let csv = to_csv(&[]);
        assert_eq!(csv.lines().count(), 1, "header only");
        let json = to_json(&[], 0);
        assert!(crate::mini_json::parse(&json).is_ok());
    }

    #[test]
    fn the_provenance_line_names_the_model_and_the_fingerprint() {
        let r = rows();
        let line = provenance_line(&r[0].provenance);
        assert!(line.contains("atto_geom::grow_film"));
        assert!(line.contains("thickness"));
        assert!(line.contains("0000000000001234"));
    }
}
