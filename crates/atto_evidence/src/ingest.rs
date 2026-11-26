// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Measured-data ingestion with a provenance chain (K-08).
//!
//! Every dataset in this tree until now was synthetic. The moment a real measurement enters,
//! the question stops being "does the model match" and becomes "is this curve what it claims
//! to be" — and a curve with no source, no instrument and no date is not evidence, it is a
//! number with a shape.
//!
//! The rule here is deliberately strict: a dataset that cannot say where it came from is
//! **refused**, not accepted with a blank field. A gate that passes an unattributed curve is
//! worse than no gate, because the curve then carries the authority of having passed.

/// Where a measured dataset came from.
#[derive(Debug, Clone, PartialEq)]
pub struct Provenance {
    /// The file or record it was read from.
    pub source: String,
    /// Who took it.
    pub operator: String,
    /// The instrument, with its model if known.
    pub instrument: String,
    /// When, as `YYYY-MM-DD`.
    pub date_iso: String,
    /// The units of the y axis, as written on the instrument.
    pub units: String,
}

/// A measured series: x and y pairs, plus where they came from.
#[derive(Debug, Clone, PartialEq)]
pub struct Dataset {
    pub name: String,
    pub provenance: Provenance,
    pub series: Vec<(f64, f64)>,
}

/// Why a dataset was refused.
#[derive(Debug, Clone, PartialEq)]
pub struct IngestError {
    pub what: String,
}

impl std::fmt::Display for IngestError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "dataset refused: {}", self.what)
    }
}

fn need(field: &str, value: &str) -> Result<(), IngestError> {
    if value.trim().is_empty() {
        return Err(IngestError { what: format!("{field} is missing — an unattributed dataset is not evidence") });
    }
    Ok(())
}

/// A `YYYY-MM-DD` date that is a real calendar date, not merely shaped like one.
fn check_date(d: &str) -> Result<(), IngestError> {
    let b = d.as_bytes();
    if b.len() != 10 || b[4] != b'-' || b[7] != b'-' {
        return Err(IngestError { what: format!("date '{d}' is not YYYY-MM-DD") });
    }
    let year: i32 = d[0..4].parse().map_err(|_| IngestError { what: format!("date '{d}' has no year") })?;
    let month: u32 = d[5..7].parse().map_err(|_| IngestError { what: format!("date '{d}' has no month") })?;
    let day: u32 = d[8..10].parse().map_err(|_| IngestError { what: format!("date '{d}' has no day") })?;
    if !(1..=12).contains(&month) {
        return Err(IngestError { what: format!("date '{d}' has month {month}") });
    }
    let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        _ => if leap { 29 } else { 28 },
    };
    if day < 1 || day > days {
        return Err(IngestError { what: format!("date '{d}' has day {day} in a {days}-day month") });
    }
    Ok(())
}

/// Validate a dataset before anything is allowed to use it.
pub fn validate(d: &Dataset) -> Result<(), IngestError> {
    need("dataset name", &d.name)?;
    need("source", &d.provenance.source)?;
    need("operator", &d.provenance.operator)?;
    need("instrument", &d.provenance.instrument)?;
    need("units", &d.provenance.units)?;
    need("date", &d.provenance.date_iso)?;
    check_date(&d.provenance.date_iso)?;
    if d.series.len() < 2 {
        return Err(IngestError { what: "fewer than two points is not a curve".into() });
    }
    for (i, (x, y)) in d.series.iter().enumerate() {
        if !x.is_finite() || !y.is_finite() {
            return Err(IngestError { what: format!("point {i} is not finite") });
        }
    }
    Ok(())
}

/// A one-line provenance string for a report, so a number can always be traced back.
pub fn provenance_line(d: &Dataset) -> String {
    format!(
        "{} — {}, {}, {} ({})",
        d.name, d.provenance.instrument, d.provenance.operator, d.provenance.date_iso, d.provenance.units
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn good() -> Dataset {
        Dataset {
            name: "SIMS B11".into(),
            provenance: Provenance {
                source: "wafer-17/sims_b11.csv".into(),
                operator: "L. Gundoju".into(),
                instrument: "CAMECA IMS 7f".into(),
                date_iso: "2026-11-17".into(),
                units: "cm^-3".into(),
            },
            series: vec![(0.0, 1.0e19), (10.0, 5.0e18), (20.0, 1.0e17)],
        }
    }

    #[test]
    fn a_labelled_dataset_is_accepted_and_carries_its_provenance() {
        let d = good();
        assert!(validate(&d).is_ok());
        let line = provenance_line(&d);
        assert!(line.contains("CAMECA"), "the instrument is in the line: {line}");
        assert!(line.contains("2026-11-17"), "and the date: {line}");
    }

    #[test]
    fn an_unattributed_dataset_is_refused() {
        for field in ["source", "operator", "instrument", "units", "date"] {
            let mut d = good();
            match field {
                "source" => d.provenance.source.clear(),
                "operator" => d.provenance.operator = "  ".into(),
                "instrument" => d.provenance.instrument.clear(),
                "units" => d.provenance.units.clear(),
                _ => d.provenance.date_iso.clear(),
            }
            assert!(validate(&d).is_err(), "a dataset with no {field} must be refused");
        }
    }

    #[test]
    fn a_date_that_is_not_a_date_is_refused() {
        for bad in ["17-11-2026", "2026-13-01", "2026-02-30", "2026-11-00", "not a date"] {
            let mut d = good();
            d.provenance.date_iso = bad.into();
            assert!(validate(&d).is_err(), "'{bad}' must be refused");
        }
        let mut leap = good();
        leap.provenance.date_iso = "2024-02-29".into();
        assert!(validate(&leap).is_ok(), "2024 is a leap year");
        let mut not_leap = good();
        not_leap.provenance.date_iso = "2026-02-29".into();
        assert!(validate(&not_leap).is_err(), "2026 is not");
    }

    #[test]
    fn a_curve_that_is_not_a_curve_is_refused() {
        let mut d = good();
        d.series.truncate(1);
        assert!(validate(&d).is_err(), "one point");
        let mut nan = good();
        nan.series[1].1 = f64::NAN;
        assert!(validate(&nan).is_err(), "not a number");
    }

    #[test]
    fn a_dataset_with_no_name_is_refused() {
        let mut d = good();
        d.name.clear();
        assert!(validate(&d).is_err(), "an unnamed dataset cannot be cited");
    }
}
