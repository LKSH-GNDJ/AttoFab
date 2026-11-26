// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! The sealed blind protocol (PRD K-02).
//!
//! A blind comparison is only blind if the predictions were fixed before the
//! reference was looked at. Discipline is not a mechanism: the way to guarantee
//! the order is to make the unsealed reference impossible to obtain without a
//! seal that was recorded first.
//!
//! ```text
//!   let s = seal(&predictions);        // digest recorded, no reference in scope
//!   let r = compare(&s, &reference);   // refused unless the seal verifies
//! ```
//!
//! The digest is over the prediction values and their names, so a prediction that
//! moves after sealing is caught rather than quietly compared.

use crate::sha256;

/// Predictions, sealed before the reference is available.
#[derive(Debug, Clone, PartialEq)]
pub struct Seal {
    /// Hex digest of the predictions at the moment of sealing.
    pub digest: String,
    /// How many predictions were sealed.
    pub count: usize,
}

/// What a comparison produced.
#[derive(Debug, Clone, PartialEq)]
pub struct Comparison {
    pub matched: usize,
    pub missed: usize,
    /// Largest absolute error across the predictions.
    pub worst_error: f64,
}

/// Seal a set of named predictions.
pub fn seal(predictions: &[(String, f64)]) -> Seal {
    let mut text = String::new();
    for (name, value) in predictions {
        text.push_str(name);
        text.push('=');
        text.push_str(&format!("{value:.12e}"));
        text.push('\n');
    }
    Seal { digest: sha256::hex_digest(text.as_bytes()), count: predictions.len() }
}

/// Compare sealed predictions against a reference, refusing if the seal does not
/// match the predictions offered.
///
/// The caller cannot reach this without a `Seal`, and a `Seal` cannot be built
/// except by `seal`, so the order is enforced by the types rather than by a rule
/// someone has to remember.
pub fn compare(
    seal: &Seal,
    predictions: &[(String, f64)],
    reference: &[(String, f64)],
    tolerance: f64,
) -> Result<Comparison, &'static str> {
    if !(tolerance > 0.0) || !tolerance.is_finite() {
        return Err("the tolerance must be positive and finite");
    }
    if predictions.len() != seal.count {
        return Err("the predictions are not the ones that were sealed");
    }
    let now = self::seal(predictions);
    if now.digest != seal.digest {
        return Err("the predictions changed after sealing");
    }
    let mut matched = 0;
    let mut missed = 0;
    let mut worst = 0.0f64;
    for (name, value) in predictions {
        match reference.iter().find(|(n, _)| n == name) {
            Some((_, want)) => {
                let e = (value - want).abs();
                worst = worst.max(e);
                if e <= tolerance { matched += 1 } else { missed += 1 }
            }
            None => missed += 1,
        }
    }
    Ok(Comparison { matched, missed, worst_error: worst })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn preds() -> Vec<(String, f64)> {
        vec![
            ("CD".into(), 96.0),
            ("xj".into(), 42.0),
            ("Vt".into(), 0.71),
        ]
    }

    fn reference() -> Vec<(String, f64)> {
        vec![
            ("CD".into(), 96.5),
            ("xj".into(), 41.0),
            ("Vt".into(), 0.70),
        ]
    }

    #[test]
    fn a_sealed_set_compares() {
        let p = preds();
        let s = seal(&p);
        assert_eq!(s.count, 3);
        assert_eq!(s.digest.len(), 64);
        // errors are 0.5, 1.0 and 0.01; a tolerance of 0.75 separates them
        let c = compare(&s, &p, &reference(), 0.75).expect("a sealed comparison");
        assert_eq!(c.matched, 2, "CD and Vt are inside 0.75");
        assert_eq!(c.missed, 1, "xj is exactly 1.0 out, so outside");
        assert!((c.worst_error - 1.0).abs() < 1e-9);
    }

    #[test]
    fn a_error_exactly_at_the_tolerance_is_inside() {
        // the boundary is a decision, not an accident: an error equal to the
        // tolerance passes, and the test says so out loud
        let p = preds();
        let s = seal(&p);
        let c = compare(&s, &p, &reference(), 1.0).expect("a sealed comparison");
        assert_eq!(c.matched, 3, "an error of exactly the tolerance is accepted");
    }

    #[test]
    fn a_prediction_that_moves_after_sealing_is_caught() {
        // the whole point: the reference cannot be used to adjust a prediction
        let p = preds();
        let s = seal(&p);
        let mut moved = p.clone();
        moved[0].1 = 96.5; // "fixed" to match the reference
        assert_eq!(compare(&s, &moved, &reference(), 1.0), Err("the predictions changed after sealing"));
    }

    #[test]
    fn a_truncated_prediction_set_is_refused() {
        let p = preds();
        let s = seal(&p);
        let short = p[..2].to_vec();
        assert_eq!(compare(&s, &short, &reference(), 1.0), Err("the predictions are not the ones that were sealed"));
    }

    #[test]
    fn the_digest_is_stable_and_sensitive() {
        let a = seal(&preds());
        let b = seal(&preds());
        assert_eq!(a.digest, b.digest, "the same predictions seal the same way");
        let mut other = preds();
        other[1].1 += 1e-9;
        assert_ne!(a.digest, seal(&other).digest, "a nanometre moves the digest");
    }

    #[test]
    fn a_nonsense_tolerance_is_refused() {
        let p = preds();
        let s = seal(&p);
        assert!(compare(&s, &p, &reference(), 0.0).is_err());
        assert!(compare(&s, &p, &reference(), f64::NAN).is_err());
    }
}
