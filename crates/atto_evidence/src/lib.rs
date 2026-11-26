// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Evidence: the sign-off gate, calibration audit and blind validation
//! (build PRD R-01 to R-06, K-01 to K-06).
//!
//! This crate is deliberately separate from the physics. A gate that shares
//! code with the thing it is judging is not a gate, and every rule here exists
//! because the audit found a way for the old one to report a better status than
//! the evidence supported.

pub mod blind;
pub mod calibration;
pub mod exact;
pub mod gate;
pub mod ingest;
pub mod metrics;
pub mod json;
pub mod lint;
pub mod sealed;
pub mod sha256;

pub use blind::{BlindDataset, BlindReport};
pub use calibration::{CalibrationAudit, Evidence, HashStatus};
pub use gate::{evaluate, CheckStatus, Finding, GateInputs, GateReport, Readiness};
pub use sha256::hex_digest;

#[cfg(test)]
mod pack_tests {
    //! Tests over the data this repository actually ships. They run without an
    //! application checkout, because the manifest records its hashes.

    use crate::blind;
    use crate::calibration::{self, Evidence, HashStatus};
    use crate::hex_digest;
    use std::path::PathBuf;

    fn root() -> PathBuf {
        let mut r = std::env::current_dir().unwrap();
        while !r.join("data").is_dir() {
            assert!(r.pop(), "could not find the data directory from {:?}", std::env::current_dir());
        }
        r
    }

    #[test]
    fn the_bundled_pack_is_unmodified_from_its_recorded_hashes() {
        let r = root();
        let text = std::fs::read_to_string(r.join("data/calibration/manifest.json")).unwrap();
        let read = |p: &str| std::fs::read(r.join(p)).ok();
        let a = calibration::audit(&text, &read, &hex_digest, None).expect("the manifest parses");
        assert!(!a.files.is_empty());
        for f in &a.files {
            assert_eq!(f.hash, HashStatus::Verified, "{} is not the file the manifest records", f.path);
        }
        // and the honest classification of what that proves
        assert_eq!(a.evidence, Evidence::Unverified, "identity is not provenance");
    }

    #[test]
    fn every_shipped_blind_dataset_is_synthetic() {
        let r = root();
        let mut seen = 0;
        for e in std::fs::read_dir(r.join("data/blind_test")).unwrap().flatten() {
            let p = e.path();
            if p.extension().map(|x| x == "json").unwrap_or(false) {
                let text = std::fs::read_to_string(&p).unwrap();
                let ds = blind::load(&text).expect("every shipped dataset loads");
                assert!(ds.synthetic, "{} claims to be measured data", ds.name);
                assert!(ds.generated_from.is_none(), "{} claims a generator it does not have", ds.name);
                seen += 1;
            }
        }
        assert!(seen >= 5, "expected the bundled datasets, saw {seen}");
    }
}
