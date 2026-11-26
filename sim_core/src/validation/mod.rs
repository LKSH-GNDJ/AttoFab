//! Verification and validation helpers for numerical process models.
//!
//! These utilities report numerical evidence; they do not certify a model or
//! substitute for experimentally traceable calibration data.

pub mod blind;
pub mod conservation;
pub mod gci;

use serde::{Deserialize, Serialize};

/// Summary state of a verification run.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum CheckStatus {
    /// All configured numerical acceptance checks passed.
    Pass,
    /// At least one configured check failed.
    Fail,
    /// The check could not be completed because inputs or evidence were missing.
    Blocked,
}

/// A portable machine-readable validation finding.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Finding {
    /// Stable short check identifier.
    pub id: String,
    /// Pass/fail/blocked state.
    pub status: CheckStatus,
    /// Human-readable detail.
    pub message: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn status_serializes_predictably() {
        assert_eq!(
            serde_json::to_string(&CheckStatus::Pass).unwrap(),
            "\"pass\""
        );
    }
}
