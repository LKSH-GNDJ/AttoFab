//! Interactive Parameter Table Grid — a spreadsheet-like matrix of numeric
//! inputs. Editing a cell fires an async recompute event rather than
//! mutating physics state directly, so the UI thread never blocks on the
//! recomputation.

use serde::{Deserialize, Serialize};
use sim_core::physics::Stage;

/// A single editable cell in the parameter table.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParamCell {
    pub key: String,
    pub label: String,
    pub value: f64,
    pub unit: String,
    pub min: Option<f64>,
    pub max: Option<f64>,
}

/// One stage's row group of parameter cells.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StageRow {
    pub stage: Stage_,
    pub cells: Vec<ParamCell>,
}

// Stage isn't Serialize in sim_core (it's a small Copy enum used for graph
// keys); mirror it locally as a serializable label so the UI layer can
// (de)serialize table state independently of the engine's internal enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Stage_ {
    SubstratePrep,
    Deposition,
    Lithography,
    Etching,
    Doping,
    FinalClean,
    Metrology,
    Cmp,
}

impl From<Stage> for Stage_ {
    fn from(s: Stage) -> Self {
        match s {
            Stage::SubstratePrep => Stage_::SubstratePrep,
            Stage::Deposition => Stage_::Deposition,
            Stage::Lithography => Stage_::Lithography,
            Stage::Etching => Stage_::Etching,
            Stage::Doping => Stage_::Doping,
            Stage::FinalClean => Stage_::FinalClean,
            Stage::Metrology => Stage_::Metrology, Stage::Cmp => Stage_::Cmp,
        }
    }
}

/// An edit event emitted when the user changes a cell — the platform shell
/// turns this into a `sim_core::engine::StageJob` and submits it to the
/// background worker.
#[derive(Debug, Clone)]
pub struct CellEditEvent {
    pub stage: Stage,
    pub key: String,
    pub new_value: f64,
}

/// The full parameter table state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ParameterTable {
    pub rows: Vec<StageRow>,
}

impl ParameterTable {
    /// Seed a table with the default parameters for every stage, using the
    /// same defaults as each `physics::*Inputs::default()`.
    pub fn default_table() -> Self {
        Self {
            rows: vec![
                StageRow {
                    stage: Stage_::SubstratePrep,
                    cells: vec![
                        cell("gamma_sg", "\u{3b3} SG", 60.0, "mN/m"),
                        cell("gamma_sl", "\u{3b3} SL", 20.0, "mN/m"),
                        cell("gamma_lg", "\u{3b3} LG", 72.8, "mN/m"),
                    ],
                },
                StageRow {
                    stage: Stage_::Deposition,
                    cells: vec![
                        cell("t_film", "Film thickness", 1e-6, "m"),
                        cell("bow_radius", "Bow radius", 10.0, "m"),
                    ],
                },
                StageRow {
                    stage: Stage_::Lithography,
                    cells: vec![
                        cell("k1", "k1", 0.3, ""),
                        cell("lambda", "\u{3bb}", 13.5e-9, "m"),
                        cell("na", "NA", 0.33, ""),
                    ],
                },
                StageRow {
                    stage: Stage_::Etching,
                    cells: vec![
                        cell("r_lateral", "R lateral", 5.0, "nm/min"),
                        cell("r_vertical", "R vertical", 50.0, "nm/min"),
                    ],
                },
                StageRow {
                    stage: Stage_::Doping,
                    cells: vec![
                        cell("dose_phi", "Dose \u{3a6}", 1e15, "ions/cm\u{b2}"),
                        cell("projected_range_nm", "Rp", 100.0, "nm"),
                    ],
                },
                StageRow {
                    stage: Stage_::Metrology,
                    cells: vec![
                        cell("incidence_deg", "Incidence angle", 70.0, "deg"),
                    ],
                },
            ],
        }
    }

    /// Apply an edit in place and return the updated cell, if found.
    pub fn apply_edit(&mut self, stage: Stage, key: &str, new_value: f64) -> Option<&ParamCell> {
        let target: Stage_ = stage.into();
        for row in self.rows.iter_mut() {
            if row.stage == target {
                if let Some(c) = row.cells.iter_mut().find(|c| c.key == key) {
                    c.value = new_value;
                    return Some(c);
                }
            }
        }
        None
    }
}

fn cell(key: &str, label: &str, value: f64, unit: &str) -> ParamCell {
    ParamCell {
        key: key.to_string(),
        label: label.to_string(),
        value,
        unit: unit.to_string(),
        min: None,
        max: None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_table_has_a_row_per_stage() {
        let table = ParameterTable::default_table();
        assert_eq!(table.rows.len(), 6);
    }

    #[test]
    fn apply_edit_updates_the_matching_cell() {
        let mut table = ParameterTable::default_table();
        let updated = table.apply_edit(Stage::Etching, "r_lateral", 12.5);
        assert_eq!(updated.unwrap().value, 12.5);
    }
}
