// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Metrology on the committed state.
//!
//! Three audit defects are fixed here: junction depth is measured from the
//! surface rather than summed from runs (F-07), step height is a difference
//! between adjacent columns rather than a percentile range that reads zero on a
//! blanket film (F-10), and critical dimension is taken from the top run rather
//! than any run containing the material (F-20). Ellipsometry handles an
//! arbitrary stack with dispersive indices (M-01, M-02).

pub mod optics;
pub mod overlay;

pub mod cdsem;
pub mod guard;
pub mod sims;
pub mod sites;

use atto_geom::column::{solid_top, ColumnGrid};
use atto_geom::svo::Material;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
/// Whether a measurement passed, failed or could not be judged.
pub enum Status {
    Pass,
    Warn,
    Fail,
    NotMeasured,
}

impl Status {
    /// The name a report uses.
    pub fn label(self) -> &'static str {
        match self {
            Status::Pass => "PASS",
            Status::Warn => "WARN",
            Status::Fail => "FAIL",
            Status::NotMeasured => "n/a",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
/// One measured quantity against its nominal and tolerance.
pub struct Measurement {
    pub name: String,
    pub unit: String,
    pub value: Option<f64>,
    pub nominal: f64,
    pub tolerance: f64,
    /// One-sigma uncertainty on the value, same unit.
    pub sigma: f64,
    pub status: Status,
    /// Why the quantity could not be measured, when it could not.
    pub reason: Option<String>,
}

impl Measurement {
    /// A measurement that could not be taken. It carries a reason rather than a
    /// zero (build PRD G-13).
    pub fn missing(name: &str, unit: &str, reason: &str) -> Self {
        Measurement {
            name: name.into(),
            unit: unit.into(),
            value: None,
            nominal: 0.0,
            tolerance: 0.0,
            sigma: 0.0,
            status: Status::NotMeasured,
            reason: Some(reason.into()),
        }
    }

    /// Judge a value against a nominal and a tolerance, with a guard band: a
    /// value within one sigma of the limit is a warning, not a pass.
    pub fn against(name: &str, unit: &str, value: f64, nominal: f64, tolerance: f64, sigma: f64) -> Self {
        let dev = (value - nominal).abs();
        let status = if dev > tolerance {
            Status::Fail
        } else if dev + sigma > tolerance {
            Status::Warn
        } else {
            Status::Pass
        };
        Measurement {
            name: name.into(),
            unit: unit.into(),
            value: Some(value),
            nominal,
            tolerance,
            sigma,
            status,
            reason: None,
        }
    }
}

/// Mean thickness of a material over the columns that carry it.
pub fn film_thickness_nm(g: &ColumnGrid, m: Material) -> Option<f64> {
    let mut sum = 0.0;
    let mut n = 0;
    for (k, col) in g.cols.iter().enumerate() {
        if !g.active[k] {
            continue;
        }
        let t: f64 = col.iter().filter(|r| r.material == m).map(|r| r.thickness()).sum();
        if t > 0.0 {
            sum += t;
            n += 1;
        }
    }
    if n == 0 { None } else { Some(sum / n as f64) }
}

/// The largest height difference between adjacent active columns. Zero on a
/// blanket film, which is the correct answer (audit F-10).
pub fn step_height_nm(g: &ColumnGrid) -> Option<f64> {
    let h = g.surface_heights();
    let mut best: f64 = 0.0;
    let mut any = false;
    for iy in 0..g.ny {
        for ix in 0..g.nx {
            let k = g.idx(ix, iy);
            if !g.active[k] {
                continue;
            }
            any = true;
            if ix + 1 < g.nx {
                let j = g.idx(ix + 1, iy);
                if g.active[j] {
                    best = best.max((h[k] - h[j]).abs());
                }
            }
            if iy + 1 < g.ny {
                let j = g.idx(ix, iy + 1);
                if g.active[j] {
                    best = best.max((h[k] - h[j]).abs());
                }
            }
        }
    }
    if any { Some(best) } else { None }
}

/// How far the surface dropped between two states: the etch depth (audit F-10
/// separates this from the step height).
pub fn etch_depth_nm(before: &ColumnGrid, after: &ColumnGrid) -> Option<f64> {
    let a = before.surface_heights();
    let b = after.surface_heights();
    if a.len() != b.len() {
        return None;
    }
    let mut sum = 0.0;
    let mut n = 0;
    for k in 0..a.len() {
        if before.active[k] && after.active[k] {
            sum += (a[k] - b[k]).max(0.0);
            n += 1;
        }
    }
    if n == 0 { None } else { Some(sum / n as f64) }
}

/// Widths of the runs of `m` that are exposed at the top of each column, along
/// the centre row. A buried blanket film does not satisfy this (audit F-20).
pub fn top_run_widths_nm(g: &ColumnGrid, m: Material) -> Vec<f64> {
    let row = g.ny / 2;
    let has = |ix: usize| {
        let k = g.idx(ix, row);
        g.active[k] && matches!(solid_top(&g.cols[k]), Some((_, mm)) if mm == m)
    };
    let mut out = Vec::new();
    let mut i = 0;
    while i < g.nx {
        if has(i) {
            let s = i;
            while i < g.nx && has(i) {
                i += 1;
            }
            if s > 0 && i < g.nx {
                out.push((i - s) as f64 * g.pitch);
            }
        } else {
            i += 1;
        }
    }
    out
}

/// Junction depth from the surface, given the deepest depth at which the net
/// doping exceeds a threshold in any column. This is the geometric counterpart
/// of the analytic profile: it takes a depth map, not a sum of runs.
pub fn junction_depth_from_map(depths_nm: &[Option<f64>]) -> Option<f64> {
    depths_nm.iter().filter_map(|d| *d).fold(None, |acc: Option<f64>, d| {
        Some(match acc {
            Some(a) => a.max(d),
            None => d,
        })
    })
}

#[cfg(test)]
mod width_tests {
    //! Audit F-20: a width must come from the material being exposed at the top
    //! of a column, not merely present somewhere inside it.
    use atto_units::{Mm, Nm, Um};
    use super::*;
    use atto_geom::column::{set_interval, ColumnGrid, ProcessWindow};
    use atto_geom::svo::{Material, Wafer};

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] };
        ColumnGrid::sample(&w, &win).unwrap()
    }

    #[test]
    fn a_blanket_film_has_no_line_width() {
        // resist everywhere, exposed at the top: no edges, so no width
        let mut g = grid();
        for k in 0..g.cols.len() {
            let top = atto_geom::column::solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top, top + 80.0, Material::Photoresist);
        }
        let w = top_run_widths_nm(&g, Material::Photoresist);
        assert!(w.is_empty(), "a blanket film has no measured width, got {w:?}");
    }

    #[test]
    fn a_buried_film_is_not_a_line() {
        // resist under an oxide cap: present, but not at the top
        let mut g = grid();
        for k in 0..g.cols.len() {
            let top = atto_geom::column::solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top, top + 80.0, Material::Photoresist);
            set_interval(&mut g.cols[k], top + 80.0, top + 160.0, Material::SiliconDioxide);
        }
        let w = top_run_widths_nm(&g, Material::Photoresist);
        assert!(w.is_empty(), "a buried film is not a line, got {w:?}");
    }

    #[test]
    fn exposed_lines_report_their_width() {
        let mut g = grid();
        // four adjacent columns of resist in the middle, every other row
        for k in 0..g.cols.len() {
            let ix = k % g.nx;
            if !(10..14).contains(&ix) {
                continue;
            }
            let top = atto_geom::column::solid_top(&g.cols[k]).map(|(t, _)| t).unwrap_or(g.z_bottom);
            set_interval(&mut g.cols[k], top, top + 80.0, Material::Photoresist);
        }
        let w = top_run_widths_nm(&g, Material::Photoresist);
        assert_eq!(w.len(), 1, "one line, got {w:?}");
        assert!((w[0] - 4.0 * g.pitch).abs() < 1e-9, "width {:.1} nm", w[0]);
    }
}

#[cfg(test)]
mod tests {
    use atto_units::{Mm, Nm, Um};
    use super::*;
    use atto_geom::column::{grow_film, set_interval, ProcessWindow};
    use atto_geom::svo::Wafer;

    fn grid() -> ColumnGrid {
        let w = Wafer::new_tile(Mm::new(300.0), Um::new(500.0), Nm::new(16.0), 7, (0.0, 0.0)).unwrap();
        ColumnGrid::sample(&w, &ProcessWindow { center: [0.0, 0.0], size: [512.0, 512.0] }).unwrap()
    }

    #[test]
    fn blanket_film_has_zero_step_height_but_a_real_thickness() {
        let g = grid();
        let (o, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(96.0), 0.0).unwrap();
        let t = film_thickness_nm(&o, Material::SiliconDioxide).unwrap();
        assert!((t - 96.0).abs() < 1e-9, "thickness {t}");
        let s = step_height_nm(&o).unwrap();
        assert!(s < 1e-9, "a blanket film has no step, got {s} nm");
    }

    #[test]
    fn a_trench_produces_a_step_and_an_etch_depth() {
        let g = grid();
        let (film, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(96.0), 0.0).unwrap();
        let mut cut = film.clone();
        for iy in 0..cut.ny {
            for ix in 10..13 {
                let k = cut.idx(ix, iy);
                let top = solid_top(&cut.cols[k]).unwrap().0;
                set_interval(&mut cut.cols[k], top - 96.0, top, Material::Air);
            }
        }
        let step = step_height_nm(&cut).unwrap();
        assert!((step - 96.0).abs() < 1e-6, "step {step}");
        let depth = etch_depth_nm(&film, &cut).unwrap();
        // 3 of 32 columns cut by 96 nm out of a 32x32 window: mean = 96*3/32 = 9
        assert!((depth - 9.0).abs() < 0.5, "mean depth over the window {depth}");
    }

    #[test]
    fn cd_is_taken_from_the_top_run_only() {
        let g = grid();
        // blanket oxide everywhere, then a resist pattern on top of it
        let (ox, _) = grow_film(&g, Material::SiliconDioxide, Nm::new(96.0), 0.0).unwrap();
        let mut pat = ox.clone();
        for iy in 0..pat.ny {
            for ix in 0..pat.nx {
                let k = pat.idx(ix, iy);
                if (4..8).contains(&ix) {
                    let top = solid_top(&pat.cols[k]).unwrap().0;
                    set_interval(&mut pat.cols[k], top, top + 64.0, Material::Photoresist);
                }
            }
        }
        let w = top_run_widths_nm(&pat, Material::Photoresist);
        assert_eq!(w.len(), 1, "one line: {w:?}");
        assert!((w[0] - 64.0).abs() < 1e-9, "line width {}", w[0]);
        // the oxide is buried, so it is not a top run and does not read as the window width
        assert!(top_run_widths_nm(&pat, Material::SiliconDioxide).is_empty());
    }

    #[test]
    fn missing_measurements_carry_a_reason() {
        let m = Measurement::missing("Junction depth", "nm", "no dopant in the window");
        assert_eq!(m.status, Status::NotMeasured);
        assert!(m.value.is_none());
        assert!(m.reason.unwrap().contains("no dopant"));
    }

    #[test]
    fn guard_band_turns_a_marginal_pass_into_a_warning() {
        let m = Measurement::against("t", "nm", 109.0, 100.0, 10.0, 2.0);
        assert_eq!(m.status, Status::Warn, "9 nm deviation with 2 nm sigma is inside the guard band");
        let ok = Measurement::against("t", "nm", 103.0, 100.0, 10.0, 2.0);
        assert_eq!(ok.status, Status::Pass);
        let bad = Measurement::against("t", "nm", 112.0, 100.0, 10.0, 2.0);
        assert_eq!(bad.status, Status::Fail);
    }

    #[test]
    fn junction_takes_the_deepest_column() {
        let d = [Some(120.0), None, Some(334.0), Some(200.0)];
        assert_eq!(junction_depth_from_map(&d), Some(334.0));
        assert_eq!(junction_depth_from_map(&[None, None]), None);
    }
}
