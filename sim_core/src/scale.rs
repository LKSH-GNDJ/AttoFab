//! Simulation-scale advisor.
//!
//! The octree is built from *cubic* cells, so one cell size serves both lateral features and film
//! thickness. A layout with 50 nm and 5 µm features, and films of 100 nm, therefore needs a cell chosen from
//! the smallest feature **and** the thinnest film, and the window is capped at [`crate::geom::MAX_GRID_AXIS`]
//! cells. This module turns those two numbers into a cell size, a tile depth and a window.

use crate::geom::MAX_GRID_AXIS;

#[derive(Debug, Clone, PartialEq)]
pub struct ScaleAdvice { pub cell_nm: f64, pub tile_levels: u32, pub window_nm: f64, pub notes: Vec<String> }

const NICE: [f64; 23] = [4.0, 5.0, 8.0, 10.0, 12.5, 16.0, 20.0, 25.0, 32.0, 40.0, 50.0, 64.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 320.0, 400.0, 500.0, 800.0, 1000.0];

/// `feature_nm`: typical smallest lateral feature (e.g. the 10th percentile of the layer);
/// `thinnest_film_nm`: the thinnest film or resist you need to resolve vertically (use infinity if none).
pub fn advise(feature_nm: f64, thinnest_film_nm: f64) -> ScaleAdvice {
    let mut notes = Vec::new();
    let feature = if feature_nm.is_finite() && feature_nm > 0.0 { feature_nm } else { 100.0 };
    let film = if thinnest_film_nm.is_finite() && thinnest_film_nm > 0.0 { thinnest_film_nm } else { f64::INFINITY };
    let want = (feature / 4.0).min(film / 2.0);
    let cell = NICE.iter().rev().find(|c| **c <= want).copied().unwrap_or(NICE[0]);
    if want < NICE[0] { notes.push(format!("the smallest feature ({feature:.0} nm) or film ({film:.0} nm) is finer than the 4 nm minimum cell; it will be under-resolved")); }
    let window = MAX_GRID_AXIS as f64 * cell * 0.95;
    let mut levels = 6u32;
    while (cell * (1u64 << levels) as f64) < window / 0.95 && levels < 12 { levels += 1; }
    if feature / cell < 4.0 { notes.push(format!("features of {feature:.0} nm span only {:.1} cells", feature / cell)); }
    if film.is_finite() && film / cell < 2.0 { notes.push(format!("a {film:.0} nm film spans only {:.1} cells: thickness will be quantised", film / cell)); }
    ScaleAdvice { cell_nm: cell, tile_levels: levels, window_nm: (window / cell).floor() * cell, notes }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_follows_the_finer_of_feature_and_film_and_the_window_fits_the_grid_cap() {
        let a = advise(50.0, 100.0);           // 50 nm features: feature/4 = 12.5
        assert_eq!(a.cell_nm, 12.5);
        let b = advise(3333.0, 100.0);         // µm features, thin film: the film sets the cell
        assert_eq!(b.cell_nm, 50.0);
        let c = advise(3333.0, f64::INFINITY); // no film constraint: features set the cell
        assert_eq!(c.cell_nm, 800.0);
        for x in [&a, &b, &c] {
            assert!(x.window_nm / x.cell_nm <= MAX_GRID_AXIS as f64 + 1e-9);
            assert!((x.cell_nm * (1u64 << x.tile_levels) as f64) >= x.window_nm, "{x:?}");
        }
        assert!(advise(f64::NAN, f64::NAN).cell_nm >= 4.0, "bad input falls back to a safe default");
    }
}
