//! Optical Lithography Solver (Phase 2, simplified).
//!
//! Mask geometry is projected straight down onto the wafer surface. This
//! module intentionally implements the simple binary-threshold model
//! (mask column open -> photoresist beneath it is exposed) rather than a
//! full aerial-image diffraction simulation. A Rayleigh-criterion-based
//! point-spread-function blur (convolving the binary mask with a
//! Gaussian/Airy kernel before thresholding) is a documented future
//! extension for modeling proximity effects near mask edges - see
//! docs/architecture.md. Get the simple case numerically right first.

use crate::grid2d::{idx, Wafer2d};
use crate::materials::Material;

/// Expose photoresist under open mask columns: every `Photoresist` voxel
/// in an exposed column becomes `PhotoresistExposed`. `mask_open` has
/// length `wafer.nx`; `true` = light passes through (exposed).
pub fn expose(wafer: &mut Wafer2d, mask_open: &[bool]) {
    assert_eq!(mask_open.len(), wafer.nx, "mask_open must have one entry per column");
    for x in 0..wafer.nx {
        if !mask_open[x] {
            continue;
        }
        for y in 0..wafer.ny {
            let i = idx(x, y, wafer.nx);
            if wafer.material[i] == Material::Photoresist {
                wafer.material[i] = Material::PhotoresistExposed;
            }
        }
    }
    wafer.process_log.push(format!(
        "expose open_columns={}",
        mask_open.iter().filter(|&&b| b).count()
    ));
}

/// Develop: strip exposed (soluble) photoresist, opening a window down to
/// whatever material sits beneath it. Modeled as converting exposed resist
/// voxels to `Void`; the etch step then acts through these openings.
pub fn develop(wafer: &mut Wafer2d) {
    let mut removed = 0usize;
    for m in wafer.material.iter_mut() {
        if *m == Material::PhotoresistExposed {
            *m = Material::Void;
            removed += 1;
        }
    }
    wafer.process_log.push(format!("develop removed_voxels={removed}"));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::materials::Dopant;

    #[test]
    fn expose_only_flips_photoresist_in_open_columns() {
        let mut w = Wafer2d::new(4, 5, 0.01, 0.01, Dopant::Boron, 1e15);
        w.spin_photoresist(3);
        let mask = [true, false, true, false];
        expose(&mut w, &mask);

        for x in 0..4 {
            for y in 0..3 {
                let i = idx(x, y, w.nx);
                let expected = if mask[x] { Material::PhotoresistExposed } else { Material::Photoresist };
                assert_eq!(w.material[i], expected, "column {x} row {y} mismatch");
            }
        }
    }

    #[test]
    fn develop_converts_exposed_resist_to_void_and_leaves_rest() {
        let mut w = Wafer2d::new(4, 5, 0.01, 0.01, Dopant::Boron, 1e15);
        w.spin_photoresist(3);
        let mask = [true, false, true, false];
        expose(&mut w, &mask);
        develop(&mut w);

        for x in 0..4 {
            for y in 0..3 {
                let i = idx(x, y, w.nx);
                let expected = if mask[x] { Material::Void } else { Material::Photoresist };
                assert_eq!(w.material[i], expected, "column {x} row {y} mismatch after develop");
            }
        }
    }
}
