//! Etch & Deposition Solver (Phase 2).
//!
//! Etch and deposition are geometric duals of each other (removing vs.
//! adding material along an advancing front), implemented here as a
//! cellular-automaton front-propagation process rather than a true
//! level-set PDE solver - CA gives visually and physically correct
//! anisotropic/isotropic behavior (including undercut) with far less
//! numerical machinery than solving a curvature-dependent Hamilton-Jacobi
//! equation, which students don't need to see to understand the concept.
//!
//! * Anisotropic etch: material is removed strictly downward (straight
//!   line from an already-open path to the surface) - models a modern
//!   directional plasma/RIE etch with no undercut.
//! * Isotropic etch: material is removed from any solid voxel adjacent to
//!   an already-open (Void) voxel, in all four directions per step - models
//!   a uniform wet chemical etch, which erodes sideways under a mask edge
//!   (undercut) exactly as real isotropic etches do.

use crate::grid2d::{idx, Wafer2d};
use crate::materials::Material;

fn is_etchable(m: Material) -> bool {
    matches!(m, Material::Silicon | Material::Oxide)
}

/// Anisotropic (directional, "straight down") etch. Only columns with an
/// existing open path from the surface (e.g. from a developed resist
/// window) advance; each step removes exactly one more voxel straight
/// down that column's current front. No lateral (sideways) erosion.
pub fn etch_anisotropic(wafer: &mut Wafer2d, depth_voxels: usize) {
    let nx = wafer.nx;
    let ny = wafer.ny;
    let mut removed = 0usize;

    for _ in 0..depth_voxels {
        for x in 0..nx {
            // Find the etch front: the first solid voxel whose immediate
            // neighbor above (or the surface itself) is already Void,
            // i.e. the etchant has a clear vertical path down to it.
            for y in 0..ny {
                let i = idx(x, y, nx);
                if wafer.material[i] == Material::Void {
                    continue; // already open - keep scanning down for the real front
                }
                // First solid voxel encountered. It can only be etched if
                // the voxel directly above it is Void (etchant has a clear
                // path down to it). Row 0 has nothing above it, so a
                // column that's still solid at row 0 is fully
                // mask-protected and can never be reached.
                let above_is_open = y > 0 && wafer.material[idx(x, y - 1, nx)] == Material::Void;
                if above_is_open && is_etchable(wafer.material[i]) {
                    wafer.material[i] = Material::Void;
                    removed += 1;
                }
                break; // whether etched or genuinely blocked, this column is done for this step
            }
        }
    }
    wafer.process_log.push(format!("etch_anisotropic depth={depth_voxels}vox removed={removed}"));
}

/// Isotropic (uniform chemical) etch. Each step, every solid voxel with
/// at least one Void 4-neighbor is removed simultaneously (computed from a
/// snapshot, so removals within a step don't cascade order-dependently).
/// This is what produces realistic undercut beneath a masking layer.
pub fn etch_isotropic(wafer: &mut Wafer2d, iterations: usize) {
    let nx = wafer.nx;
    let ny = wafer.ny;
    let mut removed = 0usize;

    for _ in 0..iterations {
        let snapshot = wafer.material.clone();
        let mut to_remove = Vec::new();
        for y in 0..ny {
            for x in 0..nx {
                let i = idx(x, y, nx);
                if !is_etchable(snapshot[i]) {
                    continue;
                }
                let neighbors_void = [
                    x > 0 && snapshot[idx(x - 1, y, nx)] == Material::Void,
                    x + 1 < nx && snapshot[idx(x + 1, y, nx)] == Material::Void,
                    y > 0 && snapshot[idx(x, y - 1, nx)] == Material::Void,
                    y + 1 < ny && snapshot[idx(x, y + 1, nx)] == Material::Void,
                ];
                if neighbors_void.iter().any(|&b| b) {
                    to_remove.push(i);
                }
            }
        }
        for i in &to_remove {
            wafer.material[*i] = Material::Void;
        }
        removed += to_remove.len();
    }
    wafer.process_log.push(format!("etch_isotropic iterations={iterations} removed={removed}"));
}

/// Deposition (CVD/PVD): the geometric inverse of etch. Blanket-fills the
/// topmost `thickness_voxels` of Void space above the existing topography
/// in every column with `material` - i.e. new material "rains down" and
/// conforms to (but does not fill beneath) whatever surface is already
/// there. Phase 2 keeps this conformal/blanket; a directional (PVD-style,
/// line-of-sight only) deposition mode is a documented future extension.
pub fn deposit(wafer: &mut Wafer2d, material: Material, thickness_voxels: usize) {
    let nx = wafer.nx;
    let ny = wafer.ny;
    let mut added = 0usize;

    for x in 0..nx {
        // Find the current surface (topmost non-Void row) in this column.
        let mut surface_row = None;
        for y in 0..ny {
            if wafer.material[idx(x, y, nx)] != Material::Void {
                surface_row = Some(y);
                break;
            }
        }
        let Some(surface_row) = surface_row else { continue };
        // Fill upward from just above the surface, thickness_voxels deep.
        let mut filled = 0;
        let mut y = surface_row;
        while filled < thickness_voxels && y > 0 {
            y -= 1;
            let i = idx(x, y, nx);
            if wafer.material[i] == Material::Void {
                wafer.material[i] = material;
                filled += 1;
                added += 1;
            } else {
                break;
            }
        }
    }
    wafer.process_log.push(format!("deposit material={material:?} thickness={thickness_voxels}vox added={added}"));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid2d::Wafer2d;
    use crate::materials::Dopant;

    /// Build a wafer with a single open (Void) column at x=2, flanked by
    /// solid Silicon everywhere else - simulates a developed resist window.
    fn wafer_with_open_window(nx: usize, ny: usize) -> Wafer2d {
        let mut w = Wafer2d::new(nx, ny, 0.01, 0.01, Dopant::Boron, 1e15);
        let open_x = nx / 2;
        for y in 0..2 {
            w.material[idx(open_x, y, nx)] = Material::Void;
        }
        w
    }

    #[test]
    fn anisotropic_etch_has_no_lateral_spread() {
        let nx = 7;
        let ny = 10;
        let mut w = wafer_with_open_window(nx, ny);
        let open_x = nx / 2;

        etch_anisotropic(&mut w, 3);

        // The open column should have etched straight down.
        for y in 0..5 {
            assert_eq!(w.material[idx(open_x, y, nx)], Material::Void, "open column row {y} should be etched");
        }
        // Every neighboring column must remain fully solid - zero undercut.
        for &x in &[open_x - 1, open_x + 1] {
            for y in 0..5 {
                assert_ne!(
                    w.material[idx(x, y, nx)],
                    Material::Void,
                    "anisotropic etch must not remove neighbor column {x} row {y}"
                );
            }
        }
    }

    #[test]
    fn isotropic_etch_undercuts_neighboring_columns() {
        let nx = 7;
        let ny = 10;
        let mut w = wafer_with_open_window(nx, ny);
        let open_x = nx / 2;

        etch_isotropic(&mut w, 3);

        // Isotropic etch should erode sideways into the columns
        // immediately adjacent to the window - i.e. produce undercut.
        let left_removed = w.material[idx(open_x - 1, 0, nx)] == Material::Void;
        let right_removed = w.material[idx(open_x + 1, 0, nx)] == Material::Void;
        assert!(left_removed || right_removed, "isotropic etch should undercut at least one neighbor");
    }

    #[test]
    fn deposit_fills_void_above_existing_surface_without_gaps() {
        let mut w = Wafer2d::new(3, 6, 0.01, 0.01, Dopant::Boron, 1e15);
        // Clear the top 3 rows to Void, simulating an etched trench.
        for x in 0..3 {
            for y in 0..3 {
                w.material[idx(x, y, w.nx)] = Material::Void;
            }
        }
        deposit(&mut w, Material::Metal, 2);

        for x in 0..3 {
            assert_eq!(w.material[idx(x, 1, w.nx)], Material::Metal);
            assert_eq!(w.material[idx(x, 2, w.nx)], Material::Metal);
            // Topmost row (0) should remain Void - only 2 voxels deposited.
            assert_eq!(w.material[idx(x, 0, w.nx)], Material::Void);
        }
    }
}
