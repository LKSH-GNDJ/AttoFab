//! Reference diffusion solver for localized silicon regions.
//!
//! Governing law: Fick's second law, discretized explicitly in 1D.
//! The solver operates only across adjacent Silicon cells; non-silicon
//! material cells block transport. Units: x in cm, D in cm^2/s, t in s,
//! concentration in ions/cm^3.

use crate::wafer::Material;

#[derive(Debug, thiserror::Error)]
pub enum DiffusionError {
    #[error("concentration/material lengths differ: {0} vs {1}")]
    LengthMismatch(usize, usize),
    #[error("dx must be finite and > 0, got {0}")]
    InvalidDx(f64),
    #[error("diffusivity must be finite and >= 0, got {0}")]
    InvalidDiffusivity(f64),
    #[error("dt must be finite and >= 0, got {0}")]
    InvalidDt(f64),
    #[error("concentration at index {0} must be finite and >= 0")]
    InvalidConcentration(usize),
    #[error("explicit scheme unstable: D*dt/dx^2 = {0} > 0.5")]
    UnstableTimestep(f64),
}

/// Advance one explicit diffusion step. Flux is allowed only between two
/// adjacent Silicon cells, so oxide/mask/air regions act as barriers.
pub fn step_1d(
    concentration_ions_cm3: &[f64],
    materials: &[Material],
    dx_cm: f64,
    diffusivity_cm2_s: f64,
    dt_s: f64,
) -> Result<Vec<f64>, DiffusionError> {
    if concentration_ions_cm3.len() != materials.len() {
        return Err(DiffusionError::LengthMismatch(
            concentration_ions_cm3.len(),
            materials.len(),
        ));
    }
    if !dx_cm.is_finite() || dx_cm <= 0.0 {
        return Err(DiffusionError::InvalidDx(dx_cm));
    }
    if !diffusivity_cm2_s.is_finite() || diffusivity_cm2_s < 0.0 {
        return Err(DiffusionError::InvalidDiffusivity(diffusivity_cm2_s));
    }
    if !dt_s.is_finite() || dt_s < 0.0 {
        return Err(DiffusionError::InvalidDt(dt_s));
    }
    for (i, &c) in concentration_ions_cm3.iter().enumerate() {
        if !c.is_finite() || c < 0.0 {
            return Err(DiffusionError::InvalidConcentration(i));
        }
    }

    let alpha = diffusivity_cm2_s * dt_s / (dx_cm * dx_cm);
    if alpha > 0.5 {
        return Err(DiffusionError::UnstableTimestep(alpha));
    }

    let mut next = concentration_ions_cm3.to_vec();
    for i in 0..concentration_ions_cm3.len() {
        if materials[i] != Material::Silicon {
            continue;
        }
        let mut laplacian = 0.0;
        if i > 0 && materials[i - 1] == Material::Silicon {
            laplacian += concentration_ions_cm3[i - 1] - concentration_ions_cm3[i];
        }
        if i + 1 < concentration_ions_cm3.len() && materials[i + 1] == Material::Silicon {
            laplacian += concentration_ions_cm3[i + 1] - concentration_ions_cm3[i];
        }
        next[i] = concentration_ions_cm3[i] + alpha * laplacian;
    }
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn one_step_stencil_matches_reference() {
        let c = [0.0, 1.0, 0.0];
        let m = [Material::Silicon; 3];
        let out = step_1d(&c, &m, 1.0, 0.2, 1.0).unwrap();
        assert!((out[0] - 0.2).abs() < 1e-12);
        assert!((out[1] - 0.6).abs() < 1e-12);
        assert!((out[2] - 0.2).abs() < 1e-12);
    }

    #[test]
    fn oxide_blocks_transport() {
        let c = [1.0, 0.0, 0.0];
        let m = [Material::Silicon, Material::SiliconDioxide, Material::Silicon];
        let out = step_1d(&c, &m, 1.0, 0.2, 1.0).unwrap();
        assert_eq!(out, c);
    }

    #[test]
    fn unstable_step_is_rejected() {
        let c = [0.0, 1.0, 0.0];
        let m = [Material::Silicon; 3];
        assert!(matches!(
            step_1d(&c, &m, 1.0, 1.0, 1.0),
            Err(DiffusionError::UnstableTimestep(_))
        ));
    }
}
