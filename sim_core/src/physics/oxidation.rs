//! Oxidation interface-geometry accounting.
//!
//! This module deliberately does not invent oxidation kinetics. Engine.txt
//! supplies the dry/wet reactions and qualitative behavior, but does not
//! supply Deal-Grove constants or a calibrated Si-consumption/oxide-growth
//! ratio. The caller therefore provides the growth amount and expansion ratio.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OxidationGeometry {
    pub silicon_consumed_nm: f64,
    pub surface_rise_nm: f64,
    pub new_interface_z_nm: f64,
    pub new_outer_surface_z_nm: f64,
}

#[derive(Debug, thiserror::Error)]
pub enum OxidationError {
    #[error("interface position must be finite")]
    InvalidInterface,
    #[error("oxide growth must be finite and >= 0, got {0}")]
    InvalidGrowth(f64),
    #[error("expansion ratio must be finite and > 0, got {0}")]
    InvalidExpansionRatio(f64),
}

/// Convert oxide growth into silicon consumption and surface/interface motion.
/// `oxide_expansion_ratio` means oxide thickness formed per unit silicon
/// thickness consumed.
pub fn evolve_interface(
    old_interface_z_nm: f64,
    oxide_growth_nm: f64,
    oxide_expansion_ratio: f64,
) -> Result<OxidationGeometry, OxidationError> {
    if !old_interface_z_nm.is_finite() {
        return Err(OxidationError::InvalidInterface);
    }
    if !oxide_growth_nm.is_finite() || oxide_growth_nm < 0.0 {
        return Err(OxidationError::InvalidGrowth(oxide_growth_nm));
    }
    if !oxide_expansion_ratio.is_finite() || oxide_expansion_ratio <= 0.0 {
        return Err(OxidationError::InvalidExpansionRatio(oxide_expansion_ratio));
    }

    let silicon_consumed_nm = oxide_growth_nm / oxide_expansion_ratio;
    let surface_rise_nm = oxide_growth_nm - silicon_consumed_nm;
    Ok(OxidationGeometry {
        silicon_consumed_nm,
        surface_rise_nm,
        new_interface_z_nm: old_interface_z_nm - silicon_consumed_nm,
        new_outer_surface_z_nm: old_interface_z_nm + surface_rise_nm,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interface_and_surface_move_consistently() {
        let out = evolve_interface(1000.0, 100.0, 2.0).unwrap();
        assert!((out.silicon_consumed_nm - 50.0).abs() < 1e-12);
        assert!((out.surface_rise_nm - 50.0).abs() < 1e-12);
        assert!((out.new_interface_z_nm - 950.0).abs() < 1e-12);
        assert!((out.new_outer_surface_z_nm - 1050.0).abs() < 1e-12);
    }

    #[test]
    fn thickness_balance_closes() {
        let out = evolve_interface(1000.0, 120.0, 2.0).unwrap();
        let consumed = 1000.0 - out.new_interface_z_nm;
        let surface_rise = out.new_outer_surface_z_nm - 1000.0;
        assert!((consumed + surface_rise - 120.0).abs() < 1e-12);
    }
}
