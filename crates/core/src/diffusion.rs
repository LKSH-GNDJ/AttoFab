//! Thermal Diffusion Solver - Fick's Second Law:
//!
//! ```text
//! dC/dt = D * d^2C/dx^2
//! ```
//!
//! solved via an explicit FTCS finite-difference scheme:
//!
//! ```text
//! C[i]^(n+1) = C[i]^n + D*dt/dx^2 * (C[i+1]^n - 2*C[i]^n + C[i-1]^n)
//! ```
//!
//! D(T) = D0 * exp(-Ea/kT) (Arrhenius, see materials.rs).
//!
//! Stability requires D*dt/dx^2 <= 0.5; the solver auto-selects a stable
//! dt (with sub-stepping) so callers just request "anneal for N minutes at
//! T degrees" without hand-tuning the numerical time step.
//!
//! Boundary conditions: zero-flux (Neumann, reflective) at both the wafer
//! surface and the deep end - standard simplifying assumption for a
//! capped/inert anneal on a wafer thick enough that the dopant front never
//! reaches the far boundary.

use crate::materials::{arrhenius, diffusion_params, Dopant};

pub struct DiffuseParams<'a> {
    pub concentration: &'a [f64],
    pub dx_um: f64,
    pub dopant: Dopant,
    pub temperature_c: f64,
    pub time_minutes: f64,
    /// Fraction of the CFL stability limit to use; must be < 0.5.
    pub safety_factor: f64,
}

pub struct DiffuseResult {
    pub concentration: Vec<f64>,
    pub d_cm2_s: f64,
    pub steps_run: u64,
    pub dt_s: f64,
}

pub fn diffusion_coefficient(dopant: Dopant, temperature_c: f64) -> f64 {
    let params = diffusion_params(dopant);
    let t_kelvin = temperature_c + 273.15;
    arrhenius(params.d0, params.ea, t_kelvin)
}

pub fn diffuse(p: &DiffuseParams) -> DiffuseResult {
    assert!(
        p.safety_factor < 0.5,
        "safety_factor must be < 0.5 for explicit FDM stability"
    );

    let d_cm2_s = diffusion_coefficient(p.dopant, p.temperature_c);
    let dx_cm = p.dx_um * 1e-4;
    let total_time_s = p.time_minutes * 60.0;

    let dt_s = p.safety_factor * dx_cm * dx_cm / d_cm2_s;
    let n_steps = ((total_time_s / dt_s).ceil() as u64).max(1);
    let actual_dt = total_time_s / n_steps as f64;

    let mut c = p.concentration.to_vec();
    let n = c.len();
    let lambda = d_cm2_s * actual_dt / (dx_cm * dx_cm);

    for _ in 0..n_steps {
        let mut c_next = vec![0.0; n];
        for i in 0..n {
            let left = if i == 0 { c[0] } else { c[i - 1] };
            let right = if i == n - 1 { c[n - 1] } else { c[i + 1] };
            c_next[i] = c[i] + lambda * (right - 2.0 * c[i] + left);
        }
        c = c_next;
    }

    DiffuseResult { concentration: c, d_cm2_s, steps_run: n_steps, dt_s: actual_dt }
}

/// Linear-interpolated crossing depth where `concentration` crosses
/// `background_concentration`. Returns None if no crossing is found.
pub fn find_junction_depth(
    concentration: &[f64],
    depth_grid_um: &[f64],
    background_concentration: f64,
) -> Option<f64> {
    for i in 0..concentration.len().saturating_sub(1) {
        let a = concentration[i] - background_concentration;
        let b = concentration[i + 1] - background_concentration;
        if (a > 0.0 && b <= 0.0) || (a < 0.0 && b >= 0.0) {
            let frac = a / (a - b);
            return Some(depth_grid_um[i] + frac * (depth_grid_um[i + 1] - depth_grid_um[i]));
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn diffusion_conserves_total_dose() {
        let n = 500;
        let dx_um = 0.01;
        let mut c = vec![0.0; n];
        // seed a narrow spike
        c[100] = 1e20;

        let before: f64 = c.iter().sum();
        let result = diffuse(&DiffuseParams {
            concentration: &c,
            dx_um,
            dopant: Dopant::Phosphorus,
            temperature_c: 1000.0,
            time_minutes: 30.0,
            safety_factor: 0.45,
        });
        let after: f64 = result.concentration.iter().sum();

        let relative_error = (after - before).abs() / before;
        assert!(relative_error < 1e-6, "mass not conserved: {}", relative_error);
    }

    #[test]
    fn diffusion_spreads_and_reduces_peak() {
        let n = 500;
        let mut c = vec![0.0; n];
        c[100] = 1e20;

        let result = diffuse(&DiffuseParams {
            concentration: &c,
            dx_um: 0.01,
            dopant: Dopant::Phosphorus,
            temperature_c: 1000.0,
            time_minutes: 30.0,
            safety_factor: 0.45,
        });

        let peak_after: f64 = result.concentration.iter().cloned().fold(0.0, f64::max);
        assert!(peak_after < 1e20);
    }

    #[test]
    fn no_negative_concentrations_appear() {
        let n = 300;
        let mut c = vec![1e15; n];
        c[50] = 1e19;

        let result = diffuse(&DiffuseParams {
            concentration: &c,
            dx_um: 0.02,
            dopant: Dopant::Boron,
            temperature_c: 1100.0,
            time_minutes: 60.0,
            safety_factor: 0.45,
        });

        assert!(result.concentration.iter().all(|&v| v >= 0.0));
    }
}
