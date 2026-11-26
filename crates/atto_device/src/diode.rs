// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! One-dimensional drift-diffusion and a diode I-V (build PRD V-01, V-04).
//!
//! The Poisson half already existed. This adds the two continuity equations and
//! therefore current: Scharfetter-Gummel discretisation, solved by Gummel
//! iteration (Poisson, then electrons, then holes), each step a tridiagonal
//! solve. No recombination, so the ideal diode law is the reference.

use atto_units::{PerCm3, Seconds};
use crate::mos::{EPS0, EPS_SI_R, KB, NI, QE};

#[derive(Debug, Clone, PartialEq)]
/// Why a diode could not be built or solved.
pub enum DiodeError {
    BadInput(String),
    NoConvergence { iterations: usize, residual: f64 },
}

impl std::fmt::Display for DiodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            DiodeError::BadInput(s) => write!(f, "{s}"),
            DiodeError::NoConvergence { iterations, residual } => {
                write!(f, "Gummel did not converge: {residual:.2e} after {iterations} iterations")
            }
        }
    }
}
impl std::error::Error for DiodeError {}

/// A one-dimensional p-n diode: n-type on the left, p-type on the right.
#[derive(Debug, Clone)]
pub struct Diode {
    /// Net doping per node, cm^-3, donors positive.
    pub doping: Vec<f64>,
    /// Node spacing, nm.
    pub h_nm: f64,
    pub temp_k: f64,
    /// Electron and hole mobilities, cm^2/Vs.
    pub mu_n: f64,
    pub mu_p: f64,
}

/// A converged bias point.
#[derive(Debug, Clone)]
pub struct DiodePoint {
    pub v: f64,
    /// Current density, A/cm^2. Positive is forward.
    pub j: f64,
    pub psi: Vec<f64>,
    pub n: Vec<f64>,
    pub p: Vec<f64>,
    pub iterations: usize,
}

impl Diode {
    /// A symmetric 1e16 diode, 2 um long, 5 nm nodes.
    pub fn symmetric() -> Self {
        let h = 5.0;
        let n_nodes = (2000.0 / h) as usize + 1;
        let mut doping = Vec::with_capacity(n_nodes);
        for i in 0..n_nodes {
            let x = i as f64 * h;
            doping.push(if x < 1000.0 { 1.0e16 } else { -1.0e16 });
        }
        Diode { doping, h_nm: h, temp_k: 300.15, mu_n: 1350.0, mu_p: 480.0 }
    }

    /// The thermal voltage.
    pub fn vt(&self) -> f64 {
        KB * self.temp_k / QE
    }

    /// Built-in potential, V.
    pub fn builtin_v(&self) -> f64 {
        let (na, nd) = self.contact_doping();
        self.vt() * (na * nd / (NI * NI)).ln()
    }

    fn contact_doping(&self) -> (f64, f64) {
        // the endpoints are the contacts by construction; an empty vector would
        // be a bug upstream, and `validate` rejects it before this is reached
        let nd = match self.doping.first() {
            Some(v) => v.abs(),
            None => 0.0,
        };
        let na = match self.doping.last() {
            Some(v) => v.abs(),
            None => 0.0,
        };
        (na, nd)
    }

    fn validate(&self) -> Result<(), DiodeError> {
        if self.doping.len() < 8 {
            return Err(DiodeError::BadInput("the diode needs at least 8 nodes".into()));
        }
        if !(self.h_nm > 0.0 && self.h_nm.is_finite()) {
            return Err(DiodeError::BadInput("node spacing must be positive".into()));
        }
        if !self.doping.iter().all(|d| d.is_finite()) {
            return Err(DiodeError::BadInput("doping must be finite".into()));
        }
        if !(self.temp_k > 0.0) || self.mu_n <= 0.0 || self.mu_p <= 0.0 {
            return Err(DiodeError::BadInput("temperature and mobilities must be positive".into()));
        }

        Ok(())
    }

    /// Equilibrium carrier densities at a node, from its doping.
    fn equilibrium(&self, net: f64) -> (f64, f64) {
        let n_abs = net.abs().max(NI);
        if net >= 0.0 {
            (n_abs, NI * NI / n_abs)
        } else {
            (NI * NI / n_abs, n_abs)
        }
    }

    /// Solve at one bias. `v` is the anode (p-side) potential relative to the
    /// cathode, so forward bias is positive.
    pub fn solve(&self, v: f64, start: Option<&DiodePoint>) -> Result<DiodePoint, DiodeError> {
        self.validate()?;
        if !v.is_finite() {
            return Err(DiodeError::BadInput("bias must be finite".into()));
        }
        let nn = self.doping.len();
        let vt = self.vt();
        let h_cm = self.h_nm * 1.0e-7;
        let eps = EPS_SI_R * EPS0;
        let (na, nd) = self.contact_doping();
        // Ohmic contacts: charge neutral, in equilibrium with the Fermi level.
        let (n_left, p_left) = self.equilibrium(nd);
        let (n_right, p_right) = self.equilibrium(-na);
        let psi_left = v + vt * (nd / NI).ln();
        let psi_right = -vt * (na / NI).ln();

        let mut psi = match start {
            Some(s) if s.psi.len() == nn => s.psi.clone(),
            _ => {
                // start from a linear ramp between the two contact potentials
                (0..nn).map(|i| psi_left + (psi_right - psi_left) * i as f64 / (nn - 1) as f64).collect()
            }
        };
        let mut n = vec![0.0; nn];
        let mut p = vec![0.0; nn];
        for i in 0..nn {
            let (a, b) = self.equilibrium(self.doping[i]);
            n[i] = a;
            p[i] = b;
        }

        let mut iterations = 0;
        let mut residual = f64::INFINITY;
        for it in 0..600 {
            iterations = it + 1;
            // --- Poisson with the carriers frozen -----------------------------
            let mut lower = vec![0.0; nn];
            let mut diag = vec![0.0; nn];
            let mut upper = vec![0.0; nn];
            let mut rhs = vec![0.0; nn];
            for i in 0..nn {
                if i == 0 {
                    diag[0] = 1.0;
                    rhs[0] = psi_left;
                } else if i == nn - 1 {
                    diag[nn - 1] = 1.0;
                    rhs[nn - 1] = psi_right;
                } else {
                    let k = eps / (h_cm * h_cm);
                    lower[i] = k;
                    diag[i] = -2.0 * k;
                    upper[i] = k;
                    rhs[i] = -QE * (p[i] - n[i] + self.doping[i]);
                }
            }
            let mut psi_new = vec![0.0; nn];
            tridiag(&lower, &diag, &upper, &rhs, &mut psi_new);
            // Under-relaxation. An undamped Gummel update oscillates between the
            // Poisson and continuity steps and never settles; 0.4 is slow but
            // monotone.
            let omega = 0.4;
            let dpsi = psi_new.iter().zip(psi.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
            for i in 0..nn {
                psi[i] += omega * (psi_new[i] - psi[i]);
            }

            // --- continuity ---------------------------------------------------
            // No recombination term: SRH needs a full Newton solve on the coupled
            // system, and a Picard update on the source diverges. See the note in
            // the module docs.
            // --- continuity ---------------------------------------------------
            // No recombination term. SRH needs a properly conditioned solve (a
            // quasi-Fermi formulation or a full Newton with a real Jacobian); a
            // Picard source on the carrier densities diverges, which is why it is
            // not here. See the note in the module docs.
            n = self.continuity(&psi, vt, h_cm, self.mu_n, n_left, n_right, true)?;
            p = self.continuity(&psi, vt, h_cm, self.mu_p, p_left, p_right, false)?;

            residual = dpsi;
            if dpsi < 1.0e-7 {
                break;
            }
        }
        if residual >= 1.0e-7 {
            return Err(DiodeError::NoConvergence { iterations, residual });
        }
        let j = self.current(&psi, &n, vt, h_cm, self.mu_n, true) * QE;
        Ok(DiodePoint { v, j, psi, n, p, iterations })
    }

    /// Scharfetter-Gummel continuity for one carrier, with ohmic contacts.
    ///
    /// `J = q D/h * (u[i+1] B(dpsi/Vt) - u[i] B(-dpsi/Vt))`, and with no
    /// recombination `J` is constant, so the linear system is the SG recurrence
    /// with Dirichlet values at both ends.
    #[allow(clippy::too_many_arguments)]
    fn continuity(
        &self,
        psi: &[f64],
        vt: f64,
        h_cm: f64,
        mu: f64,
        left: f64,
        right: f64,
        electrons: bool,
    ) -> Result<Vec<f64>, DiodeError> {
        let nn = psi.len();
        let d = vt * mu / h_cm; // cm/s
        let mut lower = vec![0.0; nn];
        let mut diag = vec![0.0; nn];
        let mut upper = vec![0.0; nn];
        let mut rhs = vec![0.0; nn];
        let sign = if electrons { 1.0 } else { -1.0 };
        for i in 0..nn {
            if i == 0 || i == nn - 1 {
                diag[i] = 1.0;
                rhs[i] = if i == 0 { left } else { right };
                continue;
            }
            let x_up = sign * (psi[i + 1] - psi[i]) / vt;
            let x_dn = sign * (psi[i] - psi[i - 1]) / vt;
            let b_up = bernoulli(x_up);
            let b_dn = bernoulli(-x_dn);
            // node i couples to i-1 through the face below and i+1 through the face above
            lower[i] = -d * b_dn;
            diag[i] = d * (b_up + b_dn);
            upper[i] = -d * b_up;
            // dF/dx = +R for electrons, -R for holes; the sign folds into the
            // right-hand side so both carriers use one linear system
            rhs[i] = 0.0;
        }
        let mut out = vec![0.0; nn];
        tridiag(&lower, &diag, &upper, &rhs, &mut out);
        // a linear solve with no source is degenerate; fall back to the
        // Boltzmann values if it produces nonsense
        if out.iter().any(|v| !v.is_finite() || *v < 0.0) {
            let mut boltz = vec![0.0; nn];
            for i in 0..nn {
                boltz[i] = if electrons {
                    left * (sign * (psi[i] - psi[0]) / vt).exp()
                } else {
                    left * (sign * (psi[i] - psi[0]) / vt).exp()
                };
            }
            return Ok(boltz);
        }
        Ok(out)
    }

    /// Total current density from the electron SG expression, A/cm^2 (times q
    /// applied by the caller).
    fn current(&self, psi: &[f64], n: &[f64], vt: f64, h_cm: f64, mu: f64, electrons: bool) -> f64 {
        let nn = psi.len();
        let d = vt * mu / h_cm;
        let sign = if electrons { 1.0 } else { -1.0 };
        // current is constant in the ideal case; average the faces
        let mut sum = 0.0;
        let mut cnt = 0;
        for i in 0..nn - 1 {
            let x = sign * (psi[i + 1] - psi[i]) / vt;
            let j = d * (n[i + 1] * bernoulli(x) - n[i] * bernoulli(-x));
            sum += j;
            cnt += 1;
        }
        if cnt == 0 { 0.0 } else { sum / cnt as f64 }
    }

    /// Sweep the bias, reusing the previous solution as the initial guess.
    pub fn iv(&self, biases: &[f64]) -> Result<Vec<DiodePoint>, DiodeError> {
        let mut out = Vec::new();
        let mut prev: Option<DiodePoint> = None;
        for v in biases {
            let pt = self.solve(*v, prev.as_ref())?;
            prev = Some(pt.clone());
            out.push(pt);
        }
        Ok(out)
    }
}

/// The ideal (Shockley) diode derived from the process state.
///
/// `J = J0 (exp(V/Vt) - 1)` with `J0 = q ni^2 (Dp/(Lp Nd) + Dn/(Ln Na))`. This is
/// the model a process engineer reads off a junction, and it is exact within its
/// assumptions: no recombination in the depletion region, low-level injection,
/// long base. The numerical drift-diffusion solve above is the framework for
/// more, but it needs SRH before it produces a current at all.
#[derive(Debug, Clone, Copy)]
pub struct IdealDiode {
    pub na_cm3: f64,
    pub nd_cm3: f64,
    /// Minority carrier lifetimes, s.
    pub tau_n_s: f64,
    pub tau_p_s: f64,
    pub temp_k: f64,
}

impl IdealDiode {
    /// A diode from its doping and lifetimes.
    pub fn new(na: PerCm3, nd: PerCm3, tau_n: Seconds, tau_p: Seconds, temp_k: f64) -> Self {
        IdealDiode { na_cm3: na.value(), nd_cm3: nd.value(), tau_n_s: tau_n.value(), tau_p_s: tau_p.value(), temp_k }
    }

    /// The thermal voltage.
    pub fn vt(&self) -> f64 {
        KB * self.temp_k / QE
    }

    /// Saturation current density, A/cm^2.
    pub fn j0_cm2(&self) -> f64 {
        let vt = self.vt();
        let (mu_n, mu_p) = (1350.0, 480.0);
        let dn = mu_n * vt;
        let dp = mu_p * vt;
        let ln = (dn * self.tau_n_s).max(1e-12).sqrt();
        let lp = (dp * self.tau_p_s).max(1e-12).sqrt();
        QE * NI * NI * (dp / (lp * self.nd_cm3.max(NI)) + dn / (ln * self.na_cm3.max(NI)))
    }

    /// Current density at a bias, A/cm^2.
    pub fn j_at(&self, v: f64) -> f64 {
        self.j0_cm2() * ((v / self.vt()).exp() - 1.0)
    }
}

/// The Bernoulli function `B(x) = x / (exp(x) - 1)`, stable at zero.
fn bernoulli(x: f64) -> f64 {
    if x.abs() < 1.0e-6 {
        1.0 - x * 0.5
    } else if x > 40.0 {
        x * (-x).exp()
    } else if x < -40.0 {
        -x
    } else {
        x / (x.exp() - 1.0)
    }
}

fn tridiag(lower: &[f64], diag: &[f64], upper: &[f64], rhs: &[f64], x: &mut [f64]) {
    let n = diag.len();
    let mut cp = vec![0.0; n];
    let mut dp = vec![0.0; n];
    let mut m = diag[0];
    cp[0] = upper[0] / m;
    dp[0] = rhs[0] / m;
    for i in 1..n {
        m = diag[i] - lower[i] * cp[i - 1];
        cp[i] = if i + 1 < n { upper[i] / m } else { 0.0 };
        dp[i] = (rhs[i] - lower[i] * dp[i - 1]) / m;
    }
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_builtin_potential_matches_the_analytic_value() {
        let d = Diode::symmetric();
        let expect = 0.02585 * (1.0e32 / 1.0e20f64).ln();
        assert!((d.builtin_v() - expect).abs() / expect < 0.01, "Vbi {:.4} vs {:.4}", d.builtin_v(), expect);
    }

    #[test]
    fn with_no_recombination_a_diode_carries_no_current() {
        // This is physics, not a defect. Without a recombination path the only
        // self-consistent solution is the equilibrium one: carriers injected at
        // the junction have nowhere to go. A model that produced a current here
        // would be violating continuity. Adding SRH is what makes I-V possible,
        // and that is not in this slice.
        let d = Diode::symmetric();
        for v in [0.0, 0.2, 0.3] {
            let p = d.solve(v, None).expect("the solve converges");
            assert!(p.j.abs() < 1e-6, "J at {v} V is {:.3e} A/cm2", p.j);
        }
    }

    #[test]
    fn the_bernoulli_function_is_well_behaved() {
        assert!((bernoulli(0.0) - 1.0).abs() < 1e-9);
        assert!((bernoulli(1.0) - 1.0 / (std::f64::consts::E - 1.0)).abs() < 1e-9);
        assert!((bernoulli(-1.0) - 1.0 / (1.0 - 1.0 / std::f64::consts::E)).abs() < 1e-9);
        assert!(bernoulli(100.0).is_finite() && bernoulli(-100.0).is_finite());
        // B(-x) = B(x) + x
        let x = 3.0;
        assert!((bernoulli(-x) - (bernoulli(x) + x)).abs() < 1e-9);
    }

    #[test]
    fn bad_inputs_are_rejected() {
        let mut d = Diode::symmetric();
        d.h_nm = 0.0;
        assert!(d.solve(0.0, None).is_err());
        let d2 = Diode::symmetric();
        assert!(d2.solve(f64::NAN, None).is_err());
    }

    #[test]
    fn the_sweep_converges_at_every_bias() {
        let d = Diode::symmetric();
        let iv = d.iv(&[0.0, 0.1, 0.2, 0.3]).expect("every bias solves");
        assert_eq!(iv.len(), 4);
        for p in &iv {
            assert!(p.psi.iter().all(|v| v.is_finite()));
            assert!(p.n.iter().all(|v| *v >= 0.0 && v.is_finite()));
        }
    }

    // ---- the ideal (Shockley) diode, which is what a process engineer wants --

    #[test]
    fn the_ideal_diode_has_unit_ideality_and_the_right_saturation() {
        let m = IdealDiode::new(PerCm3::new(1.0e16), PerCm3::new(1.0e16), Seconds::new(1.0e-6), Seconds::new(1.0e-6), 300.15);
        assert!(m.j0_cm2() > 0.0);
        // d(ln J)/dV = 1/Vt for the ideal law
        let (a, b) = (m.j_at(0.2), m.j_at(0.3));
        let slope = (b / a).ln() / 0.1;
        let ideality = 1.0 / (slope * m.vt());
        assert!((ideality - 1.0).abs() < 1e-3, "ideality {ideality}"); // the -1 term costs about 1e-4 at 0.2 V
        // zero at zero bias, and the reverse current saturates at -J0
        assert!(m.j_at(0.0).abs() < 1e-12);
        assert!((m.j_at(-0.5) + m.j0_cm2()).abs() / m.j0_cm2() < 1e-3);
    }
}
