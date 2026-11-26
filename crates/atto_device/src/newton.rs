// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A coupled Newton solve for the device equations (build PRD V-01).
//!
//! Three unknowns per node — potential, electrons, holes — solved together. The
//! Jacobian is the **analytic** Scharfetter-Gummel one, not a finite-difference
//! approximation: the difference version stalls because the carrier derivatives
//! differ from the potential derivatives by twenty orders of magnitude, and a
//! single step size cannot serve both.
//!
//! Rows are scaled so every equation is O(1) before the Newton step, which is
//! what makes a line search on the residual norm meaningful.

use crate::block::{solve, Block};
use crate::mos::{EPS0, EPS_SI_R, KB, NI, QE};

#[derive(Debug, Clone, PartialEq)]
/// Why the coupled solve stopped.
pub enum NewtonError {
    BadInput(String),
    NoConvergence { iterations: usize, residual: f64 },
    Singular,
}

impl std::fmt::Display for NewtonError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NewtonError::BadInput(s) => write!(f, "{s}"),
            NewtonError::NoConvergence { iterations, residual } => {
                write!(f, "Newton did not converge: {residual:.2e} after {iterations} iterations")
            }
            NewtonError::Singular => write!(f, "the Jacobian was singular"),
        }
    }
}
impl std::error::Error for NewtonError {}

/// `B(x) = x / (exp(x) - 1)`, the Scharfetter-Gummel Bernoulli function.
fn bern(x: f64) -> f64 {
    if x.abs() < 1.0e-6 {
        1.0 - x * 0.5 + x * x / 12.0
    } else if x > 40.0 {
        x * (-x).exp()
    } else if x < -40.0 {
        -x
    } else {
        x / (x.exp() - 1.0)
    }
}

/// `B'(x)`, with the small-x expansion so it stays accurate near zero.
fn bern_p(x: f64) -> f64 {
    if x.abs() < 1.0e-4 {
        -0.5 + x / 6.0
    } else if x > 40.0 {
        (-x).exp() * (1.0 - x)
    } else if x < -40.0 {
        -1.0
    } else {
        let e = x.exp();
        ((e - 1.0) - x * e) / ((e - 1.0) * (e - 1.0))
    }
}

/// Prefactor of the `ni(T)` power law, cm^-3 K^-3/2.
///
/// The textbook value (3.87e16) gives 8.6e10 at 300 K; this is calibrated so that
/// `ni(300 K)` reproduces the manifest's 1.0e10. Two different intrinsic
/// concentrations in one tool would be worse than a calibrated constant, and the
/// calibration is declared in the registry as an assumed parameter rather than
/// hidden in the arithmetic.
pub const NI_PREFACTOR: f64 = 4.42e15;

#[derive(Debug, Clone)]
/// A one-dimensional diode the coupled solver can run.
pub struct SmallDiode {
    pub doping: Vec<f64>,
    pub h_nm: f64,
    pub temp_k: f64,
    pub mu_n: f64,
    pub mu_p: f64,
    /// SRH lifetimes, s.
    pub tau_n_s: f64,
    pub tau_p_s: f64,
}

#[derive(Debug, Clone)]
/// A converged solution: potentials, carriers and iterations.
pub struct Solution {
    pub psi: Vec<f64>,
    pub n: Vec<f64>,
    pub p: Vec<f64>,
    pub j: f64,
    pub iterations: usize,
}

impl SmallDiode {
    /// A small device for the tests.
    pub fn small() -> Self {
        let h = 50.0;
        let nn = (2000.0 / h) as usize + 1;
        let doping = (0..nn).map(|i| if (i as f64) * h < 1000.0 { 1.0e16 } else { -1.0e16 }).collect();
        SmallDiode { doping, h_nm: h, temp_k: 300.15, mu_n: 1350.0, mu_p: 480.0, tau_n_s: 1.0e-6, tau_p_s: 1.0e-6 }
    }

    /// Shockley-Read-Hall recombination at midgap and its derivatives.
    fn srh(&self, n: f64, p: f64) -> (f64, f64, f64) {
        let ni = self.ni_at();
        let d = self.tau_p_s * (n + ni) + self.tau_n_s * (p + ni);
        if !(d > 0.0) {
            return (0.0, 0.0, 0.0);
        }
        let num = n * p - ni * ni;
        let r = num / d;
        let dr_dn = (p * d - num * self.tau_p_s) / (d * d);
        let dr_dp = (n * d - num * self.tau_n_s) / (d * d);
        (r, dr_dn, dr_dp)
    }

    /// Phonon-limited mobilities at a temperature, cm^2/Vs (build PRD V-06).
    ///
    /// Lattice scattering goes as `T^-3/2` for electrons and rather steeper for
    /// holes. A fixed room-temperature constant is the thing this replaces: at
    /// 400 K the electron mobility is about a third lower.
    pub fn mobility_at(&self, temp_k: f64) -> (f64, f64) {
        if !(temp_k > 0.0) || !temp_k.is_finite() {
            return (self.mu_n, self.mu_p);
        }
        let r = 300.0 / temp_k;
        (self.mu_n * r.powf(1.5), self.mu_p * r.powf(1.8))
    }

    /// Intrinsic carrier concentration at the device temperature, cm^-3.
    ///
    /// `ni` rises by three orders of magnitude between 300 K and 400 K, which
    /// dominates every other temperature effect: a diode heated by 100 K carries
    /// far more current even though its carriers are less mobile.
    pub fn ni_at(&self) -> f64 {
        let t = self.temp_k;
        if !(t > 0.0) || !t.is_finite() {
            return NI;
        }
        // Varshni gap for silicon, then ni(T). The prefactor is calibrated so
        // that ni(300 K) reproduces the manifest's 1.0e10: the textbook 3.87e16
        // gives 8.6e10, ten times the value the rest of this model uses, and two
        // different ni's in one tool is worse than a calibrated constant.
        let eg_ev = 1.16 - 7.02e-4 * t * t / (t + 1108.0);
        let kt_ev = KB * t / QE;
        NI_PREFACTOR * t.powf(1.5) * (-eg_ev / (2.0 * kt_ev)).exp()
    }

    /// The thermal voltage.
    pub fn vt(&self) -> f64 {
        KB * self.temp_k / QE
    }

    fn n_max(&self) -> f64 {
        self.doping.iter().fold(NI, |a, b| a.max(b.abs()))
    }

    /// Contact potentials and carrier densities for a bias.
    fn contacts(&self, v: f64) -> (f64, f64, f64, f64, f64, f64) {
        let vt = self.vt();
        let ni = self.ni_at();
        let nd = self.doping.first().copied().unwrap_or(ni).abs().max(ni);
        let na = self.doping.last().copied().unwrap_or(ni).abs().max(ni);
        let nl = nd;
        let pl = ni * ni / nd;
        let pr = na;
        let nr = ni * ni / na;
        // The bias belongs to the p-side (the anode). Applying it to the n-side
        // inverts the junction: reverse bias then forward-biases it, which is
        // why a 0.2 A/cm2 "reverse" current appeared with and without SRH.
        (vt * (nd / ni).ln(), v - vt * (na / ni).ln(), nl, pl, nr, pr)
    }

    /// Residual, scaled so every row is O(1).
    fn residual(&self, x: &[[f64; 3]], v: f64) -> Vec<[f64; 3]> {
        let nn = x.len();
        let vt = self.vt();
        let h = self.h_nm * 1.0e-7;
        let eps = EPS_SI_R * EPS0;
        let nmax = self.n_max();
        let (psl, psr, nl, pl, nr, pr) = self.contacts(v);
        let (mu_n, mu_p) = self.mobility_at(self.temp_k);
        let dn = vt * mu_n / h;
        let dp = vt * mu_p / h;
        let sp = QE * nmax;
        let sc = dn * nmax;
        let mut f = vec![[0.0; 3]; nn];
        for i in 0..nn {
            if i == 0 || i == nn - 1 {
                let (ps, ns, psn) = if i == 0 { (psl, nl, pl) } else { (psr, nr, pr) };
                f[i] = [
                    (x[i][0] - ps) / vt,
                    (x[i][1] - ns) / ns,
                    (x[i][2] - psn) / psn,
                ];
                continue;
            }
            let (pm, p0, pp) = (x[i - 1][0], x[i][0], x[i + 1][0]);
            let poisson = eps * (pm - 2.0 * p0 + pp) / (h * h) + QE * (x[i][2] - x[i][1] + self.doping[i]);
            let xn = (pp - p0) / vt;
            let xb = (p0 - pm) / vt;
            let fe = dn * (x[i + 1][1] * bern(xn) - x[i][1] * bern(-xn))
                - dn * (x[i][1] * bern(xb) - x[i - 1][1] * bern(-xb));
            let yn = -xn;
            let yb = -xb;
            let fh = dp * (x[i + 1][2] * bern(yn) - x[i][2] * bern(-yn))
                - dp * (x[i][2] * bern(yb) - x[i - 1][2] * bern(-yb));
            let (r, _, _) = self.srh(x[i][1], x[i][2]);
            f[i] = [poisson / sp, (fe - r * h) / sc, (fh + r * h) / sc];
        }
        f
    }

    /// The analytic Jacobian, in 3x3 blocks.
    fn jacobian(&self, x: &[[f64; 3]], v: f64) -> (Vec<Block>, Vec<Block>, Vec<Block>) {
        let nn = x.len();
        let vt = self.vt();
        let h = self.h_nm * 1.0e-7;
        let eps = EPS_SI_R * EPS0;
        let nmax = self.n_max();
        let (_, _, nl, pl, nr, pr) = self.contacts(v);
        let (mu_n, mu_p) = self.mobility_at(self.temp_k);
        let dn = vt * mu_n / h;
        let dp = vt * mu_p / h;
        let sp = QE * nmax;
        let sc = dn * nmax;
        let mut lower = vec![[[0.0; 3]; 3]; nn];
        let mut diag = vec![[[0.0; 3]; 3]; nn];
        let mut upper = vec![[[0.0; 3]; 3]; nn];
        for i in 0..nn {
            if i == 0 || i == nn - 1 {
                let (ns, ps) = if i == 0 { (nl, pl) } else { (nr, pr) };
                diag[i] = [[1.0 / vt, 0.0, 0.0], [0.0, 1.0 / ns, 0.0], [0.0, 0.0, 1.0 / ps]];
                continue;
            }
            let (pm, p0, pp) = (x[i - 1][0], x[i][0], x[i + 1][0]);
            let k = eps / (h * h);
            // Poisson row: psi and the carriers
            lower[i][0][0] = k / sp;
            diag[i][0][0] = -2.0 * k / sp;
            upper[i][0][0] = k / sp;
            diag[i][0][1] = -QE / sp;
            diag[i][0][2] = QE / sp;

            let xn = (pp - p0) / vt;
            let xb = (p0 - pm) / vt;
            let (nm, n0, np) = (x[i - 1][1], x[i][1], x[i + 1][1]);
            // electron continuity row
            lower[i][1][1] = dn * bern(-xb) / sc;
            diag[i][1][1] = -dn * (bern(-xn) + bern(xb)) / sc;
            upper[i][1][1] = dn * bern(xn) / sc;
            diag[i][1][0] = -(dn / vt)
                * (np * bern_p(xn) + n0 * bern_p(-xn) + n0 * bern_p(xb) + nm * bern_p(-xb))
                / sc;
            upper[i][1][0] = dn * (np * bern_p(xn) + n0 * bern_p(-xn)) / vt / sc;
            lower[i][1][0] = dn * (n0 * bern_p(xb) + nm * bern_p(-xb)) / vt / sc;

            let yn = -xn;
            let yb = -xb;
            let (pmv, p0v, ppv) = (x[i - 1][2], x[i][2], x[i + 1][2]);
            // hole continuity row
            lower[i][2][2] = dp * bern(-yb) / sc;
            diag[i][2][2] = -dp * (bern(-yn) + bern(yb)) / sc;
            upper[i][2][2] = dp * bern(yn) / sc;
            diag[i][2][0] = -(dp / vt)
                * (ppv * bern_p(yn) + p0v * bern_p(-yn) + p0v * bern_p(yb) + pmv * bern_p(-yb))
                / sc
                * -1.0;
            upper[i][2][0] = -dp * (ppv * bern_p(yn) + p0v * bern_p(-yn)) / vt / sc;
            lower[i][2][0] = dp * (p0v * bern_p(yb) + pmv * bern_p(-yb)) / vt / sc * -1.0;

            // SRH couples the two carrier equations to each other
            let (_, drdn, drdp) = self.srh(x[i][1], x[i][2]);
            diag[i][1][1] -= drdn * h / sc;
            diag[i][1][2] -= drdp * h / sc;
            diag[i][2][1] += drdn * h / sc;
            diag[i][2][2] += drdp * h / sc;
        }
        (lower, diag, upper)
    }

    /// Solve at a bias, warm-starting from a previous solution.
    pub fn solve(&self, v: f64, start: Option<&Solution>) -> Result<Solution, NewtonError> {
        let nn = self.doping.len();
        if nn < 8 || !v.is_finite() || !(self.h_nm > 0.0) {
            return Err(NewtonError::BadInput("invalid device or bias".into()));
        }
        let vt = self.vt();
        let (psl, psr, nl, pl, nr, pr) = self.contacts(v);
        let mut x: Vec<[f64; 3]> = match start {
            Some(s) if s.psi.len() == nn => (0..nn).map(|i| [s.psi[i], s.n[i], s.p[i]]).collect(),
            _ => {
                let mut vv: Vec<[f64; 3]> = (0..nn)
                    .map(|i| {
                        let net = self.doping[i];
                        let a = net.abs().max(NI);
                        let phi = vt * (a / NI).ln();
                        let ps = if net >= 0.0 { phi } else { -phi };
                        [ps, a, NI * NI / a]
                    })
                    .collect();
                vv[0] = [psl, nl, pl];
                vv[nn - 1] = [psr, nr, pr];
                vv
            }
        };
        let mut residual = f64::INFINITY;
        let mut iterations = 0;
        for it in 0..80 {
            iterations = it + 1;
            let f = self.residual(&x, v);
            let norm: f64 = f.iter().flat_map(|b| b.iter()).map(|z| z * z).sum::<f64>().sqrt();
            residual = norm;
            if norm < 1e-9 {
                break;
            }
            let (lower, diag, upper) = self.jacobian(&x, v);
            let rhs: Vec<[f64; 3]> = f.iter().map(|b| [-b[0], -b[1], -b[2]]).collect();
            let Some(delta) = solve(&lower, &diag, &upper, &rhs) else {
                return Err(NewtonError::Singular);
            };
            let mut lambda = 1.0;
            let mut improved = false;
            for _ in 0..24 {
                let mut trial = x.clone();
                for i in 0..nn {
                    trial[i][0] += lambda * delta[i][0];
                    trial[i][1] = (trial[i][1] + lambda * delta[i][1]).max(1.0);
                    trial[i][2] = (trial[i][2] + lambda * delta[i][2]).max(1.0);
                }
                let ft = self.residual(&trial, v);
                let nt: f64 = ft.iter().flat_map(|b| b.iter()).map(|z| z * z).sum::<f64>().sqrt();
                if nt.is_finite() && nt < norm {
                    x = trial;
                    improved = true;
                    break;
                }
                lambda *= 0.5;
            }
            if !improved {
                break;
            }
        }
        if residual >= 1e-9 {
            return Err(NewtonError::NoConvergence { iterations, residual });
        }
        let psi: Vec<f64> = x.iter().map(|b| b[0]).collect();
        let n: Vec<f64> = x.iter().map(|b| b[1]).collect();
        let p: Vec<f64> = x.iter().map(|b| b[2]).collect();
        let j = self.current(&psi, &n, &p);
        Ok(Solution { psi, n, p, j, iterations })
    }

    /// Total current density, A/cm^2.
    ///
    /// With recombination, the electron current alone is not constant across the
    /// device — only the sum of both carrier currents is. Averaging the electron
    /// flux on its own reports the recombination profile rather than the current.
    fn current(&self, psi: &[f64], n: &[f64], p: &[f64]) -> f64 {
        let vt = self.vt();
        let h = self.h_nm * 1.0e-7;
        let (mu_n, mu_p) = self.mobility_at(self.temp_k);
        let dn = vt * mu_n / h;
        let dp = vt * mu_p / h;
        let mut sum = 0.0;
        let mut cnt = 0;
        for i in 0..psi.len() - 1 {
            let x = (psi[i + 1] - psi[i]) / vt;
            let jn = dn * (n[i + 1] * bern(x) - n[i] * bern(-x));
            let jp = dp * (p[i + 1] * bern(-x) - p[i] * bern(x));
            sum += jn - jp;
            cnt += 1;
        }
        if cnt == 0 { 0.0 } else { -QE * sum / cnt as f64 }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_bernoulli_derivative_matches_a_finite_difference() {
        for x in [-8.0, -3.0, -1.0, -0.1, 0.0, 0.1, 1.0, 3.0, 8.0] {
            let h = 1e-6;
            let fd = (bern(x + h) - bern(x - h)) / (2.0 * h);
            assert!((bern_p(x) - fd).abs() < 1e-3, "B-prime({x}) analytic {:.6} vs fd {:.6}", bern_p(x), fd);
        }
    }

    #[test]
    fn the_analytic_jacobian_matches_finite_differences() {
        // This is the test that matters: a derivative bug is invisible in the
        // solution but fatal to convergence.
        let d = SmallDiode::small();
        let nn = d.doping.len();
        let vt = d.vt();
        let (psl, psr, nl, pl, nr, pr) = d.contacts(0.2);
        let x: Vec<[f64; 3]> = (0..nn)
            .map(|i| {
                if i == 0 {
                    return [psl, nl, pl];
                }
                if i == nn - 1 {
                    return [psr, nr, pr];
                }
                let s = i as f64 / (nn - 1) as f64;
                let ps = psl + (psr - psl) * s;
                [ps, nl * ((ps - psl) / vt).exp(), pl * (-(ps - psl) / vt).exp()]
            })
            .collect();
        let (lower, diag, upper) = d.jacobian(&x, 0.2);
        let f = d.residual(&x, 0.2);
        for i in [3usize, 10, 20, nn - 5] {
            for u in 0..3 {
                let mut xp = x.clone();
                let step = match u {
                    0 => 1e-7,
                    1 => xp[i][1] * 1e-6,
                    _ => xp[i][2] * 1e-6,
                };
                xp[i][u] += step;
                let fp = d.residual(&xp, 0.2);
                for j in (i - 1)..=(i + 1) {
                    for r in 0..3 {
                        let fd = (fp[j][r] - f[j][r]) / step;
                        let an = match j as isize - i as isize {
                            -1 => lower[i][r][u],
                            0 => diag[i][r][u],
                            _ => upper[i][r][u],
                        };
                        // Entries that are tiny in both are finite-difference
                        // cancellation noise, not derivatives.
                        if an.abs() < 1e-12 && fd.abs() < 1e-12 {
                            continue;
                        }
                        let scale = an.abs().max(fd.abs());
                        assert!(
                            (an - fd).abs() / scale < 1e-2,
                            "block {i} row {r} col {u}: analytic {an:.6e} vs fd {fd:.6e}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn it_converges_at_equilibrium_and_carries_no_current() {
        let d = SmallDiode::small();
        let s = d.solve(0.0, None).expect("equilibrium converges");
        assert!(s.iterations < 40, "took {} iterations", s.iterations);
        assert!(s.j.abs() < 1e-6, "J at equilibrium {:.3e}", s.j);
    }

    #[test]
    fn a_forward_bias_carries_an_exponential_current() {
        let d = SmallDiode::small();
        let mut prev = None;
        let mut pts = Vec::new();
        for v in [0.1, 0.2, 0.3, 0.4] {
            let s = d.solve(v, prev.as_ref()).expect("forward converges");
            assert!(s.j > 0.0, "forward current must be positive at {v}: {:.3e}", s.j);
            pts.push((v, s.j));
            prev = Some(s);
        }
        for w in pts.windows(2) {
            assert!(w[1].1 > w[0].1, "J must rise: {:.3e} -> {:.3e}", w[0].1, w[1].1);
        }
        // d(ln J)/dV = 1/(n Vt); with midgap SRH the ideality sits near one
        let (a, b) = (pts[1], pts[3]);
        let ideality = 1.0 / (((b.1 / a.1).ln() / (b.0 - a.0)) * d.vt());
        assert!((0.8..1.5).contains(&ideality), "ideality {ideality:.3}");
    }

    #[test]
    fn a_reverse_bias_carries_a_small_current() {
        let d = SmallDiode::small();
        for v in [-0.2, -0.5] {
            let s = d.solve(v, None).expect("reverse converges");
            assert!(s.j.abs() < 1e-6, "reverse current at {v} is {:.3e} A/cm2", s.j);
        }
    }

    #[test]
    fn the_builtin_potential_shows_in_the_solution() {
        let d = SmallDiode::small();
        let s = d.solve(0.0, None).expect("converges");
        let drop = s.psi[0] - s.psi[s.psi.len() - 1];
        let expect = d.vt() * (1.0e32 / 1.0e20f64).ln();
        assert!((drop - expect).abs() / expect < 0.02, "drop {drop:.4} vs {expect:.4}");
    }

    #[test]
    fn a_forward_sweep_converges_and_stays_physical() {
        let d = SmallDiode::small();
        let mut prev = None;
        for v in [0.0, 0.1, 0.2, 0.3] {
            let s = d.solve(v, prev.as_ref()).expect("every bias converges");
            assert!(s.n.iter().all(|x| *x >= 1.0), "carriers stay positive at {v}");
            assert!(s.p.iter().all(|x| *x >= 1.0));
            assert!(s.psi.iter().all(|x| x.is_finite()));
            prev = Some(s);
        }
    }

    #[test]
    fn mobility_falls_with_temperature() {
        let d = SmallDiode::small();
        let (n300, p300) = d.mobility_at(300.0);
        assert!((n300 - d.mu_n).abs() < 1e-9, "300 K is the reference");
        let (n400, p400) = d.mobility_at(400.0);
        assert!(n400 < n300 && p400 < p300, "mobility must fall: {n300:.0} -> {n400:.0}");
        // T^-3/2 for electrons: (300/400)^1.5 = 0.650
        assert!((n400 / n300 - (0.75f64).powf(1.5)).abs() < 1e-9);
        // holes fall faster
        assert!(p400 / p300 < n400 / n300, "holes scatter more steeply");
        // and a bad temperature is not a number
        assert_eq!(d.mobility_at(0.0), (d.mu_n, d.mu_p));
    }

    #[test]
    fn a_hotter_diode_carries_more_current() {
        let mut cold = SmallDiode::small();
        cold.temp_k = 300.0;
        let mut hot = SmallDiode::small();
        hot.temp_k = 400.0;
        let a = cold.solve(0.3, None).expect("cold");
        let b = hot.solve(0.3, None).expect("hot");
        assert!(b.j > a.j, "a hotter junction carries more: {:.3e} vs {:.3e}", a.j, b.j);
    }

    #[test]
    fn bad_input_is_rejected() {
        let d = SmallDiode::small();
        assert!(d.solve(f64::NAN, None).is_err());
    }
}
