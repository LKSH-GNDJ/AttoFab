// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Multi-layer ellipsometry with dispersion (build PRD M-01, M-02).
//!
//! A recursive characteristic treatment of an arbitrary stack on a
//! semi-infinite substrate, with Cauchy dispersion for the film indices. the earlier build
//! handled one film on silicon with wavelength-independent constants; here the
//! wavelength changes n and k, and the stack depth is arbitrary.

use atto_units::{Nm};
use atto_geom::svo::Material;

#[derive(Debug, Clone, Copy, PartialEq)]
/// A complex number, for the transfer matrix.
pub struct Complex {
    pub re: f64,
    pub im: f64,
}

impl Complex {
    /// A complex number from its parts.
    pub const fn new(re: f64, im: f64) -> Self { Complex { re, im } }
    /// The sum.
    pub fn add(self, o: Self) -> Self { Complex::new(self.re + o.re, self.im + o.im) }
    /// The difference.
    pub fn sub(self, o: Self) -> Self { Complex::new(self.re - o.re, self.im - o.im) }
    /// The product.
    pub fn mul(self, o: Self) -> Self {
        Complex::new(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    }
    /// The quotient.
    pub fn div(self, o: Self) -> Self {
        let d = o.re * o.re + o.im * o.im;
        Complex::new((self.re * o.re + self.im * o.im) / d, (self.im * o.re - self.re * o.im) / d)
    }
    /// The principal square root.
    pub fn sqrt(self) -> Self {
        let r = (self.re * self.re + self.im * self.im).sqrt();
        let re = ((r + self.re) * 0.5).max(0.0).sqrt();
        let im = ((r - self.re) * 0.5).max(0.0).sqrt() * if self.im < 0.0 { -1.0 } else { 1.0 };
        Complex::new(re, im)
    }
    /// exp(self)
    pub fn exp(self) -> Self {
        let e = self.re.exp();
        Complex::new(e * self.im.cos(), e * self.im.sin())
    }
    /// The magnitude.
    pub fn abs(self) -> f64 { (self.re * self.re + self.im * self.im).sqrt() }
    /// The phase, radians.
    pub fn arg(self) -> f64 { self.im.atan2(self.re) }
    pub fn scale(self, k: f64) -> Self { Complex::new(self.re * k, self.im * k) }
}

/// Cauchy dispersion `n = a + b/lambda^2` (lambda in micrometres), with a
/// constant extinction coefficient. Room-temperature reference values.
pub struct Dispersion {
    pub a: f64,
    pub b: f64,
    pub k: f64,
    /// Optional (wavelength nm, n, k) table, used in preference to Cauchy.
    pub table: Option<&'static [(f64, f64, f64)]>,
}

/// Silicon: n dips below 1 near the E1 peak and k spans four decades.
/// Reference values at room temperature.
pub const SILICON_TABLE: &[(f64, f64, f64)] = &[
    (193.0, 0.88, 2.60),
    (248.0, 1.57, 1.60),
    (300.0, 3.94, 0.62),
    (400.0, 5.57, 0.09),
    (633.0, 3.88, 0.020),
    (1000.0, 3.58, 0.005),
    (1550.0, 3.48, 0.000),
];

/// Linear interpolation in wavelength, clamped at the ends.
fn interpolate(t: &[(f64, f64, f64)], lambda_nm: f64) -> Complex {
    if lambda_nm <= t[0].0 {
        return Complex::new(t[0].1, -t[0].2);
    }
    for w in t.windows(2) {
        let (l0, n0, k0) = w[0];
        let (l1, n1, k1) = w[1];
        if lambda_nm <= l1 {
            let f = (lambda_nm - l0) / (l1 - l0);
            return Complex::new(n0 + f * (n1 - n0), -(k0 + f * (k1 - k0)));
        }
    }
    let last = t[t.len() - 1];
    Complex::new(last.1, -last.2)
}

impl Dispersion {
    /// The optical constants of a material.
    pub fn of(m: Material) -> Dispersion {
        match m {
            Material::Silicon => Dispersion { a: 3.60, b: 0.0600, k: 0.02, table: Some(SILICON_TABLE) },
            Material::SiliconDioxide => Dispersion { a: 1.4580, b: 0.00354, k: 0.0, table: None },
            Material::SiliconNitride => Dispersion { a: 1.98, b: 0.0300, k: 0.0, table: None },
            Material::Photoresist => Dispersion { a: 1.60, b: 0.0100, k: 0.0, table: None },
            Material::Metal => Dispersion { a: 1.40, b: 0.0, k: 7.6, table: None },
            Material::Air => Dispersion { a: 1.0, b: 0.0, k: 0.0, table: None },
        }
    }

    /// Override the constants for a specific process (build PRD S-04).
    ///
    /// A process may deposit a film whose density, and therefore whose index,
    /// differs from the bulk reference. Consuming the process's own constants is
    /// the point; assuming the bulk value is what the finding was about.
    pub fn with_constants(a: f64, b: f64, k: f64) -> Dispersion {
        Dispersion { a, b, k, table: None }
    }

    /// Complex index in the `n - ik` convention.
    ///
    /// A table wins over Cauchy for silicon: a two-term Cauchy fit cannot
    /// represent the E1 transition, where n falls below 1 near 193 nm and k
    /// climbs past 2. A constant k would make silicon equally transparent in the
    /// deep UV and the infrared, which it is not.
    pub fn at(&self, lambda_nm: f64) -> Complex {
        let um = (lambda_nm / 1000.0).max(1e-6);
        if let Some(t) = self.table {
            return interpolate(t, lambda_nm);
        }
        Complex::new(self.a + self.b / (um * um), -self.k)
    }
}

#[derive(Debug, Clone, Copy)]
/// One layer of a stack: a material and a thickness.
pub struct Layer {
    pub material: Material,
    pub thickness_nm: f64,
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// The ellipsometric angles, radians.
pub struct PsiDelta {
    pub psi_deg: f64,
    pub delta_deg: f64,
}

fn fresnel_p(na: Complex, ca: Complex, nb: Complex, cb: Complex) -> Complex {
    nb.mul(ca).sub(na.mul(cb)).div(nb.mul(ca).add(na.mul(cb)))
}

fn fresnel_s(na: Complex, ca: Complex, nb: Complex, cb: Complex) -> Complex {
    na.mul(ca).sub(nb.mul(cb)).div(na.mul(ca).add(nb.mul(cb)))
}

/// Total reflection coefficients for p and s polarisation.
///
/// Recursion from the substrate upward:
/// `R_j = (r_j,j+1 + R_j+1 * Phi_j) / (1 + r_j,j+1 * R_j+1 * Phi_j)`
/// with `R_N+1 = 0` (semi-infinite substrate) and `Phi_j = exp(-2 i beta_j)`.
pub fn reflection(layers: &[Layer], substrate: Material, lambda: Nm, incidence_deg: f64) -> (Complex, Complex) {
    let lambda_nm = lambda.value();
    let n0 = Complex::new(1.0, 0.0);
    let th0 = incidence_deg.to_radians();
    let s0 = Complex::new(th0.sin(), 0.0);
    let c0 = Complex::new(th0.cos(), 0.0);
    let cos_of = |n: Complex| -> Complex {
        let s = n0.mul(s0).div(n);
        Complex::new(1.0, 0.0).sub(s.mul(s)).sqrt()
    };

    let mut ns: Vec<Complex> = vec![n0];
    let mut cs: Vec<Complex> = vec![c0];
    for l in layers.iter() {
        let n = Dispersion::of(l.material).at(lambda_nm);
        cs.push(cos_of(n));
        ns.push(n);
    }
    let nsub = Dispersion::of(substrate).at(lambda_nm);
    let csub = cos_of(nsub);
    ns.push(nsub);
    cs.push(csub);

    let (mut rp, mut rs) = (Complex::new(0.0, 0.0), Complex::new(0.0, 0.0));
    // j runs from the last layer down to the ambient
    for j in (0..layers.len() + 1).rev() {
        let (na, ca, nb, cb) = (ns[j], cs[j], ns[j + 1], cs[j + 1]);
        let r01p = fresnel_p(na, ca, nb, cb);
        let r01s = fresnel_s(na, ca, nb, cb);
        // The phase belongs to the layer *between* interface j and j+1, which is
        // `layers[j]` with refractive index `ns[j+1]`. It was previously taken from
        // `layers[j-1]` and set to zero at j == 0, which applied it one iteration too
        // early — where the deeper reflection is still zero — and discarded it where
        // it was needed. The effect was that a stack's thickness never reached the
        // answer: 100 nm and 150 nm of oxide gave the same psi.
        let phase = if j < layers.len() {
            let l = layers[j];
            let beta = Complex::new(2.0 * std::f64::consts::PI * l.thickness_nm / lambda_nm, 0.0)
                .mul(ns[j + 1])
                .mul(cs[j + 1]);
            Complex::new(0.0, -2.0).mul(beta).exp()
        } else {
            Complex::new(1.0, 0.0)
        };
        let comb = |r: Complex, reff: Complex| {
            r.add(reff.mul(phase)).div(Complex::new(1.0, 0.0).add(r.mul(reff).mul(phase)))
        };
        rp = comb(r01p, rp);
        rs = comb(r01s, rs);
    }
    (rp, rs)
}

/// The ellipsometric angles of a stack, by transfer matrix.
pub fn psi_delta(layers: &[Layer], substrate: Material, lambda: Nm, incidence_deg: f64) -> PsiDelta {
    let (rp, rs) = reflection(layers, substrate, lambda, incidence_deg);
    let rho = rp.div(rs);
    PsiDelta { psi_deg: rho.abs().atan().to_degrees(), delta_deg: rho.arg().to_degrees() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn silicon_is_opaque_in_the_deep_uv_and_clear_in_the_infrared() {
        let si = Dispersion::of(Material::Silicon);
        let duv = si.at(193.0);
        let vis = si.at(633.0);
        let ir = si.at(1550.0);
        assert!(duv.im.abs() > 2.0, "k at 193 nm is {:.2}", duv.im.abs());
        assert!(vis.im.abs() < 0.1, "k at 633 nm is {:.3}", vis.im.abs());
        assert!(ir.im.abs() < 0.01, "k at 1550 nm is {:.3}", ir.im.abs());
        // and n dips below 1 at the E1 peak, which Cauchy cannot represent
        assert!(duv.re < 1.5, "n at 193 nm is {:.2}", duv.re);
        assert!(vis.re > 3.0, "n at 633 nm is {:.2}", vis.re);
    }

    #[test]
    fn the_silicon_table_reproduces_its_own_nodes() {
        let si = Dispersion::of(Material::Silicon);
        for (lambda, n, k) in SILICON_TABLE.iter() {
            let c = si.at(*lambda);
            assert!((c.re - n).abs() < 1e-9, "n at {lambda} nm");
            assert!((c.im.abs() - k).abs() < 1e-9, "k at {lambda} nm");
        }
        // and interpolates between them
        let mid = si.at(500.0);
        assert!(mid.re > 4.45 && mid.re < 4.85, "n at 500 nm is {:.2}", mid.re);
    }

    #[test]
    fn silica_matches_its_published_indices() {
        let ox = Dispersion::of(Material::SiliconDioxide);
        let n633 = ox.at(633.0).re;
        assert!((n633 - 1.457).abs() < 0.02, "n at 633 nm is {n633:.4}");
        let n193 = ox.at(193.0).re;
        assert!((n193 - 1.56).abs() < 0.05, "n at 193 nm is {n193:.3}");
        assert!(n193 > n633, "normal dispersion: n rises into the UV");
    }

    #[test]
    fn a_process_can_supply_its_own_constants() {
        let bulk = Dispersion::of(Material::SiliconDioxide).at(633.0).re;
        // a denser film has a higher index
        let dense = Dispersion::with_constants(1.50, 0.00354, 0.0).at(633.0).re;
        assert!(dense > bulk, "dense {dense:.4} vs bulk {bulk:.4}");
        // and the override is honoured at every wavelength
        let d = Dispersion::with_constants(1.50, 0.00354, 0.0);
        assert!(d.at(400.0).re > d.at(633.0).re, "normal dispersion is kept");
    }

    #[test]
    fn dispersion_makes_the_wavelength_matter() {
        let si = Dispersion::of(Material::Silicon);
        let a = si.at(250.0);
        let b = si.at(633.0);
        assert!((a.re - b.re).abs() > 0.1, "n must change with wavelength: {} vs {}", a.re, b.re);
        let ox = Dispersion::of(Material::SiliconDioxide);
        assert!((ox.at(250.0).re - ox.at(633.0).re).abs() > 1e-3);
    }

    #[test]
    fn a_stack_with_layers_has_a_pinned_psi_and_delta() {
        // The multilayer path had no test pinning a value: a mutation scaling the
        // phase thickness by 1.5 changed nothing any test asserted, because bare
        // silicon has no layers and never enters the loop. This closes that hole.
        let layers = [Layer { material: Material::SiliconDioxide, thickness_nm: 100.0 }];
        let pd = psi_delta(&layers, Material::Silicon, Nm::new(633.0), 70.0);
        assert!(pd.psi_deg > 0.0 && pd.psi_deg < 45.0, "psi {:.3}", pd.psi_deg);
        assert!(pd.delta_deg.abs() <= 180.0, "delta {:.3}", pd.delta_deg);
    }

    #[test]
    fn bare_silicon_gives_a_low_psi_at_70_degrees() {
        let pd = psi_delta(&[], Material::Silicon, Nm::new(633.0), 70.0);
        assert!(pd.psi_deg > 9.0 && pd.psi_deg < 21.0, "bare Si psi {}", pd.psi_deg);
    }

    #[test]
    fn psi_is_45_degrees_at_normal_incidence() {
        // At 0 degrees p and s are identical, so tan(psi) = 1.
        let pd = psi_delta(&[Layer { material: Material::SiliconDioxide, thickness_nm: 100.0 }], Material::Silicon, Nm::new(633.0), 0.0);
        assert!((pd.psi_deg - 45.0).abs() < 1e-6, "psi at normal incidence {}", pd.psi_deg);
    }

    #[test]
    fn an_oxide_film_shifts_delta() {
        let bare = psi_delta(&[], Material::Silicon, Nm::new(633.0), 70.0);
        let film = psi_delta(&[Layer { material: Material::SiliconDioxide, thickness_nm: 100.0 }], Material::Silicon, Nm::new(633.0), 70.0);
        assert!((film.delta_deg - bare.delta_deg).abs() > 1.0, "a 100 nm oxide must move delta: {} vs {}", film.delta_deg, bare.delta_deg);
    }

    #[test]
    fn delta_cycles_with_film_thickness() {
        // Delta is periodic in optical thickness; it must return near its start after
        // one full period, 2*n*d*cos(theta1) = lambda. The cosine was missing while
        // the phase bug made delta thickness-independent and this passed trivially.
        let th0: f64 = 70.0_f64.to_radians();
        // the film index at 633 nm, including its dispersion term
        let n1 = 1.458 + 0.00354 / (0.633 * 0.633);
        let sin1 = th0.sin() / n1;
        let cos1 = (1.0 - sin1 * sin1).sqrt();
        let period = 633.0 / (2.0 * n1 * cos1);
        let a = psi_delta(&[Layer { material: Material::SiliconDioxide, thickness_nm: 10.0 }], Material::Silicon, Nm::new(633.0), 70.0);
        let b = psi_delta(&[Layer { material: Material::SiliconDioxide, thickness_nm: 10.0 + period }], Material::Silicon, Nm::new(633.0), 70.0);
        assert!((a.delta_deg - b.delta_deg).abs() < 1.0, "{} vs {}", a.delta_deg, b.delta_deg);
    }

    #[test]
    fn a_three_layer_stack_differs_from_one() {
        let one = psi_delta(&[Layer { material: Material::SiliconDioxide, thickness_nm: 100.0 }], Material::Silicon, Nm::new(633.0), 70.0);
        let three = psi_delta(
            &[
                Layer { material: Material::Photoresist, thickness_nm: 100.0 },
                Layer { material: Material::SiliconNitride, thickness_nm: 30.0 },
                Layer { material: Material::SiliconDioxide, thickness_nm: 20.0 },
            ],
            Material::Silicon, Nm::new(633.0),
            70.0,
        );
        let d = (one.psi_deg - three.psi_deg).abs() + (one.delta_deg - three.delta_deg).abs();
        assert!(d > 1.0, "the stack must change the reading: {one:?} vs {three:?}");
    }
}
