// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Three-dimensional field solve (PRD V-07).
//!
//! # Scope
//!
//! V-07 asks for three-dimensional device simulation "once the geometry representation
//! supports it". The geometry is a sparse voxel octree, which is three-dimensional, so
//! the condition the PRD names is met. What that work needs first is a **three-dimensional
//! Poisson solve**, and that is what this module provides: a uniform-grid solver with
//! Dirichlet boundaries, verified against a function whose Laplacian is known exactly.
//!
//! It is a first slice, not the whole of V-07: it solves for a potential from a charge
//! density, which is the core of the device problem. It does not carry carrier continuity
//! in three dimensions — that is stated here rather than left for a reader to discover.

/// A uniform grid over a box, with a value at each node.
#[derive(Debug, Clone, PartialEq)]
pub struct Volume3D {
    /// Nodes per axis.
    pub n: usize,
    /// Node spacing, nm.
    pub h_nm: f64,
    /// Potential at each node, indexed `z * n * n + y * n + x`.
    pub phi: Vec<f64>,
}

impl Volume3D {
    /// An empty volume of `n` nodes per axis at spacing `h`.
    pub fn new(n: usize, h_nm: f64) -> Option<Self> {
        if n < 3 || !(h_nm > 0.0) || !h_nm.is_finite() {
            return None;
        }
        Some(Volume3D { n, h_nm, phi: vec![0.0; n * n * n] })
    }

    fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        z * self.n * self.n + y * self.n + x
    }

    /// The value at a node, or `None` if the node is outside the box.
    pub fn at(&self, x: usize, y: usize, z: usize) -> Option<f64> {
        if x >= self.n || y >= self.n || z >= self.n {
            return None;
        }
        Some(self.phi[self.idx(x, y, z)])
    }

    /// Set a boundary or interior value.
    pub fn set(&mut self, x: usize, y: usize, z: usize, v: f64) -> bool {
        if x >= self.n || y >= self.n || z >= self.n || !v.is_finite() {
            return false;
        }
        let i = self.idx(x, y, z);
        self.phi[i] = v;
        true
    }

    /// Solve `∇²φ = f` with the boundary held fixed.
    ///
    /// `charge` gives `f` at each interior node, in the same indexing as `phi`. The
    /// boundary values must already be set; they are not touched. Jacobi iteration is
    /// used because it is unconditionally stable and its convergence is easy to reason
    /// about — a wrong answer from an over-clever scheme is worse than a slow one.
    ///
    /// Returns the number of sweeps taken, or `None` if the charge field is the wrong
    /// length or contains a non-finite value.
    pub fn solve(&mut self, charge: &[f64], sweeps: usize, tol: f64) -> Option<usize> {
        let n = self.n;
        if charge.len() != self.phi.len() || charge.iter().any(|c| !c.is_finite()) {
            return None;
        }
        let h2 = self.h_nm * self.h_nm;
        let mut next = self.phi.clone();
        for sweep in 0..sweeps {
            let mut worst = 0.0f64;
            for z in 1..n - 1 {
                for y in 1..n - 1 {
                    for x in 1..n - 1 {
                        let i = self.idx(x, y, z);
                        // the six-point stencil: the average of the neighbours, less the
                        // source term scaled by h^2/6
                        let sum = self.phi[self.idx(x - 1, y, z)]
                            + self.phi[self.idx(x + 1, y, z)]
                            + self.phi[self.idx(x, y - 1, z)]
                            + self.phi[self.idx(x, y + 1, z)]
                            + self.phi[self.idx(x, y, z - 1)]
                            + self.phi[self.idx(x, y, z + 1)];
                        let v = sum / 6.0 - charge[i] * h2 / 6.0;
                        worst = worst.max((v - self.phi[i]).abs());
                        next[i] = v;
                    }
                }
            }
            self.phi.copy_from_slice(&next);
            if worst < tol {
                return Some(sweep + 1);
            }
        }
        Some(sweeps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;

    /// The exact solution used below: sin(pi x) sin(pi y) sin(pi z) on the unit box,
    /// whose Laplacian is -3 pi^2 times itself. Zero on every boundary face.
    fn exact(x: f64, y: f64, z: f64) -> f64 {
        (PI * x).sin() * (PI * y).sin() * (PI * z).sin()
    }

    /// Build the unit box with the exact solution on the boundary and the matching charge.
    fn box_with_exact(n: usize) -> (Volume3D, Vec<f64>) {
        let h = 1.0 / (n - 1) as f64;
        let mut v = Volume3D::new(n, h).expect("a valid box");
        let mut charge = vec![0.0; n * n * n];
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let (fx, fy, fz) = (x as f64 * h, y as f64 * h, z as f64 * h);
                    let e = exact(fx, fy, fz);
                    let boundary = x == 0 || y == 0 || z == 0 || x == n - 1 || y == n - 1 || z == n - 1;
                    if boundary {
                        v.set(x, y, z, e);
                    } else {
                        let i = z * n * n + y * n + x;
                        charge[i] = -3.0 * PI * PI * e;
                    }
                }
            }
        }
        (v, charge)
    }

    #[test]
    fn the_solver_reproduces_a_function_whose_laplacian_is_known() {
        // the test that makes this a solver rather than a smoother: a field with a
        // non-zero right-hand side, checked against the exact answer
        let n = 17;
        let (mut v, charge) = box_with_exact(n);
        let sweeps = v.solve(&charge, 4000, 1e-12).expect("solves");
        assert!(sweeps > 0);
        let h = 1.0 / (n - 1) as f64;
        let mut worst = 0.0f64;
        for z in 1..n - 1 {
            for y in 1..n - 1 {
                for x in 1..n - 1 {
                    let want = exact(x as f64 * h, y as f64 * h, z as f64 * h);
                    worst = worst.max((v.at(x, y, z).unwrap() - want).abs());
                }
            }
        }
        assert!(worst < 0.02, "worst interior error {worst:.4}");
    }

    #[test]
    fn a_flat_boundary_with_no_charge_gives_a_flat_interior() {
        // the maximum principle: a harmonic function with constant boundary is constant
        let n = 9;
        let mut v = Volume3D::new(n, 1.0).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    let boundary = x == 0 || y == 0 || z == 0 || x == n - 1 || y == n - 1 || z == n - 1;
                    if boundary {
                        v.set(x, y, z, 5.0);
                    }
                }
            }
        }
        v.solve(&vec![0.0; n * n * n], 500, 1e-12).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    assert!((v.at(x, y, z).unwrap() - 5.0).abs() < 1e-9);
                }
            }
        }
    }

    #[test]
    fn the_boundary_is_never_disturbed_by_the_solve() {
        let n = 9;
        let mut v = Volume3D::new(n, 0.125).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    if x == 0 || y == 0 || z == 0 || x == n - 1 || y == n - 1 || z == n - 1 {
                        v.set(x, y, z, (x + y + z) as f64);
                    }
                }
            }
        }
        let before = v.phi.clone();
        v.solve(&vec![1.0; n * n * n], 20, 0.0).unwrap();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n {
                    if x == 0 || y == 0 || z == 0 || x == n - 1 || y == n - 1 || z == n - 1 {
                        assert_eq!(v.at(x, y, z).unwrap(), before[z * n * n + y * n + x]);
                    }
                }
            }
        }
    }

    #[test]
    fn a_box_too_small_or_badly_spaced_is_refused() {
        assert!(Volume3D::new(2, 1.0).is_none(), "two nodes have no interior");
        assert!(Volume3D::new(8, 0.0).is_none());
        assert!(Volume3D::new(8, -1.0).is_none());
        assert!(Volume3D::new(8, f64::NAN).is_none());
    }

    #[test]
    fn a_charge_field_of_the_wrong_shape_is_refused() {
        let mut v = Volume3D::new(5, 1.0).unwrap();
        assert!(v.solve(&[0.0; 10], 10, 1e-9).is_none());
        assert!(v.solve(&[f64::NAN; 125], 10, 1e-9).is_none());
    }

    #[test]
    fn nodes_outside_the_box_are_not_addressable() {
        let mut v = Volume3D::new(5, 1.0).unwrap();
        assert!(v.at(5, 0, 0).is_none());
        assert!(!v.set(0, 0, 5, 1.0));
        assert!(!v.set(0, 0, 0, f64::INFINITY));
    }
}
