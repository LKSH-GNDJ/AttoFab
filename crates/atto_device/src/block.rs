// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! A block-tridiagonal solver for coupled systems (build PRD V-01).
//!
//! The device equations couple three unknowns per node — potential and both
//! carrier densities — so the Jacobian is tridiagonal in 3x3 blocks rather than
//! in scalars. This solves that system by block Thomas elimination.

/// One 3x3 block.
pub type Block = [[f64; 3]; 3];

/// Solve a block-tridiagonal system `A_i x_{i-1} + B_i x_i + C_i x_{i+1} = r_i`
/// for `n` blocks of 3 unknowns. Returns `None` if a pivot block is singular.
pub fn solve(lower: &[Block], diag: &[Block], upper: &[Block], rhs: &[[f64; 3]]) -> Option<Vec<[f64; 3]>> {
    let n = diag.len();
    if n == 0 || lower.len() != n || upper.len() != n || rhs.len() != n {
        return None;
    }
    let mut cp: Vec<Block> = vec![[[0.0; 3]; 3]; n];
    let mut dp: Vec<[f64; 3]> = vec![[0.0; 3]; n];

    // first block
    let mut b = diag[0];
    let mut r = rhs[0];
    cp[0] = upper[0];
    // dp[0] = B^-1 r, and cp[0] = B^-1 C
    let (binv, det) = inv3(&b)?;
    if det.abs() < 1e-300 {
        return None;
    }
    dp[0] = mul3(binv, r);
    cp[0] = mul_block(binv, upper[0]);
    let _ = &mut b;
    let _ = &mut r;

    for i in 1..n {
        // M = A_i * cp[i-1];  B' = B_i - M;  r' = r_i - A_i * dp[i-1]
        let m = mul_block(lower[i], cp[i - 1]);
        let mut bi = diag[i];
        for a in 0..3 {
            for c in 0..3 {
                bi[a][c] -= m[a][c];
            }
        }
        let mut ri = rhs[i];
        let ad = mul3(lower[i], dp[i - 1]);
        for a in 0..3 {
            ri[a] -= ad[a];
        }
        let (biinv, det) = inv3(&bi)?;
        if det.abs() < 1e-300 {
            return None;
        }
        dp[i] = mul3(biinv, ri);
        cp[i] = if i + 1 < n { mul_block(biinv, upper[i]) } else { [[0.0; 3]; 3] };
    }

    let mut x: Vec<[f64; 3]> = vec![[0.0; 3]; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        let cx = mul3(cp[i], x[i + 1]);
        for a in 0..3 {
            x[i][a] = dp[i][a] - cx[a];
        }
    }
    Some(x)
}

fn mul3(m: Block, v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
    ]
}

fn mul_block(a: Block, b: Block) -> Block {
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    out
}

/// Inverse of a 3x3 block by cofactors. Returns the inverse and the determinant.
fn inv3(m: &Block) -> Option<(Block, f64)> {
    let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    if !det.is_finite() || det.abs() < 1e-300 {
        return None;
    }
    let mut inv = [[0.0; 3]; 3];
    inv[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
    inv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
    inv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
    inv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
    inv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
    inv[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
    inv[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
    inv[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
    inv[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
    Some((inv, det))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eye(k: f64) -> Block {
        [[k, 0.0, 0.0], [0.0, k, 0.0], [0.0, 0.0, k]]
    }

    #[test]
    fn a_single_block_is_just_a_3x3_solve() {
        let b = [[2.0, 1.0, 0.0], [1.0, 3.0, 1.0], [0.0, 1.0, 2.0]];
        let r = [[1.0, 2.0, 3.0]];
        let x = solve(&[eye(0.0)], &[b], &[eye(0.0)], &r).expect("solvable");
        // verify by substitution
        let back = mul3(b, x[0]);
        for i in 0..3 {
            assert!((back[i] - r[0][i]).abs() < 1e-12, "residual {i}");
        }
    }

    #[test]
    fn two_coupled_blocks_solve_and_substitute() {
        let a = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let lower = [a, [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]];
        let diag = [
            [[4.0, 1.0, 0.0], [1.0, 4.0, 1.0], [0.0, 1.0, 4.0]],
            [[3.0, 0.5, 0.0], [0.5, 3.0, 0.5], [0.0, 0.5, 3.0]],
        ];
        let upper = [eye(1.0), eye(0.0)];
        let rhs = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]];
        let x = solve(&lower, &diag, &upper, &rhs).expect("solvable");
        // check both rows
        let r0 = mul3(diag[0], x[0]);
        let u0 = mul3(upper[0], x[1]);
        for i in 0..3 {
            assert!((r0[i] + u0[i] - rhs[0][i]).abs() < 1e-10, "row 0, {i}");
        }
        let l1 = mul3(lower[1], x[0]);
        let r1 = mul3(diag[1], x[1]);
        for i in 0..3 {
            assert!((l1[i] + r1[i] - rhs[1][i]).abs() < 1e-10, "row 1, {i}");
        }
    }

    #[test]
    fn a_singular_block_is_reported_not_papered_over() {
        let zero = [[0.0; 3]; 3];
        assert!(solve(&[zero], &[zero], &[zero], &[[1.0, 1.0, 1.0]]).is_none());
    }

    #[test]
    fn mismatched_lengths_are_rejected() {
        let z = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(solve(&[z], &[z, z], &[z, z], &[[0.0; 3], [0.0; 3]]).is_none());
    }
}
