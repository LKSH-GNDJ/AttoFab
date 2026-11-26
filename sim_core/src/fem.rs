//! 2D plane-strain finite-element stress solver on the SVO cross-section.
//!
//! * Bilinear quadrilateral (Q4) elements, 2×2 Gauss integration, rectangular cells.
//! * Matrix-free Jacobi-preconditioned conjugate gradients, so the whole
//!   high-resolution tile (≈10⁴ elements) solves in well under a second in release.
//! * Intrinsic film stress enters as an **eigenstrain** `e` acting along the wafer
//!   plane (x and the out-of-plane axis), not through the thickness. A film that
//!   is laterally constrained and free on top then carries `σ = −E·e / (1 − ν)`.
//! * The mesh is built from the same column stacks the process operators use, so
//!   the stress field is computed on the geometry that was actually fabricated.
//!
//! Scope: this is a *local* cross-section model (tile bottom clamped, side walls
//! on rollers = laterally periodic). Wafer-scale bow is still the Stoney
//! relation in `physics::deposition`. The solver is verified against Stoney
//! curvature with a mesh-refinement study (see tests); it is not calibrated.

use crate::geom::{material_at_z, ColumnGrid, GeomError, ProcessWindow};
use crate::wafer::{Material, WaferState};
use thiserror::Error;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Elastic { pub e_pa: f64, pub nu: f64 }

impl Elastic {
    pub fn lame(&self) -> (f64, f64) {
        let (e, nu) = (self.e_pa, self.nu);
        (e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu)), e / (2.0 * (1.0 + nu)))
    }
    /// Biaxial modulus `E / (1 − ν)`: the stress of a film that is laterally constrained in both
    /// in-plane directions and free through its thickness, per unit in-plane eigenstrain.
    pub fn eigen_modulus(&self) -> f64 { self.e_pa / (1.0 - self.nu) }
}

/// Isotropic room-temperature textbook values (GPa → Pa). Not calibrated.
pub fn elastic(m: Material) -> Option<Elastic> {
    let (e, nu) = match m {
        Material::Air => return None,
        Material::Silicon | Material::Dopant(_) => (130.0, 0.28),
        Material::SiliconDioxide => (70.0, 0.17),
        Material::SiliconCarbide => (410.0, 0.19),
        Material::GalliumNitride => (210.0, 0.25),
        Material::Photoresist => (3.0, 0.35),
        Material::Metal => (70.0, 0.35),
    };
    Some(Elastic { e_pa: e * 1.0e9, nu })
}

#[derive(Debug, Error, PartialEq)]
pub enum FemError {
    #[error("mesh has no solid elements")]
    Empty,
    #[error("non-finite input")]
    NonFinite,
    #[error("conjugate gradients did not converge: residual {residual:.2e} after {iterations} iterations")]
    NotConverged { iterations: usize, residual: f64 },
    #[error("{0}")]
    Geometry(String),
}

impl From<GeomError> for FemError { fn from(e: GeomError) -> Self { FemError::Geometry(e.to_string()) } }

#[derive(Debug, Clone)]
pub struct Mesh2D {
    pub nx: usize,
    pub ny: usize,
    pub dx: f64,
    pub dy: f64,
    /// Row-major, `j * nx + i`, `j = 0` at the bottom; `None` is void.
    pub cells: Vec<Option<Elastic>>,
    pub material: Vec<Material>,
    /// In-plane eigenstrain per cell (dimensionless; negative = film wants to shrink = tensile).
    pub eigenstrain: Vec<f64>,
}

impl Mesh2D {
    pub fn new(nx: usize, ny: usize, dx: f64, dy: f64) -> Self {
        Self { nx, ny, dx, dy, cells: vec![None; nx * ny], material: vec![Material::Air; nx * ny], eigenstrain: vec![0.0; nx * ny] }
    }
    pub fn set(&mut self, i: usize, j: usize, m: Material, eigen: f64) {
        let k = j * self.nx + i;
        self.cells[k] = elastic(m);
        self.material[k] = m;
        self.eigenstrain[k] = if self.cells[k].is_some() { eigen } else { 0.0 };
    }
    pub fn node(&self, i: usize, j: usize) -> usize { j * (self.nx + 1) + i }
    pub fn dofs(&self) -> usize { 2 * (self.nx + 1) * (self.ny + 1) }
}

#[derive(Debug, Clone)]
pub struct Constraints { pub fixed: Vec<bool> }

impl Constraints {
    pub fn none(mesh: &Mesh2D) -> Self { Self { fixed: vec![false; mesh.dofs()] } }
    /// Bottom edge clamped (embedded in the thick substrate) and both side walls on rollers.
    pub fn tile(mesh: &Mesh2D) -> Self {
        let mut c = Self::none(mesh);
        for i in 0..=mesh.nx { let n = mesh.node(i, 0); c.fixed[2 * n] = true; c.fixed[2 * n + 1] = true; }
        for j in 0..=mesh.ny { c.fixed[2 * mesh.node(0, j)] = true; c.fixed[2 * mesh.node(mesh.nx, j)] = true; }
        c
    }
    /// Symmetric half-strip: roller on the symmetry plane `x = 0`, one vertical pin (no rigid-body motion).
    pub fn symmetric_half(mesh: &Mesh2D) -> Self {
        let mut c = Self::none(mesh);
        for j in 0..=mesh.ny { c.fixed[2 * mesh.node(0, j)] = true; }
        c.fixed[2 * mesh.node(0, 0) + 1] = true;
        c
    }
}

/// 8×8 plane-strain Q4 stiffness for a `dx × dy` rectangle (unit depth).
fn ke(dx: f64, dy: f64, el: Elastic) -> [[f64; 8]; 8] {
    let (l, m) = el.lame();
    let d = [[l + 2.0 * m, l, 0.0], [l, l + 2.0 * m, 0.0], [0.0, 0.0, m]];
    let g = 1.0 / 3f64.sqrt();
    let mut k = [[0.0; 8]; 8];
    for (xi, eta) in [(-g, -g), (g, -g), (g, g), (-g, g)] {
        let b = bmat(dx, dy, xi, eta);
        let w = dx * dy / 4.0; // |J| for a rectangle, Gauss weights = 1
        for a in 0..8 { for c in 0..8 {
            let mut s = 0.0;
            for p in 0..3 { for q in 0..3 { s += b[p][a] * d[p][q] * b[q][c]; } }
            k[a][c] += s * w;
        } }
    }
    k
}

fn bmat(dx: f64, dy: f64, xi: f64, eta: f64) -> [[f64; 8]; 3] {
    let sx = [-1.0, 1.0, 1.0, -1.0];
    let sy = [-1.0, -1.0, 1.0, 1.0];
    let mut b = [[0.0; 8]; 3];
    for n in 0..4 {
        let dnx = sx[n] * (1.0 + sy[n] * eta) / (4.0) * (2.0 / dx);
        let dny = sy[n] * (1.0 + sx[n] * xi) / (4.0) * (2.0 / dy);
        b[0][2 * n] = dnx;
        b[1][2 * n + 1] = dny;
        b[2][2 * n] = dny;
        b[2][2 * n + 1] = dnx;
    }
    b
}

#[derive(Debug, Clone)]
pub struct FemResult {
    pub ux: Vec<f64>,
    pub uy: Vec<f64>,
    /// Per-cell in-plane stresses and von Mises (3D, with the out-of-plane σyy), MPa.
    pub sxx_mpa: Vec<f64>,
    pub szz_mpa: Vec<f64>,
    pub von_mises_mpa: Vec<f64>,
    pub iterations: usize,
    pub residual: f64,
}

impl FemResult {
    pub fn max_von_mises_mpa(&self) -> f64 { self.von_mises_mpa.iter().cloned().fold(0.0, f64::max) }
    pub fn uy_at(&self, mesh: &Mesh2D, i: usize, j: usize) -> f64 { self.uy[mesh.node(i, j)] }
}

pub fn solve(mesh: &Mesh2D, bc: &Constraints, rel_tol: f64, max_iter: usize) -> Result<FemResult, FemError> {
    if !mesh.eigenstrain.iter().all(|v| v.is_finite()) || !(mesh.dx > 0.0 && mesh.dy > 0.0) { return Err(FemError::NonFinite); }
    if mesh.cells.iter().all(|c| c.is_none()) { return Err(FemError::Empty); }
    let n = mesh.dofs();
    // Per-material stiffness cache (few distinct materials).
    let mut kes: Vec<(Elastic, [[f64; 8]; 8])> = Vec::new();
    let mut kidx = vec![usize::MAX; mesh.cells.len()];
    for (c, cell) in mesh.cells.iter().enumerate() {
        if let Some(el) = cell {
            let id = match kes.iter().position(|(e, _)| e == el) { Some(i) => i, None => { kes.push((*el, ke(mesh.dx, mesh.dy, *el))); kes.len() - 1 } };
            kidx[c] = id;
        }
    }
    let enodes = |i: usize, j: usize| [mesh.node(i, j), mesh.node(i + 1, j), mesh.node(i + 1, j + 1), mesh.node(i, j + 1)];
    // Nodes touched by no solid element have no stiffness: pin them.
    let mut fixed = bc.fixed.clone();
    let mut touched = vec![false; n / 2];
    for j in 0..mesh.ny { for i in 0..mesh.nx { if kidx[j * mesh.nx + i] != usize::MAX { for nd in enodes(i, j) { touched[nd] = true; } } } }
    for (nd, t) in touched.iter().enumerate() { if !t { fixed[2 * nd] = true; fixed[2 * nd + 1] = true; } }

    // Load vector from eigenstress s* = [(2λ+2μ)e, 2λe, 0].
    let mut f = vec![0.0; n];
    let g = 1.0 / 3f64.sqrt();
    for j in 0..mesh.ny { for i in 0..mesh.nx {
        let c = j * mesh.nx + i;
        let Some(el) = mesh.cells[c] else { continue };
        let e = mesh.eigenstrain[c];
        if e == 0.0 { continue; }
        let (l, m) = el.lame();
        let s = [(2.0 * l + 2.0 * m) * e, 2.0 * l * e, 0.0];
        let nodes = enodes(i, j);
        for (xi, eta) in [(-g, -g), (g, -g), (g, g), (-g, g)] {
            let b = bmat(mesh.dx, mesh.dy, xi, eta);
            let w = mesh.dx * mesh.dy / 4.0;
            for a in 0..8 { let mut v = 0.0; for p in 0..3 { v += b[p][a] * s[p]; } f[2 * nodes[a / 2] + a % 2] += v * w; }
        }
    } }
    for d in 0..n { if fixed[d] { f[d] = 0.0; } }

    let matvec = |x: &[f64], y: &mut [f64]| {
        for v in y.iter_mut() { *v = 0.0; }
        for j in 0..mesh.ny { for i in 0..mesh.nx {
            let id = kidx[j * mesh.nx + i];
            if id == usize::MAX { continue; }
            let k = &kes[id].1;
            let nodes = enodes(i, j);
            let mut xe = [0.0; 8];
            for a in 0..8 { let d = 2 * nodes[a / 2] + a % 2; xe[a] = if fixed[d] { 0.0 } else { x[d] }; }
            for a in 0..8 { let mut s = 0.0; for c in 0..8 { s += k[a][c] * xe[c]; } y[2 * nodes[a / 2] + a % 2] += s; }
        } }
        for d in 0..n { if fixed[d] { y[d] = x[d]; } }
    };
    // Jacobi preconditioner.
    let mut diag = vec![0.0; n];
    for j in 0..mesh.ny { for i in 0..mesh.nx {
        let id = kidx[j * mesh.nx + i];
        if id == usize::MAX { continue; }
        let nodes = enodes(i, j);
        for a in 0..8 { diag[2 * nodes[a / 2] + a % 2] += kes[id].1[a][a]; }
    } }
    for d in 0..n { if fixed[d] || diag[d] <= 0.0 { diag[d] = 1.0; } }

    let mut x = vec![0.0; n];
    let mut r = f.clone();
    let mut z: Vec<f64> = r.iter().zip(&diag).map(|(a, d)| a / d).collect();
    let mut p = z.clone();
    let mut rz: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
    let bnorm = f.iter().map(|v| v * v).sum::<f64>().sqrt().max(1e-300);
    let mut ap = vec![0.0; n];
    let (mut it, mut res) = (0usize, 1.0);
    while it < max_iter {
        matvec(&p, &mut ap);
        let pap: f64 = p.iter().zip(&ap).map(|(a, b)| a * b).sum();
        if pap.abs() < 1e-300 { break; }
        let alpha = rz / pap;
        for d in 0..n { x[d] += alpha * p[d]; r[d] -= alpha * ap[d]; }
        res = r.iter().map(|v| v * v).sum::<f64>().sqrt() / bnorm;
        it += 1;
        if res < rel_tol { break; }
        for d in 0..n { z[d] = r[d] / diag[d]; }
        let rz_new: f64 = r.iter().zip(&z).map(|(a, b)| a * b).sum();
        let beta = rz_new / rz; rz = rz_new;
        for d in 0..n { p[d] = z[d] + beta * p[d]; }
    }
    if f.iter().all(|v| *v == 0.0) { res = 0.0; }
    if res >= rel_tol { return Err(FemError::NotConverged { iterations: it, residual: res }); }

    // Stress recovery at the cell centre.
    let mut sxx = vec![0.0; mesh.cells.len()];
    let mut szz = vec![0.0; mesh.cells.len()];
    let mut vm = vec![0.0; mesh.cells.len()];
    for j in 0..mesh.ny { for i in 0..mesh.nx {
        let c = j * mesh.nx + i;
        let Some(el) = mesh.cells[c] else { continue };
        let (l, m) = el.lame();
        let nodes = enodes(i, j);
        let b = bmat(mesh.dx, mesh.dy, 0.0, 0.0);
        let mut eps = [0.0; 3];
        for a in 0..8 { let u = x[2 * nodes[a / 2] + a % 2]; for p in 0..3 { eps[p] += b[p][a] * u; } }
        let e = mesh.eigenstrain[c];
        let sx = (l + 2.0 * m) * eps[0] + l * eps[1] - (2.0 * l + 2.0 * m) * e;
        let sz = l * eps[0] + (l + 2.0 * m) * eps[1] - 2.0 * l * e;
        let tau = m * eps[2];
        let sy = l * (eps[0] + eps[1]) - 2.0 * l * e - 2.0 * m * e; // out-of-plane (plane strain, eigenstrain e along y)
        let v = (0.5 * ((sx - sy).powi(2) + (sy - sz).powi(2) + (sz - sx).powi(2)) + 3.0 * tau * tau).sqrt();
        sxx[c] = sx / 1.0e6; szz[c] = sz / 1.0e6; vm[c] = v / 1.0e6;
    } }
    let (ux, uy) = (0..n / 2).map(|nd| (x[2 * nd], x[2 * nd + 1])).unzip();
    Ok(FemResult { ux, uy, sxx_mpa: sxx, szz_mpa: szz, von_mises_mpa: vm, iterations: it, residual: res })
}

/// Mesh of the x–z cross-section through the window's centre row, straight from the SVO.
/// `film_stress_mpa` lists tensile (+) / compressive (−) intrinsic stress per film material.
pub fn mesh_from_wafer(wafer: &WaferState, window: &ProcessWindow, film_stress_mpa: &[(Material, f64)]) -> Result<Mesh2D, FemError> {
    let grid = ColumnGrid::sample(wafer, window)?;
    let row = grid.ny / 2;
    let nz = ((grid.z_top - grid.z_bottom) / grid.pitch).round() as usize;
    let mut mesh = Mesh2D::new(grid.nx, nz, grid.pitch, grid.pitch);
    for i in 0..grid.nx {
        let k = grid.idx(i, row);
        if !grid.active[k] { continue; }
        for j in 0..nz {
            let z = grid.z_bottom + (j as f64 + 0.5) * grid.pitch;
            let m = material_at_z(&grid.cols[k], z);
            // lint:allow no matching film in the stack, so no stress contribution
            let eigen = film_stress_mpa.iter().rev().find(|(fm, _)| *fm == m).and_then(|(fm, s)| elastic(*fm).map(|el| -s * 1.0e6 / el.eigen_modulus())).unwrap_or(0.0);
            mesh.set(i, j, m, eigen);
        }
    }
    Ok(mesh)
}

#[derive(Debug, Clone, Copy)]
pub struct StressSummary { pub max_von_mises_mpa: f64, pub max_film_von_mises_mpa: f64, pub uniform_film_mpa: f64, pub iterations: usize }

/// Solve the tile cross-section for the stress state left by the recorded film stresses.
pub fn tile_stress(wafer: &WaferState, window: &ProcessWindow) -> Result<Option<StressSummary>, FemError> {
    let films: Vec<(Material, f64)> = wafer.meta.film_stress.iter().map(|f| (f.material, f.stress_mpa)).collect();
    if films.is_empty() { return Ok(None); }
    let mesh = mesh_from_wafer(wafer, window, &films)?;
    let res = solve(&mesh, &Constraints::tile(&mesh), 1e-9, 20_000)?;
    let film_max = res.von_mises_mpa.iter().zip(&mesh.material).filter(|(_, m)| films.iter().any(|(f, _)| f == *m)).map(|(v, _)| *v).fold(0.0, f64::max);
    // lint:allow no matching film in the stack, so no stress contribution
    let uniform = films.last().map(|f| f.1.abs()).unwrap_or(0.0);
    Ok(Some(StressSummary { max_von_mises_mpa: res.max_von_mises_mpa(), max_film_von_mises_mpa: film_max, uniform_film_mpa: uniform, iterations: res.iterations }))
}

/// Closed-form and strip benchmarks used by the tests and the V&V suite.
pub mod benchmark {
    use super::*;

    /// Film on a substrate strip (half model). Returns the FEM curvature from the bottom-edge deflection.
    pub fn strip_curvature(hs: f64, hf: f64, length: f64, film_cells: usize, sigma_f_pa: f64) -> f64 {
        let film = elastic(Material::SiliconDioxide).unwrap();
        let dy = hf / film_cells as f64; // the film thickness is represented exactly at every refinement level
        let cells_across_sub = (hs / dy).round() as usize;
        let nf = film_cells;
        let nx = ((length / 2.0) / dy).round() as usize;
        let mut mesh = Mesh2D::new(nx, cells_across_sub + nf, dy, dy);
        let e = -sigma_f_pa / film.eigen_modulus();
        for j in 0..mesh.ny { for i in 0..nx { if j < cells_across_sub { mesh.set(i, j, Material::Silicon, 0.0) } else { mesh.set(i, j, Material::SiliconDioxide, e) } } }
        let res = solve(&mesh, &Constraints::symmetric_half(&mesh), 1e-11, 200_000).expect("benchmark strip solves");
        // Fit uy(x) = c0 + c2 x² over the central 60% of the half-length (away from the free end).
        let (mut s4, mut s2, mut s0, mut t2, mut t0) = (0.0, 0.0, 0.0, 0.0, 0.0);
        let last = (nx as f64 * 0.6) as usize;
        for i in 0..=last {
            let x = i as f64 * dy; let u = res.uy_at(&mesh, i, 0);
            s4 += x.powi(4); s2 += x * x; s0 += 1.0; t2 += x * x * u; t0 += u;
        }
        2.0 * (t2 * s0 - t0 * s2) / (s4 * s0 - s2 * s2)
    }

    /// Thin-film (Stoney) limit, plane-strain substrate modulus.
    pub fn stoney(hs: f64, hf: f64, sigma_f: f64) -> f64 {
        let s = elastic(Material::Silicon).unwrap();
        6.0 * sigma_f * hf / (s.e_pa / (1.0 - s.nu * s.nu) * hs * hs)
    }

    /// Exact Timoshenko two-layer beam (plane strain, in-plane eigenstrain in the film only).
    pub fn timoshenko(hs: f64, hf: f64, sigma_f: f64) -> f64 {
        let (s, f) = (elastic(Material::Silicon).unwrap(), elastic(Material::SiliconDioxide).unwrap());
        let (ms, mf) = (s.e_pa / (1.0 - s.nu * s.nu), f.e_pa / (1.0 - f.nu * f.nu));
        let eps = sigma_f / (f.e_pa / (1.0 - f.nu)) * (1.0 + f.nu); // effective in-plane free strain magnitude
        let d = ms * ms * hs.powi(4) + 4.0 * mf * ms * hf * hs.powi(3) + 6.0 * mf * ms * hf * hf * hs * hs + 4.0 * mf * ms * hf.powi(3) * hs + mf * mf * hf.powi(4);
        6.0 * eps * mf * ms * hf * hs * (hf + hs) / d
    }


}

#[cfg(test)]
mod tests {
    use super::*;
    use super::benchmark::{stoney, strip_curvature, timoshenko};
    use crate::validation::gci::{calculate_gci, GridResult};

    #[test]
    fn free_expansion_of_a_uniform_block_is_stress_free_in_plane() {
        let mut mesh = Mesh2D::new(12, 8, 10.0, 10.0);
        for j in 0..8 { for i in 0..12 { mesh.set(i, j, Material::Silicon, -1.0e-3); } }
        let r = solve(&mesh, &Constraints::symmetric_half(&mesh), 1e-12, 10_000).unwrap();
        let scale = elastic(Material::Silicon).unwrap().eigen_modulus() * 1.0e-3 / 1.0e6;
        for k in 0..96 { assert!(r.sxx_mpa[k].abs() < 1e-6 * scale && r.szz_mpa[k].abs() < 1e-6 * scale, "cell {k}: {} {}", r.sxx_mpa[k], r.szz_mpa[k]); }
    }

    #[test]
    fn constrained_film_carries_exactly_the_prescribed_stress() {
        // Rollers on both sides, clamped bottom: εxx = 0 everywhere, so σxx = −(2λ+2μ)·e = σ_f in the film.
        let (nx, ny) = (10, 6);
        let mut mesh = Mesh2D::new(nx, ny, 5.0, 5.0);
        let film = elastic(Material::SiliconDioxide).unwrap();
        let sigma = 300.0e6;
        for j in 0..ny { for i in 0..nx { mesh.set(i, j, Material::SiliconDioxide, -sigma / film.eigen_modulus()); } }
        let mut bc = Constraints::tile(&mesh);
        for j in 0..=ny { for i in 0..=nx { let n = mesh.node(i, j); bc.fixed[2 * n] = true; } }
        let r = solve(&mesh, &bc, 1e-12, 10_000).unwrap();
        for k in 0..nx * ny { assert!((r.sxx_mpa[k] - 300.0).abs() < 1e-6, "{}", r.sxx_mpa[k]); }
    }

    #[test]
    fn matches_the_exact_two_layer_beam_solution_and_converges_under_mesh_refinement() {
        let (hs, hf, l, sigma) = (100.0, 10.0, 1200.0, 200.0e6);
        let exact = timoshenko(hs, hf, sigma);
        let k = [1usize, 2, 3].map(|n| strip_curvature(hs, hf, l, n, sigma));
        assert!(k[2] > 0.0, "a tensile film on top curves the strip concave toward the film");
        let errs: Vec<f64> = k.iter().map(|v| (v - exact).abs() / exact).collect();
        assert!(errs[2] < 0.05, "FEM {:.4e} vs exact beam {:.4e}: errors by mesh {errs:?}", k[2], exact);
        assert!(errs[2] <= errs[0] + 1e-9, "refinement must not worsen the error: {errs:?}");
        // Thin-film Stoney is the limit of the exact solution (hf/hs = 5%): within the expected O(hf/hs) band.
        assert!((stoney(hs, hf, sigma) - exact).abs() / exact < 0.25, "Stoney is only the thin-film limit (hf/hs = 10%)");
        let gci = calculate_gci(GridResult { spacing: hf, value: k[0] }, GridResult { spacing: hf / 2.0, value: k[1] }, GridResult { spacing: hf / 3.0, value: k[2] }, 1.25);
        if let Ok(g) = gci { assert!(g.fine_gci_percent < 6.0, "fine-grid GCI {:.2}%", g.fine_gci_percent); }
    }

    #[test]
    fn film_edges_concentrate_stress_in_a_patterned_tile() {
        use crate::pipeline::{run_chain};
        use crate::params::StageParams;
        use crate::physics::Stage;
        let w = WaferState::new_wafer_tile(300.0, 500.0, 16.0, 7, (0.0, 0.0)).unwrap();
        let win = ProcessWindow::default();
        let chain = run_chain(&w, &win, &[StageParams::default_for(Stage::Deposition), StageParams::default_for(Stage::Lithography), StageParams::default_for(Stage::Etching)], 1).unwrap();
        let s = tile_stress(&chain[2].wafer, &win).unwrap().expect("film stress recorded");
        assert!(s.max_film_von_mises_mpa > 0.0 && s.max_von_mises_mpa.is_finite());
        assert!(s.max_film_von_mises_mpa > 0.5 * s.uniform_film_mpa, "edges and lines carry stress: {s:?}");
    }

    #[test]
    fn empty_mesh_and_non_finite_loads_are_rejected() {
        let m = Mesh2D::new(4, 4, 1.0, 1.0);
        assert_eq!(solve(&m, &Constraints::none(&m), 1e-9, 10).unwrap_err(), FemError::Empty);
        let mut m2 = Mesh2D::new(2, 2, 1.0, 1.0);
        m2.set(0, 0, Material::Silicon, f64::NAN);
        assert_eq!(solve(&m2, &Constraints::symmetric_half(&m2), 1e-9, 10).unwrap_err(), FemError::NonFinite);
    }
}
