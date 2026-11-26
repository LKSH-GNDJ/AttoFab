//! Metrology on the committed wafer state (PRD §2 stage 7).
//!
//! Every value is *measured from geometry* (column stacks of the SVO) or from
//! recorded process metadata; nothing is echoed back from stage inputs. The
//! report carries the fingerprint of the exact wafer it measured.

use crate::geom::{solid_top, ColumnGrid, ProcessWindow};
use crate::params::{MetrologyParams, Quantity, Target};
use crate::physics::metrology as ellipso;
use crate::wafer::{Material, WaferState};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Status { Pass, Fail, NotMeasured }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Measurement {
    pub name: String,
    pub unit: String,
    pub value: Option<f64>,
    pub nominal: f64,
    pub tolerance: f64,
    pub deviation: Option<f64>,
    pub status: Status,
    /// Per-site values (row-major `sites × sites`), where the quantity is mapped.
    pub site_values: Vec<Option<f64>>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct WaferMap {
    pub label: String,
    pub n: usize,
    pub values: Vec<Option<f64>>,
    pub nominal: Option<f64>,
    pub tolerance: Option<f64>,
}

impl WaferMap {
    pub fn out_of_spec(&self, i: usize) -> bool {
        match (self.values.get(i).copied().flatten(), self.nominal, self.tolerance) {
            (Some(v), Some(n), Some(t)) => (v - n).abs() > t,
            _ => false,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Ellipsometry { pub psi_deg: f64, pub delta_deg: f64, pub film_nm: f64, pub film: Material }

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct MetrologyReport {
    /// Fingerprint of the wafer that was measured.
    pub wafer_fingerprint: u64,
    pub window: ProcessWindow,
    pub measurements: Vec<Measurement>,
    pub wafer_map: WaferMap,
    pub optical: Option<Ellipsometry>,
    /// `None` when no target could be measured.
    pub pass: Option<bool>,
}

// --- small deterministic RNG for seeded measurement noise ------------------
fn splitmix(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}
pub fn gaussian(seed: u64, idx: u64) -> f64 {
    let a = (splitmix(seed ^ idx.wrapping_mul(2)) >> 11) as f64 / (1u64 << 53) as f64;
    let b = (splitmix(seed ^ idx.wrapping_mul(2).wrapping_add(1)) >> 11) as f64 / (1u64 << 53) as f64;
    (-2.0 * a.max(1e-300).ln()).sqrt() * (2.0 * std::f64::consts::PI * b).cos()
}

fn thickness_of(col: &[crate::wafer::Run], m: Material) -> f64 { col.iter().filter(|r| r.material == m).map(|r| r.z1 - r.z0).sum() }

/// Blocks of the grid for an `n × n` site map.
fn site_ranges(g: &ColumnGrid, n: usize) -> Vec<(usize, usize, usize, usize)> {
    let mut out = Vec::new();
    for sj in 0..n { for si in 0..n {
        out.push((si * g.nx / n, (si + 1) * g.nx / n, sj * g.ny / n, (sj + 1) * g.ny / n));
    } }
    out
}

fn site_mean(g: &ColumnGrid, r: (usize, usize, usize, usize), f: impl Fn(usize) -> Option<f64>) -> Option<f64> {
    let (mut s, mut c) = (0.0, 0);
    for iy in r.2..r.3 { for ix in r.0..r.1 {
        let k = g.idx(ix, iy);
        if g.active[k] { if let Some(v) = f(k) { s += v; c += 1; } }
    } }
    if c > 0 { Some(s / c as f64) } else { None }
}

/// Run lengths (nm) of columns containing `m` along the window's centre row,
/// falling back to the centre column. Runs touching the window edge are ignored.
pub fn line_widths(g: &ColumnGrid, m: Material) -> Vec<f64> {
    let scan = |along_x: bool| -> Vec<f64> {
        let (len, other) = if along_x { (g.nx, g.ny / 2) } else { (g.ny, g.nx / 2) };
        let has = |i: usize| { let k = if along_x { g.idx(i, other) } else { g.idx(other, i) }; g.active[k] && g.cols[k].iter().any(|r| r.material == m) };
        let mut out = Vec::new();
        let mut i = 0;
        while i < len {
            if has(i) {
                let s = i;
                while i < len && has(i) { i += 1; }
                if s > 0 && i < len { out.push((i - s) as f64 * g.pitch); }
            } else { i += 1; }
        }
        out
    };
    let x = scan(true);
    if !x.is_empty() { x } else { scan(false) }
}

// --- thin-film ellipsometry (single film on silicon) -----------------------
#[derive(Clone, Copy)]
struct C(f64, f64);
impl C {
    fn add(self, o: C) -> C { C(self.0 + o.0, self.1 + o.1) }
    fn sub(self, o: C) -> C { C(self.0 - o.0, self.1 - o.1) }
    fn mul(self, o: C) -> C { C(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0) }
    fn div(self, o: C) -> C { let d = o.0 * o.0 + o.1 * o.1; C((self.0 * o.0 + self.1 * o.1) / d, (self.1 * o.0 - self.0 * o.1) / d) }
    fn exp(self) -> C { let e = self.0.exp(); C(e * self.1.cos(), e * self.1.sin()) }
    fn sqrt(self) -> C {
        let r = (self.0 * self.0 + self.1 * self.1).sqrt();
        let re = ((r + self.0) * 0.5).max(0.0).sqrt();
        let im = ((r - self.0) * 0.5).max(0.0).sqrt() * if self.1 < 0.0 { -1.0 } else { 1.0 };
        C(re, im)
    }
}

/// Complex refractive index N = n + i·(−k) convention (absorbing k > 0 → Im < 0).
pub fn refractive_index(m: Material) -> (f64, f64) {
    match m {
        Material::Silicon | Material::Dopant(_) => (3.88, 0.02), Material::SiliconDioxide => (1.457, 0.0), Material::SiliconCarbide => (2.65, 0.0),
        Material::GalliumNitride => (2.38, 0.0), Material::Photoresist => (1.64, 0.0), Material::Metal => (1.4, 7.6), Material::Air => (1.0, 0.0),
    }
}

fn nc(m: Material) -> C { let (n, k) = refractive_index(m); C(n, -k) }

/// Reflection coefficients (r_p, r_s) of ambient/film/silicon at wavelength λ.
pub fn film_on_silicon_reflection(film: Material, d_nm: f64, wavelength_nm: f64, incidence_deg: f64) -> ((f64, f64), (f64, f64)) {
    let th0 = incidence_deg.to_radians();
    let (n0, n1, n2) = (C(1.0, 0.0), nc(film), nc(Material::Silicon));
    let sin0 = C(th0.sin(), 0.0);
    let cosf = |n: C| -> C { let s = n0.mul(sin0).div(n); C(1.0, 0.0).sub(s.mul(s)).sqrt() };
    let (c0, c1, c2) = (C(th0.cos(), 0.0), cosf(n1), cosf(n2));
    let rs = |na: C, ca: C, nb: C, cb: C| na.mul(ca).sub(nb.mul(cb)).div(na.mul(ca).add(nb.mul(cb)));
    let rp = |na: C, ca: C, nb: C, cb: C| nb.mul(ca).sub(na.mul(cb)).div(nb.mul(ca).add(na.mul(cb)));
    if d_nm <= 0.0 || film == Material::Silicon || matches!(film, Material::Dopant(_)) {
        let (p, s) = (rp(n0, c0, n2, c2), rs(n0, c0, n2, c2));
        return ((p.0, p.1), (s.0, s.1));
    }
    let beta = C(2.0 * std::f64::consts::PI * d_nm / wavelength_nm, 0.0).mul(n1).mul(c1);
    let phase = C(0.0, -2.0).mul(beta).exp();
    let comb = |r01: C, r12: C| r01.add(r12.mul(phase)).div(C(1.0, 0.0).add(r01.mul(r12).mul(phase)));
    let p = comb(rp(n0, c0, n1, c1), rp(n1, c1, n2, c2));
    let s = comb(rs(n0, c0, n1, c1), rs(n1, c1, n2, c2));
    ((p.0, p.1), (s.0, s.1))
}

fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() { return 0.0; }
    sorted[((sorted.len() - 1) as f64 * q).round() as usize]
}

/// Measure the wafer inside `window`.
pub fn measure(wafer: &WaferState, window: &ProcessWindow, p: &MetrologyParams, seed: u64) -> Result<MetrologyReport, crate::geom::GeomError> {
    let g = ColumnGrid::sample(wafer, window)?;
    let n = p.sites.clamp(1, 15);
    let ranges = site_ranges(&g, n);
    let surface = wafer.surface_z_nm;
    let heights = g.surface_heights();
    let active_h: Vec<f64> = (0..heights.len()).filter(|k| g.active[*k]).map(|k| heights[k]).collect();
    let mut sorted = active_h.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mean_h = active_h.iter().sum::<f64>() / active_h.len().max(1) as f64;
    let mut noise_idx = 0u64;
    let mut noisy = |v: f64, unit_is_nm: bool| -> f64 { noise_idx += 1; if unit_is_nm { v + p.noise_nm * gaussian(seed, noise_idx) } else { v } };

    let mut measurements = Vec::new();
    for t in &p.targets {
        let Target { quantity, nominal, tolerance } = *t;
        let (value, sites): (Option<f64>, Vec<Option<f64>>) = match quantity {
            Quantity::FilmThickness(m) => {
                let th: Vec<f64> = g.cols.iter().enumerate().filter(|(k, _)| g.active[*k]).map(|(_, c)| thickness_of(c, m)).filter(|t| *t > 0.0).collect();
                let v = if th.is_empty() { None } else { Some(noisy(th.iter().sum::<f64>() / th.len() as f64, true)) };
                // Mean thickness over the columns that actually carry the film: a block that straddles a line edge
                // or trench must not read as a thin film.
                let s = ranges.iter().map(|r| site_mean(&g, *r, |k| { let t = thickness_of(&g.cols[k], m); if t > 0.0 { Some(t) } else { None } })).collect();
                (v, s)
            }
            Quantity::LineWidth(m) => {
                let w = line_widths(&g, m);
                (if w.is_empty() { None } else { Some(noisy(w.iter().sum::<f64>() / w.len() as f64, true)) }, vec![])
            }
            Quantity::StepHeight => {
                // audit F-10: adjacent columns, not a 10-90 percentile spread
                // `heights` is indexed by column; `active_h` is compacted to the
                // active ones and would not line up with g.active
                let v = atto_bridge::step_height_nm(&heights, &g.active, g.nx, g.ny)
                    .map(|d| noisy(d, true));
                (v, vec![])
            }
            Quantity::JunctionDepth => {
                let d = g.cols.iter().enumerate().filter(|(k, _)| g.active[*k])
                    .map(|(_, c)| c.iter().filter(|r| matches!(r.material, Material::Dopant(_))).map(|r| r.z1 - r.z0).sum::<f64>()).fold(0.0, f64::max);
                (if d > 0.0 { Some(noisy(d, true)) } else { None }, vec![])
            }
            Quantity::TopographyRms => {
                let var = active_h.iter().map(|h| (h - mean_h) * (h - mean_h)).sum::<f64>() / active_h.len().max(1) as f64;
                (Some(noisy(var.sqrt(), true)), vec![])
            }
            Quantity::FilmStress => (wafer.meta.film_stress.last().map(|f| f.stress_mpa), vec![]),
            // lint:allow no FEM tile covers this window; zero peak stress is the reported state
            Quantity::PeakFilmStress => (crate::fem::tile_stress(wafer, window).ok().flatten().map(|s| s.max_film_von_mises_mpa), vec![]),
            Quantity::ResistResidue => {
                let v: f64 = g.cols.iter().enumerate().filter(|(k, _)| g.active[*k]).map(|(_, c)| thickness_of(c, Material::Photoresist)).sum::<f64>() / g.active_count().max(1) as f64;
                (Some(v), vec![])
            }
        };
        let (deviation, status) = match value {
            Some(v) => (Some(v - nominal), if (v - nominal).abs() <= tolerance { Status::Pass } else { Status::Fail }),
            None => (None, Status::NotMeasured),
        };
        measurements.push(Measurement { name: quantity.label(), unit: quantity.unit().into(), value, nominal, tolerance, deviation, status, site_values: sites });
    }

    // Wafer map: height above the original substrate surface per site.
    let film = p.targets.iter().find_map(|t| if let Quantity::FilmThickness(_) = t.quantity { Some(*t) } else { None });
    let film_sites = measurements.iter().zip(&p.targets).find_map(|(m, t)| if matches!(t.quantity, Quantity::FilmThickness(_)) && !m.site_values.is_empty() { Some(m.site_values.clone()) } else { None });
    // With a film-thickness target the map shows that film's thickness where it exists (so it can be judged
    // against the target); otherwise it shows the mean height above the original surface.
    let (label, map_vals): (String, Vec<Option<f64>>) = match (&film, film_sites) {
        (Some(t), Some(v)) => (format!("{} (nm), mean where present", t.quantity.label()), v),
        _ => ("Height above original surface (nm)".into(), ranges.iter().map(|r| site_mean(&g, *r, |k| Some(heights[k] - surface))).collect()),
    };
    let have_film_map = film.is_some() && label.contains("mean where present");
    let wafer_map = WaferMap { label, n, values: map_vals, nominal: if have_film_map { film.map(|t| t.nominal) } else { None }, tolerance: if have_film_map { film.map(|t| t.tolerance) } else { None } };

    // Optical: film stack at the window centre.
    let optical = {
        let (cx, cy) = (window.center[0], window.center[1]);
        g.cell_at(cx, cy).and_then(|(ix, iy)| {
            let col = &g.cols[g.idx(ix, iy)];
            let (top, tm) = solid_top(col)?;
            let film_m = if matches!(tm, Material::Silicon | Material::Dopant(_)) { Material::Silicon } else { tm };
            let d = if film_m == Material::Silicon { 0.0 } else { top - surface };
            let (rp, rs) = film_on_silicon_reflection(film_m, d.max(0.0), p.wavelength_nm, p.incidence_deg);
            ellipso::run(&ellipso::MetrologyInputs { r_p_re: rp.0, r_p_im: rp.1, r_s_re: rs.0, r_s_im: rs.1, incidence_deg: p.incidence_deg }).ok()
                .map(|o| Ellipsometry { psi_deg: o.psi_deg, delta_deg: o.delta_deg, film_nm: d.max(0.0), film: film_m })
        })
    };

    let measured: Vec<&Measurement> = measurements.iter().filter(|m| m.status != Status::NotMeasured).collect();
    let pass = if measured.is_empty() { None } else { Some(measured.iter().all(|m| m.status == Status::Pass)) };
    Ok(MetrologyReport { wafer_fingerprint: wafer.fingerprint(), window: *window, measurements, wafer_map, optical, pass })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bare_silicon_ellipsometry_matches_the_known_brewster_region() {
        // Bare Si at 70°, 633 nm: Psi ≈ 10–12°, well away from 45°.
        let (rp, rs) = film_on_silicon_reflection(Material::Silicon, 0.0, 633.0, 70.0);
        let o = ellipso::run(&ellipso::MetrologyInputs { r_p_re: rp.0, r_p_im: rp.1, r_s_re: rs.0, r_s_im: rs.1, incidence_deg: 70.0 }).unwrap();
        assert!(o.psi_deg > 5.0 && o.psi_deg < 20.0, "psi {}", o.psi_deg);
    }

    #[test]
    fn oxide_film_changes_delta_monotonically_for_thin_films() {
        let delta = |d: f64| {
            let (rp, rs) = film_on_silicon_reflection(Material::SiliconDioxide, d, 633.0, 70.0);
            ellipso::run(&ellipso::MetrologyInputs { r_p_re: rp.0, r_p_im: rp.1, r_s_re: rs.0, r_s_im: rs.1, incidence_deg: 70.0 }).unwrap().delta_deg
        };
        assert!((delta(0.0) - delta(100.0)).abs() > 1.0, "a 100 nm oxide must shift Delta");
    }

    #[test]
    fn seeded_noise_is_deterministic() {
        assert_eq!(gaussian(7, 3), gaussian(7, 3));
        assert_ne!(gaussian(7, 3), gaussian(8, 3));
    }
}
