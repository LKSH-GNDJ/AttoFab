// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Design of experiments, sensitivity and Monte Carlo (build PRD W-01, W-02).
//!
//! The registry in `atto_params` declares an uncertainty for every constant, so
//! a study can ask which of them actually matters. Everything here is
//! deterministic: the same seed gives the same sample, so a reported spread can
//! be reproduced by a reviewer.

pub mod fit;
pub mod identify;
pub mod surface;

use atto_dopant::Type;
use atto_device::mos::{threshold_from, MosParams};
use atto_units::{Kelvin, Nm, PerCm2};
use atto_params::{Ambient, Species};
use atto_process::anneal::{diffuse, AnnealConfig};
use atto_process::implant::Profile;
use atto_process::oxidation_nm;

/// The knobs a study turns. Every one of them is a process input, not a
/// constant: the constants are what a Monte Carlo perturbs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Knobs {
    pub implant_energy_kev: f64,
    pub implant_dose_cm2: f64,
    pub anneal_temp_c: f64,
    pub anneal_min: f64,
    pub ox_temp_c: f64,
    pub ox_min: f64,
    pub bulk_cm3: f64,
    pub qf_cm2: f64,
}

impl Default for Knobs {
    fn default() -> Self {
        Knobs {
            implant_energy_kev: 60.0,
            implant_dose_cm2: 5.0e14,
            anneal_temp_c: 1000.0,
            anneal_min: 10.0,
            ox_temp_c: 1000.0,
            ox_min: 40.0,
            bulk_cm3: 1.0e15,
            qf_cm2: 1.0e11,
        }
    }
}

/// What a run reports.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Outcomes {
    /// Junction depth from the surface, nm.
    pub xj_nm: f64,
    /// Gate oxide thickness, nm.
    pub tox_nm: f64,
    /// Threshold voltage, V.
    pub vt_v: f64,
    /// Implant dose that is electrically active after the anneal, cm^-2.
    ///
    /// Not the same as the implanted dose: everything above the solid solubility limit
    /// precipitates and donates nothing, and it is the active dose that sheet resistance and
    /// every device metric see (D-14).
    pub active_dose_cm2: f64,
}

/// Run the chain: implant, anneal, oxidise, extract.
///
/// The anneal is the real one — a 1-D diffusion of the as-implanted profile at
/// leaf resolution — so a study measures the model, not a fitted curve.
pub fn run(k: &Knobs) -> Option<Outcomes> {
    let profile = Profile::from_species(Species::Phosphorus, k.implant_energy_kev).ok()?;
    // 1 nm leaves: the junction is found by scanning the profile, so a coarse
    // leaf quantises xj and a finite difference below that resolution reads zero.
    let leaf_nm = 1.0;
    let depth_nm = profile.rp_nm + 40.0 * profile.drp_nm;
    let n = (depth_nm / leaf_nm).ceil() as usize + 2;
    let mut c: Vec<f64> = (0..n)
        .map(|i| profile.concentration_cm3(Nm::new((i as f64 + 0.5) * leaf_nm), k.implant_dose_cm2))
        .collect();
    let cfg = AnnealConfig::new(Species::Phosphorus, Type::Donor, Kelvin::from_celsius(k.anneal_temp_c), k.anneal_min);
    diffuse(&mut c, Nm::new(leaf_nm), &cfg).ok()?;
    // first depth where the profile falls below the background
    let mut xj = None;
    for i in 0..n {
        if c[i] <= k.bulk_cm3 {
            xj = Some((i as f64 + 0.5) * leaf_nm);
            break;
        }
    }
    let tox_nm = oxidation_nm(Ambient::Dry, Kelvin::from_celsius(k.ox_temp_c), k.ox_min, 0.0).ok()?;
    let vt = threshold_from(Nm::new(tox_nm), -k.bulk_cm3, PerCm2::new(k.qf_cm2), &MosParams::default())?.vt;
    // What the anneal made electrically active. The peak of the diffused profile is the
    // concentration the solubility limit acts on; the dose keeps the same active fraction.
    let peak = c.iter().cloned().fold(0.0f64, f64::max);
    let act = atto_dopant::activation::activate(
        atto_units::PerCm3(peak),
        Species::Phosphorus,
        Kelvin::from_celsius(k.anneal_temp_c),
    )
    .ok()?;
    Some(Outcomes { xj_nm: xj?, tox_nm, vt_v: vt, active_dose_cm2: k.implant_dose_cm2 * act.fraction })
}

/// One knob, by name, so a table can label a column.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Knob {
    ImplantEnergy,
    ImplantDose,
    AnnealTemp,
    AnnealMinutes,
    OxTemp,
    OxMinutes,
    Bulk,
    Qf,
}

impl Knob {
    /// Every knob a study can move.
    pub const ALL: [Knob; 8] = [
        Knob::ImplantEnergy,
        Knob::ImplantDose,
        Knob::AnnealTemp,
        Knob::AnnealMinutes,
        Knob::OxTemp,
        Knob::OxMinutes,
        Knob::Bulk,
        Knob::Qf,
    ];

    /// The name a report uses.
    pub fn label(self) -> &'static str {
        match self {
            Knob::ImplantEnergy => "implant energy (keV)",
            Knob::ImplantDose => "implant dose (cm-2)",
            Knob::AnnealTemp => "anneal temperature (C)",
            Knob::AnnealMinutes => "anneal time (min)",
            Knob::OxTemp => "oxidation temperature (C)",
            Knob::OxMinutes => "oxidation time (min)",
            Knob::Bulk => "bulk doping (cm-3)",
            Knob::Qf => "fixed charge (cm-2)",
        }
    }

    /// A knob's current value.
    pub fn get(self, k: &Knobs) -> f64 {
        match self {
            Knob::ImplantEnergy => k.implant_energy_kev,
            Knob::ImplantDose => k.implant_dose_cm2,
            Knob::AnnealTemp => k.anneal_temp_c,
            Knob::AnnealMinutes => k.anneal_min,
            Knob::OxTemp => k.ox_temp_c,
            Knob::OxMinutes => k.ox_min,
            Knob::Bulk => k.bulk_cm3,
            Knob::Qf => k.qf_cm2,
        }
    }

    /// Set a knob's value.
    pub fn set(self, k: &mut Knobs, v: f64) {
        match self {
            Knob::ImplantEnergy => k.implant_energy_kev = v,
            Knob::ImplantDose => k.implant_dose_cm2 = v,
            Knob::AnnealTemp => k.anneal_temp_c = v,
            Knob::AnnealMinutes => k.anneal_min = v,
            Knob::OxTemp => k.ox_temp_c = v,
            Knob::OxMinutes => k.ox_min = v,
            Knob::Bulk => k.bulk_cm3 = v,
            Knob::Qf => k.qf_cm2 = v,
        }
    }
}

/// A one-dimensional sweep.
#[derive(Debug, Clone)]
pub struct SweepPoint {
    pub value: f64,
    pub outcomes: Outcomes,
}

/// Run a recipe across a knob's range.
pub fn sweep(base: &Knobs, knob: Knob, values: &[f64]) -> Vec<SweepPoint> {
    let mut out = Vec::with_capacity(values.len());
    for v in values {
        let mut k = *base;
        knob.set(&mut k, *v);
        if let Some(o) = run(&k) {
            out.push(SweepPoint { value: *v, outcomes: o });
        }
    }
    out
}

/// Finite-difference sensitivity of every output to every knob, normalised as
/// a relative change: `(dy/y) / (dx/x)`.
pub fn sensitivity(base: &Knobs) -> Vec<(Knob, Outcomes)> {
    let mut out = Vec::new();
    let Some(y0) = run(base) else { return out };
    for knob in Knob::ALL {
        let x0 = knob.get(base);
        if x0 == 0.0 {
            continue;
        }
        // 1%: the outputs are quantised (xj to the leaf, tox to the lattice), so
        // a smaller step measures the grid rather than the model.
        let h = x0.abs() * 1.0e-2;
        let mut up = *base;
        knob.set(&mut up, x0 + h);
        let mut dn = *base;
        knob.set(&mut dn, x0 - h);
        let (Some(a), Some(b)) = (run(&up), run(&dn)) else { continue };
        let rel = |ya: f64, yb: f64, y: f64| -> f64 {
            if y == 0.0 {
                0.0
            } else {
                ((ya - yb) / (2.0 * h)) * (x0 / y)
            }
        };
        out.push((
            knob,
            Outcomes {
                xj_nm: rel(a.xj_nm, b.xj_nm, y0.xj_nm),
                tox_nm: rel(a.tox_nm, b.tox_nm, y0.tox_nm),
                vt_v: rel(a.vt_v, b.vt_v, y0.vt_v),
                active_dose_cm2: rel(a.active_dose_cm2, b.active_dose_cm2, y0.active_dose_cm2),
            },
        ));
    }
    out
}

/// A constant with an uncertainty, as the registry declares it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Perturbation {
    pub knob: Knob,
    /// Relative standard deviation, e.g. 0.05 for 5%.
    pub rel_sigma: f64,
    pub label: &'static str,
}

/// Deterministic sample generator: xorshift64*, so a seed reproduces a run.
pub struct Rng(u64);

impl Rng {
    /// A seeded sampler, so a run can be repeated.
    pub fn new(seed: u64) -> Self {
        Rng(if seed == 0 { 0x9E37_79B9_7F4A_7C15 } else { seed })
    }
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Standard normal by Box-Muller.
    pub fn normal(&mut self) -> f64 {
        let u1 = ((self.next_u64() >> 11) as f64 + 1.0) / (1u64 << 53) as f64;
        let u2 = ((self.next_u64() >> 11) as f64 + 1.0) / (1u64 << 53) as f64;
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// How much one knob contributes to a target's variance.
pub struct Contribution {
    pub label: &'static str,
    /// Correlation between this perturbation and the output.
    pub correlation: f64,
    /// Share of the output variance, |corr|^2 normalised.
    pub share: f64,
}

#[derive(Debug, Clone)]
/// A Monte Carlo study over declared input distributions.
pub struct MonteCarlo {
    pub samples: usize,
    pub mean: Outcomes,
    pub sd: Outcomes,
    /// Per-constant contribution, largest first.
    pub contributions: Vec<Contribution>,
}

/// Sample the declared distributions and report a yield.
pub fn monte_carlo(base: &Knobs, pert: &[Perturbation], samples: usize, seed: u64) -> Option<MonteCarlo> {
    if samples < 2 || pert.is_empty() {
        return None;
    }
    let mut rng = Rng::new(seed);
    let mut xs: Vec<Vec<f64>> = vec![Vec::with_capacity(samples); pert.len()];
    let mut ys: Vec<Outcomes> = Vec::with_capacity(samples);
    for _ in 0..samples {
        let mut k = *base;
        for (i, p) in pert.iter().enumerate() {
            let x0 = p.knob.get(base);
            let z = rng.normal();
            xs[i].push(z);
            p.knob.set(&mut k, x0 * (1.0 + p.rel_sigma * z));
        }
        ys.push(run(&k)?);
    }
    let n = samples as f64;
    let mean = Outcomes {
        xj_nm: ys.iter().map(|o| o.xj_nm).sum::<f64>() / n,
        tox_nm: ys.iter().map(|o| o.tox_nm).sum::<f64>() / n,
        vt_v: ys.iter().map(|o| o.vt_v).sum::<f64>() / n,
        active_dose_cm2: ys.iter().map(|o| o.active_dose_cm2).sum::<f64>() / n,
    };
    let sd = Outcomes {
        xj_nm: std_dev(ys.iter().map(|o| o.xj_nm), mean.xj_nm),
        tox_nm: std_dev(ys.iter().map(|o| o.tox_nm), mean.tox_nm),
        vt_v: std_dev(ys.iter().map(|o| o.vt_v), mean.vt_v),
        active_dose_cm2: std_dev(ys.iter().map(|o| o.active_dose_cm2), mean.active_dose_cm2),
    };
    let y: Vec<f64> = ys.iter().map(|o| o.vt_v).collect();
    let mut contributions = Vec::new();
    for (i, p) in pert.iter().enumerate() {
        let c = correlation(&xs[i], &y);
        contributions.push(Contribution { label: p.label, correlation: c, share: c * c });
    }
    let total: f64 = contributions.iter().map(|c| c.share).sum();
    if total > 0.0 {
        for c in contributions.iter_mut() {
            c.share /= total;
        }
    }
    contributions.sort_by(|a, b| b.share.partial_cmp(&a.share).unwrap_or(std::cmp::Ordering::Equal));
    Some(MonteCarlo { samples, mean, sd, contributions })
}

fn std_dev(it: impl Iterator<Item = f64>, mean: f64) -> f64 {
    let mut n = 0.0;
    let mut s = 0.0;
    for v in it {
        s += (v - mean) * (v - mean);
        n += 1.0;
    }
    if n < 2.0 {
        0.0
    } else {
        (s / (n - 1.0)).sqrt()
    }
}

fn correlation(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    let ma = a.iter().sum::<f64>() / n;
    let mb = b.iter().sum::<f64>() / n;
    let (mut num, mut da, mut db) = (0.0, 0.0, 0.0);
    for i in 0..a.len() {
        let x = a[i] - ma;
        let y = b[i] - mb;
        num += x * y;
        da += x * x;
        db += y * y;
    }
    if da <= 0.0 || db <= 0.0 {
        0.0
    } else {
        num / (da.sqrt() * db.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_chain_produces_plausible_numbers() {
        let o = run(&Knobs::default()).expect("the chain runs");
        assert!(o.xj_nm > 240.0 && o.xj_nm < 420.0, "xj {:.1} nm", o.xj_nm);
        assert!(o.tox_nm > 38.0 && o.tox_nm < 75.0, "tox {:.1} nm", o.tox_nm);
        assert!(o.vt_v.abs() < 2.0, "Vt {:.3} V", o.vt_v);
    }

    #[test]
    fn a_sweep_is_monotone_in_anneal_time() {
        let base = Knobs::default();
        let pts = sweep(&base, Knob::AnnealMinutes, &[1.0, 5.0, 20.0, 60.0]);
        assert_eq!(pts.len(), 4);
        for w in pts.windows(2) {
            assert!(w[1].outcomes.xj_nm >= w[0].outcomes.xj_nm, "{:?} -> {:?}", w[0].outcomes.xj_nm, w[1].outcomes.xj_nm);
        }
    }

    #[test]
    fn sensitivity_has_the_right_signs() {
        let s = sensitivity(&Knobs::default());
        let get = |k: Knob| s.iter().find(|(a, _)| *a == k).map(|(_, o)| *o).unwrap();
        assert!(get(Knob::AnnealMinutes).xj_nm > 0.0, "anneal time deepens xj");
        assert!(get(Knob::OxMinutes).tox_nm > 0.0, "oxidation time thickens tox");
        assert!(get(Knob::Bulk).vt_v > 0.0, "heavier doping raises Vt");
        assert!(get(Knob::Qf).vt_v < 0.0, "positive fixed charge lowers Vt");
    }

    #[test]
    fn the_oxide_sensitivity_of_vt_flips_with_fixed_charge() {
        // Vt = phi_ms + 2 phi_F + (Qdep - Qf)/Cox. A thicker oxide raises Vt only
        // while Qf is small: here Qf = 1e11 cm^-2 is comparable to Qdep, so the
        // sign inverts. The sweep found this; the test pins it down.
        let vt_of = |qf: f64| {
            let mut k = Knobs::default();
            k.qf_cm2 = qf;
            sensitivity(&k).iter().find(|(a, _)| *a == Knob::OxMinutes).map(|(_, o)| o.vt_v).unwrap()
        };
        assert!(vt_of(0.0) > 0.0, "no fixed charge: thicker oxide raises Vt");
        assert!(vt_of(1.0e11) < 0.0, "Qf comparable to Qdep: the sign inverts");
        // and the crossover is where Qf = Qdep
        let base = run(&Knobs::default()).unwrap();
        let qdep_sheet = {
            let cox = atto_device::mos::EPS_OX_R * atto_device::mos::EPS0 / (base.tox_nm * 1e-7);
            let phi_f = 0.02585 * (1.0e15f64 / 1.0e10).ln();
            let qdep = (2.0 * atto_device::mos::EPS_SI_R * atto_device::mos::EPS0 * atto_device::mos::QE * 1.0e15 * 2.0 * phi_f).sqrt();
            let _ = cox;
            qdep / atto_device::mos::QE
        };
        assert!(vt_of(qdep_sheet * 0.9) > 0.0);
        assert!(vt_of(qdep_sheet * 1.1) < 0.0);
    }

    #[test]
    fn monte_carlo_is_reproducible_and_centred() {
        let base = Knobs::default();
        let pert = [
            Perturbation { knob: Knob::AnnealTemp, rel_sigma: 0.02, label: "anneal T" },
            Perturbation { knob: Knob::OxMinutes, rel_sigma: 0.05, label: "oxidation time" },
        ];
        let a = monte_carlo(&base, &pert, 200, 42).unwrap();
        let b = monte_carlo(&base, &pert, 200, 42).unwrap();
        assert_eq!(a.sd.vt_v.to_bits(), b.sd.vt_v.to_bits(), "same seed, same spread");
        assert!(a.sd.vt_v > 0.0);
        let base_vt = run(&base).unwrap().vt_v;
        assert!((a.mean.vt_v - base_vt).abs() < 0.2 * a.sd.vt_v.max(1e-9) + 0.05, "mean {:.3} vs base {:.3}", a.mean.vt_v, base_vt);
    }

    #[test]
    fn zero_uncertainty_gives_zero_spread() {
        let base = Knobs::default();
        let pert = [Perturbation { knob: Knob::AnnealTemp, rel_sigma: 0.0, label: "anneal T" }];
        let mc = monte_carlo(&base, &pert, 50, 7).unwrap();
        assert!(mc.sd.vt_v < 1e-12, "sd {:.3e}", mc.sd.vt_v);
    }

    #[test]
    fn the_dominant_contributor_is_the_one_that_was_perturbed_hardest() {
        let base = Knobs::default();
        let pert = [
            Perturbation { knob: Knob::OxMinutes, rel_sigma: 0.20, label: "oxidation time (wide)" },
            Perturbation { knob: Knob::AnnealTemp, rel_sigma: 0.001, label: "anneal T (narrow)" },
        ];
        let mc = monte_carlo(&base, &pert, 400, 3).unwrap();
        assert_eq!(mc.contributions[0].label, "oxidation time (wide)", "{:?}", mc.contributions);
        assert!(mc.contributions[0].share > 0.8, "share {:.2}", mc.contributions[0].share);
    }

    #[test]
    fn the_rng_is_deterministic_and_roughly_normal() {
        let mut a = Rng::new(1);
        let mut b = Rng::new(1);
        for _ in 0..10 {
            assert_eq!(a.normal().to_bits(), b.normal().to_bits());
        }
        let mut r = Rng::new(99);
        let xs: Vec<f64> = (0..2000).map(|_| r.normal()).collect();
        let m = xs.iter().sum::<f64>() / xs.len() as f64;
        let s = (xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (xs.len() as f64 - 1.0)).sqrt();
        assert!(m.abs() < 0.1, "mean {m:.3}");
        assert!((s - 1.0).abs() < 0.1, "sd {s:.3}");
    }
}

#[cfg(test)]
mod active_dose_tests {
    use super::*;

    #[test]
    fn a_heavy_implant_reports_less_active_dose_than_it_put_in() {
        // D-14 in the study chain: the dose that survives the anneal is below the implant dose
        let k = Knobs { implant_dose_cm2: 5.0e15, ..Default::default() };
        let o = run(&k).expect("the chain runs");
        assert!(o.active_dose_cm2 > 0.0, "some of it is active: {}", o.active_dose_cm2);
        assert!(o.active_dose_cm2 < k.implant_dose_cm2, "and not all of it: {} vs {}", o.active_dose_cm2, k.implant_dose_cm2);
    }

    #[test]
    fn the_active_dose_follows_the_implant_dose() {
        let a = run(&Knobs { implant_dose_cm2: 1.0e14, ..Default::default() }).expect("a");
        let b = run(&Knobs { implant_dose_cm2: 2.0e14, ..Default::default() }).expect("b");
        assert!(b.active_dose_cm2 > a.active_dose_cm2, "more in, more active: {} vs {}", b.active_dose_cm2, a.active_dose_cm2);
    }
}
