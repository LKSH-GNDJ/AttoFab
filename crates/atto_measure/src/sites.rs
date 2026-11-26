// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Site sampling and measurement noise (build PRD M-03, audit F-19).
//!
//! the earlier build added noise to the reported mean and left the site values untouched, so
//! the number moved and its spread did not — which makes a uniformity figure
//! meaningless and hides the fact that a mean of many sites is more certain than
//! any one of them. Here every site carries its own noise, and the uncertainty
//! of the mean falls as sigma/sqrt(n) because that is what averaging does.

/// Where the sites sit on the wafer.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SitePlan {
    /// One site: the measurement is what it is.
    Single,
    /// Five sites: centre plus four at 45 degrees.
    Cross5,
    /// Nine sites: a 3x3 grid.
    Grid9,
    /// Nine sites on a circle plus the centre: the standard radial plan for a
    /// rotating-stage instrument (PRD M-05).
    Radial9,
    /// Five sites along the two axes: the cross plan, which catches astigmatism
    /// that a grid averages away.
    CrossAxis5,
    /// A declared number of pseudo-random sites, reproducible from a seed.
    Random(u64, usize),
}

impl SitePlan {
    /// Site offsets in units of the site radius, y up.
    ///
    /// Random sites are reproducible from the seed: a plan that moves between runs
    /// makes every measurement of it unrepeatable.
    pub fn offsets(self) -> Vec<(f64, f64)> {
        match self {
            SitePlan::Single => vec![(0.0, 0.0)],
            SitePlan::Cross5 => vec![(0.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)],
            SitePlan::Grid9 => {
                let mut v = Vec::with_capacity(9);
                for j in -1..=1 {
                    for i in -1..=1 {
                        v.push((i as f64, j as f64));
                    }
                }
                v
            }
            SitePlan::Radial9 => {
                let mut v = vec![(0.0, 0.0)];
                for k in 0..8 {
                    let a = std::f64::consts::TAU * k as f64 / 8.0;
                    v.push((a.cos(), a.sin()));
                }
                v
            }
            SitePlan::CrossAxis5 => vec![(0.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)],
            SitePlan::Random(seed, n) => {
                let mut state = seed | 1;
                let mut v = Vec::with_capacity(n);
                for _ in 0..n {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    let x = ((state & 0xffff) as f64 / 65535.0) * 2.0 - 1.0;
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    let y = ((state & 0xffff) as f64 / 65535.0) * 2.0 - 1.0;
                    v.push((x, y));
                }
                v
            }
        }
    }

    /// How many sites the plan places.
    pub fn count(self) -> usize {
        match self {
            SitePlan::Single => 1,
            SitePlan::Cross5 => 5,
            SitePlan::Grid9 => 9,
            SitePlan::Radial9 => 9,
            SitePlan::CrossAxis5 => 5,
            SitePlan::Random(_, n) => n,
        }
    }

    /// The name used in a report.
    pub fn label(self) -> &'static str {
        match self {
            SitePlan::Single => "1 site",
            SitePlan::Cross5 => "5-site cross",
            SitePlan::Grid9 => "9-site grid",
            SitePlan::Radial9 => "9-site radial",
            SitePlan::CrossAxis5 => "5-site axis cross",
            // the count and seed are carried by the plan, not the label
            SitePlan::Random(_, _) => "random",
        }
    }
}

/// A small deterministic generator, so a reported spread can be reproduced.
#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    /// A noise source with a reproducible seed.
    pub fn new(seed: u64) -> Self {
        Rng(if seed == 0 { 0x9E37_79B9_7F4A_7C15 } else { seed })
    }
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// Standard normal by Box-Muller.
    pub fn normal(&mut self) -> f64 {
        let u1 = ((self.next() >> 11) as f64 + 1.0) / (1u64 << 53) as f64;
        let u2 = ((self.next() >> 11) as f64 + 1.0) / (1u64 << 53) as f64;
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

/// What a multi-site measurement reports.
#[derive(Debug, Clone, PartialEq)]
pub struct SiteResult {
    pub values: Vec<f64>,
    pub mean: f64,
    /// Standard deviation of the site values themselves.
    pub sigma: f64,
    /// Standard error of the mean, `sigma / sqrt(n)`.
    pub sigma_mean: f64,
    /// `(max - min) / (max + min)`, as a percentage.
    pub uniformity_pct: f64,
}

/// Measure one quantity at every site, each with its own noise.
pub fn sample(value: f64, sigma: f64, plan: SitePlan, seed: u64) -> SiteResult {
    let n = plan.count();
    let mut rng = Rng::new(seed);
    let mut values = Vec::with_capacity(n);
    for _ in 0..n {
        let z = if sigma > 0.0 { rng.normal() } else { 0.0 };
        values.push(value + sigma * z);
    }
    summarise(values)
}

/// Summarise site values. A single site has no spread and no standard error.
pub fn summarise(values: Vec<f64>) -> SiteResult {
    let n = values.len();
    if n == 0 {
        return SiteResult { values, mean: 0.0, sigma: 0.0, sigma_mean: 0.0, uniformity_pct: 0.0 };
    }
    let mean = values.iter().sum::<f64>() / n as f64;
    let sigma = if n < 2 {
        0.0
    } else {
        (values.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / (n as f64 - 1.0)).sqrt()
    };
    let sigma_mean = if n < 2 { 0.0 } else { sigma / (n as f64).sqrt() };
    let mx = values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mn = values.iter().cloned().fold(f64::INFINITY, f64::min);
    let uniformity_pct = if (mx + mn).abs() > 0.0 { 100.0 * (mx - mn) / (mx + mn).abs() } else { 0.0 };
    SiteResult { values, mean, sigma, sigma_mean, uniformity_pct }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_site_carries_its_own_noise() {
        let r = sample(100.0, 1.0, SitePlan::Grid9, 42);
        assert_eq!(r.values.len(), 9);
        let distinct = r.values.windows(2).filter(|w| (w[1] - w[0]).abs() > 1e-9).count();
        assert!(distinct >= 7, "the sites must differ from each other, got {:?}", r.values);
        assert!(r.sigma > 0.0, "and the spread must be real");
    }

    #[test]
    fn the_radial_plan_puts_eight_sites_on_a_circle() {
        let o = SitePlan::Radial9.offsets();
        assert_eq!(o.len(), 9);
        assert_eq!(o[0], (0.0, 0.0), "the centre first");
        for p in &o[1..] {
            let r = (p.0 * p.0 + p.1 * p.1).sqrt();
            assert!((r - 1.0).abs() < 1e-12, "every ring site is at radius one: {r:.6}");
        }
    }

    #[test]
    fn the_cross_plan_covers_both_axes() {
        let o = SitePlan::CrossAxis5.offsets();
        assert_eq!(o.len(), 5);
        assert!(o.contains(&(1.0, 0.0)) && o.contains(&(-1.0, 0.0)));
        assert!(o.contains(&(0.0, 1.0)) && o.contains(&(0.0, -1.0)));
    }

    #[test]
    fn the_random_plan_is_reproducible_from_its_seed() {
        let a = SitePlan::Random(0x1234, 20).offsets();
        let b = SitePlan::Random(0x1234, 20).offsets();
        let c = SitePlan::Random(0x9999, 20).offsets();
        assert_eq!(a, b, "the same seed gives the same sites");
        assert_ne!(a, c, "a different seed gives different sites");
        assert_eq!(a.len(), 20);
        assert!(a.iter().all(|(x, y)| x.abs() <= 1.0 && y.abs() <= 1.0), "inside the unit square");
    }

    #[test]
    fn every_plan_reports_its_own_count() {
        for p in [SitePlan::Single, SitePlan::Cross5, SitePlan::Grid9, SitePlan::Radial9, SitePlan::CrossAxis5, SitePlan::Random(7, 13)] {
            assert_eq!(p.count(), p.offsets().len(), "{:?} disagrees with itself", p.label());
            assert!(!p.label().is_empty());
        }
    }

    #[test]
    fn the_mean_of_many_sites_is_more_certain_than_one() {
        let one = sample(100.0, 1.0, SitePlan::Single, 7);
        let nine = sample(100.0, 1.0, SitePlan::Grid9, 7);
        assert_eq!(one.sigma_mean, 0.0, "a single site has no standard error");
        assert!(nine.sigma_mean > 0.0);
        // sigma/sqrt(9) = sigma/3, to within sampling error
        let ratio = nine.sigma_mean / (nine.sigma / 3.0);
        assert!((ratio - 1.0).abs() < 1e-9, "ratio {ratio}");
    }

    #[test]
    fn the_same_seed_gives_the_same_sites() {
        let a = sample(100.0, 1.0, SitePlan::Cross5, 99);
        let b = sample(100.0, 1.0, SitePlan::Cross5, 99);
        assert_eq!(a, b);
        let c = sample(100.0, 1.0, SitePlan::Cross5, 100);
        assert_ne!(a.values, c.values, "a different seed gives different sites");
    }

    #[test]
    fn zero_noise_gives_identical_sites_and_zero_uniformity() {
        let r = sample(100.0, 0.0, SitePlan::Grid9, 1);
        assert!(r.values.iter().all(|v| (v - 100.0).abs() < 1e-12));
        assert_eq!(r.sigma, 0.0);
        assert_eq!(r.uniformity_pct, 0.0);
    }

    #[test]
    fn uniformity_is_the_range_over_the_mean() {
        let r = summarise(vec![90.0, 100.0, 110.0]);
        assert!((r.mean - 100.0).abs() < 1e-12);
        // the industry convention is (max - min) / (max + min): 20/200 = 10 %,
        // not range over mean, which would give 20 %
        assert!((r.uniformity_pct - 10.0).abs() < 1e-9, "uniformity {:.3}", r.uniformity_pct);
    }

    #[test]
    fn an_empty_plan_does_not_divide_by_zero() {
        let r = summarise(vec![]);
        assert_eq!(r.mean, 0.0);
        assert_eq!(r.sigma, 0.0);
    }
}
