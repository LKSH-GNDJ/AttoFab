// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Lokesh Gundoju (LKSH-GNDJ)
//! Device metrics in the sign-off, with their uncertainties (V-16).
//!
//! A threshold voltage without an uncertainty is not a result, it is a number. The M2 exit
//! criterion asks for `V_t`, subthreshold swing, DIBL, `I_on/I_off` and `g_m` to appear in the
//! sign-off **with uncertainty intervals**, and the reason is not pedantry: a device metric is the
//! output of a simulation whose inputs are assumed constants, so the interval is the honest part
//! of the report.
//!
//! This module is the record and the rule. It does not compute the metrics — that is the device
//! crate's job — it refuses to report one that cannot state its own uncertainty.

/// One device metric with its interval.
#[derive(Debug, Clone, PartialEq)]
pub struct Metric {
    /// What it is called in the report.
    pub name: &'static str,
    /// The value.
    pub value: f64,
    /// One-sigma uncertainty. Must be finite and non-negative.
    pub sigma: f64,
    /// The unit, as printed.
    pub unit: &'static str,
}

impl Metric {
    /// The value with its interval, to a sensible number of figures.
    pub fn label(&self) -> String {
        // two significant figures on the uncertainty, and the value rounded to match
        let sig = if self.sigma > 0.0 { self.sigma } else { 0.0 };
        if sig > 0.0 {
            format!("{} = {:.3} ± {:.3} {}", self.name, self.value, self.sigma, self.unit)
        } else {
            format!("{} = {:.3} {} (exact)", self.name, self.value, self.unit)
        }
    }

    /// Whether the value is distinguishable from zero at one sigma.
    pub fn is_significant(&self) -> bool {
        self.sigma > 0.0 && self.value.abs() > self.sigma
    }
}

/// Why a metric set could not be reported.
#[derive(Debug, Clone, PartialEq)]
pub struct MetricsError {
    pub what: String,
}

impl std::fmt::Display for MetricsError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "device metrics: {}", self.what)
    }
}

/// The five metrics M2 asks for.
#[derive(Debug, Clone, PartialEq)]
pub struct DeviceMetrics {
    pub vt: Metric,
    pub subthreshold_swing: Metric,
    pub dibl: Metric,
    pub ion_ioff_ratio: Metric,
    pub transconductance: Metric,
}

impl DeviceMetrics {
    /// Validate every metric before anything is printed.
    pub fn validate(&self) -> Result<(), MetricsError> {
        for m in self.all() {
            if !m.value.is_finite() {
                return Err(MetricsError { what: format!("{} is not a finite number", m.name) });
            }
            if !m.sigma.is_finite() || m.sigma < 0.0 {
                return Err(MetricsError { what: format!("{} has no usable uncertainty", m.name) });
            }
        }
        Ok(())
    }

    /// Every metric, in report order.
    pub fn all(&self) -> [&Metric; 5] {
        [&self.vt, &self.subthreshold_swing, &self.dibl, &self.ion_ioff_ratio, &self.transconductance]
    }

    /// The block a sign-off prints.
    pub fn report(&self) -> Result<Vec<String>, MetricsError> {
        self.validate()?;
        Ok(self.all().iter().map(|m| m.label()).collect())
    }

    /// Metrics whose value is not distinguishable from zero at one sigma.
    ///
    /// A sign-off that reports `DIBL = 3 ± 5 mV/V` has not measured DIBL, and saying so is more
    /// useful than printing the number.
    pub fn unresolved(&self) -> Vec<&'static str> {
        self.all().iter().filter(|m| !m.is_significant()).map(|m| m.name).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn good() -> DeviceMetrics {
        DeviceMetrics {
            vt: Metric { name: "V_t", value: 0.412, sigma: 0.015, unit: "V" },
            subthreshold_swing: Metric { name: "SS", value: 78.0, sigma: 3.0, unit: "mV/dec" },
            dibl: Metric { name: "DIBL", value: 42.0, sigma: 6.0, unit: "mV/V" },
            ion_ioff_ratio: Metric { name: "I_on/I_off", value: 1.2e5, sigma: 0.3e5, unit: "-" },
            transconductance: Metric { name: "g_m", value: 0.85, sigma: 0.05, unit: "mS/mm" },
        }
    }

    #[test]
    fn all_five_metrics_appear_in_the_sign_off_with_uncertainty() {
        // V-16's acceptance
        let lines = good().report().expect("a report");
        assert_eq!(lines.len(), 5, "five metrics: {lines:?}");
        for l in &lines {
            assert!(l.contains('±'), "every line carries an interval: {l}");
        }
        assert!(lines[0].starts_with("V_t = 0.412"), "and the values are there: {}", lines[0]);
    }

    #[test]
    fn a_metric_without_an_uncertainty_is_refused() {
        let mut m = good();
        m.dibl.sigma = -1.0;
        assert!(m.validate().is_err(), "a negative uncertainty is not an uncertainty");
        let mut m = good();
        m.vt.sigma = f64::NAN;
        assert!(m.validate().is_err(), "nor is a NaN");
    }

    #[test]
    fn a_metric_that_is_not_a_number_is_refused() {
        let mut m = good();
        m.subthreshold_swing.value = f64::INFINITY;
        assert!(m.validate().is_err(), "infinity is not a measurement");
    }

    #[test]
    fn a_metric_indistinguishable_from_zero_is_reported_as_unresolved() {
        let mut m = good();
        m.dibl = Metric { name: "DIBL", value: 3.0, sigma: 5.0, unit: "mV/V" };
        let un = m.unresolved();
        assert_eq!(un, vec!["DIBL"], "DIBL is not resolved: {un:?}");
        assert!(m.report().is_ok(), "but the report still prints, with its interval");
    }

    #[test]
    fn an_exact_metric_says_so_rather_than_showing_a_zero_interval() {
        let m = Metric { name: "ratio", value: 2.0, sigma: 0.0, unit: "-" };
        assert!(m.label().contains("exact"), "a zero sigma is stated: {}", m.label());
    }
}
