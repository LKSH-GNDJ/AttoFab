//! Live Analytical Plotting Pane — instantaneous mathematical curves that
//! update as table parameters change: Gaussian dopant depth profiles,
//! Stoney film-stress thresholds, and ellipsometry spectra.

use sim_core::physics::{deposition, doping, metrology};

/// A single named curve of (x, y) points, ready to hand to egui_plot's
/// `Line::new(PlotPoints::from(...))`.
#[derive(Debug, Clone)]
pub struct Curve {
    pub name: String,
    pub points: Vec<(f64, f64)>,
}

/// A horizontal/vertical threshold line to overlay on a plot (e.g. a yield
/// strength limit).
#[derive(Debug, Clone)]
pub struct ThresholdLine {
    pub name: String,
    pub value: f64,
}

#[derive(Debug, Clone)]
pub struct PlotPane {
    pub title: String,
    pub x_label: String,
    pub y_label: String,
    pub curves: Vec<Curve>,
    pub thresholds: Vec<ThresholdLine>,
}

/// Build the dopant depth-profile plot from the current doping parameters.
pub fn dopant_profile_plot(inputs: &doping::DopingInputs) -> Result<PlotPane, doping::DopingError> {
    let x_start = (inputs.projected_range_nm - 4.0 * inputs.straggle_nm).max(0.0);
    let x_end = inputs.projected_range_nm + 4.0 * inputs.straggle_nm;
    let samples = doping::sample_profile(inputs, x_start, x_end, 200)?;

    Ok(PlotPane {
        title: "Dopant Profile N(x)".into(),
        x_label: "Depth x (nm)".into(),
        y_label: "Concentration (ions/cm\u{b3})".into(),
        curves: vec![Curve { name: "N(x)".into(), points: samples }],
        thresholds: vec![],
    })
}

/// Build a Stoney stress-vs-film-thickness sweep, with the material yield
/// strength overlaid as a threshold line so the user can see at a glance
/// where delamination risk begins.
pub fn stoney_stress_sweep_plot(
    base: &deposition::DepositionInputs,
    t_film_min: f64,
    t_film_max: f64,
    n_samples: usize,
) -> Result<PlotPane, deposition::DepositionError> {
    let mut points = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let t = t_film_min + (t_film_max - t_film_min) * (i as f64 / (n_samples.max(2) - 1) as f64);
        let inputs = deposition::DepositionInputs { t_film: t, ..*base };
        let out = deposition::run(&inputs)?;
        points.push((t, out.stress_pa));
    }

    Ok(PlotPane {
        title: "Stoney Film Stress vs. Thickness".into(),
        x_label: "Film thickness (m)".into(),
        y_label: "Stress \u{3c3} (Pa)".into(),
        curves: vec![Curve { name: "\u{3c3}(t_f)".into(), points }],
        thresholds: vec![ThresholdLine { name: "Yield strength".into(), value: base.yield_strength }],
    })
}

/// Build an ellipsometry spectrum plot (Psi, Delta) swept across incidence
/// angle, holding the complex reflection coefficients fixed. A full
/// implementation would instead sweep wavelength through a dispersion model
/// per material; this exposes the same Psi/Delta output shape for the UI.
pub fn ellipsometry_spectrum_plot(
    base: &metrology::MetrologyInputs,
    angle_min_deg: f64,
    angle_max_deg: f64,
    n_samples: usize,
) -> Result<PlotPane, metrology::MetrologyError> {
    let mut psi_points = Vec::with_capacity(n_samples);
    let mut delta_points = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let angle = angle_min_deg
            + (angle_max_deg - angle_min_deg) * (i as f64 / (n_samples.max(2) - 1) as f64);
        let inputs = metrology::MetrologyInputs { incidence_deg: angle, ..*base };
        let out = metrology::run(&inputs)?;
        psi_points.push((angle, out.psi_deg));
        delta_points.push((angle, out.delta_deg));
    }

    Ok(PlotPane {
        title: "Ellipsometry Spectrum".into(),
        x_label: "Incidence angle (deg)".into(),
        y_label: "Psi / Delta (deg)".into(),
        curves: vec![
            Curve { name: "\u{3a8} (Psi)".into(), points: psi_points },
            Curve { name: "\u{394} (Delta)".into(), points: delta_points },
        ],
        thresholds: vec![],
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dopant_plot_has_samples() {
        let inputs = doping::DopingInputs::default();
        let plot = dopant_profile_plot(&inputs).unwrap();
        assert!(!plot.curves[0].points.is_empty());
    }

    #[test]
    fn stress_sweep_reflects_thickness_dependence() {
        let base = deposition::DepositionInputs::default();
        let plot = stoney_stress_sweep_plot(&base, 0.5e-6, 2e-6, 10).unwrap();
        let first = plot.curves[0].points.first().unwrap().1;
        let last = plot.curves[0].points.last().unwrap().1;
        assert!(last.abs() < first.abs());
    }
}
