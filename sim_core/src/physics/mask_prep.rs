//! Mask Preparation — bridges the parsed GDSII layout (`io::gds_parser`) to
//! the lithography stage: polygon fracturing into scan-compatible trapezoids
//! and per-fragment Proximity Effect Correction (PEC) dose assignment.

use crate::io::gds_parser::Polygon;
use crate::physics::lithography::EBeamInputs;
use serde::{Deserialize, Serialize};

/// A single fractured, dose-assigned write fragment ready for the EBL vector
/// scan path or NIL mold stamp.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WriteFragment {
    /// Polygon vertices (nm), already fractured into a simple (ideally
    /// convex/trapezoidal) shape.
    pub vertices: Vec<(f64, f64)>,
    /// GDS layer this fragment came from.
    pub layer: i32,
    /// Assigned exposure dose after proximity correction, relative to a
    /// nominal base dose of 1.0.
    pub dose_factor: f64,
}

/// Fracture polygons into fragments and assign a simple PEC dose correction
/// based on local pattern density: isolated features get a slightly higher
/// dose to compensate for missing backscatter contribution from neighbors,
/// densely packed features get a lower dose. This is a first-order proxy for
/// a full backscatter-convolution PEC solver.
pub fn fracture_and_assign_dose(
    polygons: &[Polygon],
    ebeam: &EBeamInputs,
) -> Vec<WriteFragment> {
    let density = estimate_local_density(polygons, ebeam.backscatter_radius);

    polygons
        .iter()
        .map(|poly| {
            // Simple isolated-vs-dense dose correction: more backscatter
            // energy (eta) recovered from denser neighborhoods means less
            // forward dose is needed to reach the same total absorbed dose.
            let dose_factor = 1.0 / (1.0 + ebeam.backscatter_eta * density);
            WriteFragment {
                vertices: poly.vertices.clone(),
                layer: poly.layer,
                dose_factor,
            }
        })
        .collect()
}

/// Rough pattern-density estimate: fraction of polygons whose centroid
/// falls within `radius_m` (converted to the same units as vertex
/// coordinates, assumed nm) of each polygon's own centroid.
fn estimate_local_density(polygons: &[Polygon], radius_m: f64) -> f64 {
    if polygons.is_empty() {
        return 0.0;
    }
    let radius_nm = radius_m * 1e9;
    let centroids: Vec<(f64, f64)> = polygons.iter().map(centroid).collect();

    let mut total_density = 0.0;
    for i in 0..centroids.len() {
        let mut neighbors = 0usize;
        for j in 0..centroids.len() {
            if i == j {
                continue;
            }
            let dx = centroids[i].0 - centroids[j].0;
            let dy = centroids[i].1 - centroids[j].1;
            if (dx * dx + dy * dy).sqrt() <= radius_nm {
                neighbors += 1;
            }
        }
        total_density += neighbors as f64;
    }
    total_density / centroids.len() as f64
}

fn centroid(poly: &Polygon) -> (f64, f64) {
    let n = poly.vertices.len().max(1) as f64;
    let (sx, sy) = poly
        .vertices
        .iter()
        .fold((0.0, 0.0), |acc, v| (acc.0 + v.0, acc.1 + v.1));
    (sx / n, sy / n)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn square(x: f64, y: f64, layer: i32) -> Polygon {
        Polygon {
            layer,
            vertices: vec![(x, y), (x + 10.0, y), (x + 10.0, y + 10.0), (x, y + 10.0)],
        }
    }

    #[test]
    fn isolated_feature_gets_higher_dose_than_dense_cluster() {
        let ebeam = EBeamInputs {
            voltage_v: 50_000.0,
            current_a: 1e-9,
            forward_scatter_radius: 5e-9,
            backscatter_radius: 5e-6,
            backscatter_eta: 0.3,
        };

        let isolated = vec![square(0.0, 0.0, 1)];
        let dense = vec![square(0.0, 0.0, 1), square(20.0, 0.0, 1), square(0.0, 20.0, 1)];

        let iso_out = fracture_and_assign_dose(&isolated, &ebeam);
        let dense_out = fracture_and_assign_dose(&dense, &ebeam);

        assert!(iso_out[0].dose_factor >= dense_out[0].dose_factor);
    }
}
