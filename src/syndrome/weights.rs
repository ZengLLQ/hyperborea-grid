//! Edge weight calculation for spacetime syndrome defect graphs.
//!
//! Converts physical Pauli error probabilities and measurement readout flip rates into
//! log-likelihood shortest path weights according to $w = -\ln(p / (1-p))$.

use crate::code::lattice::SurfaceLattice;
use crate::syndrome::extraction::DefectEvent;
use serde::{Deserialize, Serialize};

/// Metric calculator for edge weights in spacetime syndrome graphs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricWeights {
    /// Space-like edge weight for data qubit Pauli errors.
    pub space_weight: f64,
    /// Time-like edge weight for measurement readout bit-flips.
    pub time_weight: f64,
    /// Integer scaling factor for Blossom algorithm numerical stability.
    pub scale_factor: f64,
}

impl MetricWeights {
    /// Constructs metric weights from data error probability $p$ and measurement error probability $p_{\text{meas}}$.
    pub fn from_probabilities(p: f64, p_meas: f64) -> Self {
        let p_clamped = p.clamp(1e-7, 0.499);
        let p_meas_clamped = p_meas.clamp(1e-7, 0.499);

        let space_weight = ((1.0 - p_clamped) / p_clamped).ln().max(0.1);
        let time_weight = ((1.0 - p_meas_clamped) / p_meas_clamped).ln().max(0.1);

        Self {
            space_weight,
            time_weight,
            scale_factor: 1000.0,
        }
    }

    /// Constructs uniform unit weights (Manhattan distance metric).
    pub fn uniform() -> Self {
        Self {
            space_weight: 1.0,
            time_weight: 1.0,
            scale_factor: 100.0,
        }
    }

    /// Computes the spacetime distance between two defect events as a float.
    pub fn distance(&self, lattice: &SurfaceLattice, d1: &DefectEvent, d2: &DefectEvent) -> f64 {
        let d_space = lattice.spatial_distance(d1.ancilla_id, d2.ancilla_id) as f64;
        let d_time = (d1.round as isize - d2.round as isize).unsigned_abs() as f64;
        d_space * self.space_weight + d_time * self.time_weight
    }

    /// Computes the distance from a defect event to its relevant code boundary.
    pub fn boundary_distance(&self, lattice: &SurfaceLattice, d: &DefectEvent) -> f64 {
        let d_bnd = lattice.distance_to_boundary(d.ancilla_id) as f64;
        d_bnd * self.space_weight
    }

    /// Computes integer-scaled weight for Blossom algorithm.
    pub fn scaled_distance(&self, lattice: &SurfaceLattice, d1: &DefectEvent, d2: &DefectEvent) -> i64 {
        (self.distance(lattice, d1, d2) * self.scale_factor).round() as i64
    }

    /// Computes integer-scaled distance to boundary for Blossom algorithm.
    pub fn scaled_boundary_distance(&self, lattice: &SurfaceLattice, d: &DefectEvent) -> i64 {
        (self.boundary_distance(lattice, d) * self.scale_factor).round() as i64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_metric_weights_positive() {
        let weights = MetricWeights::from_probabilities(0.01, 0.01);
        assert!(weights.space_weight > 0.0);
        assert!(weights.time_weight > 0.0);
    }
}
