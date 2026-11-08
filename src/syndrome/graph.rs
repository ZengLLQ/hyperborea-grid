//! Spacetime syndrome defect graph construction with virtual boundary pairing.
//!
//! Builds complete weighted matching graphs from detected defect events, incorporating
//! spacetime metrics and virtual boundary vertices for minimum-weight perfect matching.

use crate::code::lattice::SurfaceLattice;
use crate::syndrome::extraction::DefectEvent;
use crate::syndrome::weights::MetricWeights;
use serde::{Deserialize, Serialize};

/// Spacetime syndrome defect matching graph.
///
/// Contains $2m$ vertices where $m$ is the number of detected defects:
/// - Vertices $0 \le i < m$: Real defect events $D_i$.
/// - Vertices $m \le i < 2m$: Virtual boundary companion nodes $B_{i - m}$.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpacetimeDefectGraph {
    /// List of detected defect events ($m$ items).
    pub defects: Vec<DefectEvent>,
    /// Metric calculator for weights.
    pub weights: MetricWeights,
    /// Total vertices in augmented graph ($2m$).
    pub total_vertices: usize,
    /// Pairwise weight matrix of size $(2m) \times (2m)$.
    pub weight_matrix: Vec<Vec<i64>>,
}

impl SpacetimeDefectGraph {
    /// Constructs an augmented spacetime defect matching graph.
    pub fn build(
        lattice: &SurfaceLattice,
        defects: Vec<DefectEvent>,
        weights: MetricWeights,
    ) -> Self {
        let m = defects.len();
        let total_vertices = 2 * m;

        // Initialize weight matrix with a very large sentinel value:
        let infinity = 1_000_000_000i64;
        let mut weight_matrix = vec![vec![infinity; total_vertices]; total_vertices];

        for i in 0..total_vertices {
            weight_matrix[i][i] = 0;
        }

        // 1. Defect-to-defect edges (D_i, D_j):
        for i in 0..m {
            for j in (i + 1)..m {
                let w = weights.scaled_distance(lattice, &defects[i], &defects[j]);
                weight_matrix[i][j] = w;
                weight_matrix[j][i] = w;
            }
        }

        // 2. Defect-to-boundary edges (D_i, B_i):
        for i in 0..m {
            let b_idx = m + i;
            let w_bnd = weights.scaled_boundary_distance(lattice, &defects[i]);
            weight_matrix[i][b_idx] = w_bnd;
            weight_matrix[b_idx][i] = w_bnd;
        }

        // 3. Boundary-to-boundary edges (B_i, B_j) with zero weight:
        for i in 0..m {
            let bi = m + i;
            for j in (i + 1)..m {
                let bj = m + j;
                weight_matrix[bi][bj] = 0;
                weight_matrix[bj][bi] = 0;
            }
        }

        Self {
            defects,
            weights,
            total_vertices,
            weight_matrix,
        }
    }

    /// Number of detected defect events $m$.
    #[inline]
    pub fn num_defects(&self) -> usize {
        self.defects.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::lattice::{StabilizerType, SurfaceLattice};

    #[test]
    fn test_spacetime_defect_graph_construction() {
        let lattice = SurfaceLattice::new(3);
        let weights = MetricWeights::uniform();

        let z0 = lattice.z_ancillae[0];
        let z1 = lattice.z_ancillae[1];
        let d1 = DefectEvent {
            id: 0,
            ancilla_id: z0,
            round: 0,
            kind: StabilizerType::Z,
            coord: lattice.ancillae[z0].coord,
        };
        let d2 = DefectEvent {
            id: 1,
            ancilla_id: z1,
            round: 0,
            kind: StabilizerType::Z,
            coord: lattice.ancillae[z1].coord,
        };

        let graph = SpacetimeDefectGraph::build(&lattice, vec![d1, d2], weights);
        assert_eq!(graph.num_defects(), 2);
        assert_eq!(graph.total_vertices, 4);

        // Check boundary-to-boundary weight is 0:
        assert_eq!(graph.weight_matrix[2][3], 0);
        assert_eq!(graph.weight_matrix[3][2], 0);

        // Check defect-to-defect weight is positive:
        assert!(graph.weight_matrix[0][1] > 0);
    }
}
