//! Dijkstra shortest path and distance weighting for surface code lattices.
//!
//! Computes all-pairs shortest paths and path reconstructions across the dual spacetime lattice.

use crate::code::lattice::SurfaceLattice;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// State for Dijkstra priority queue.
#[derive(Copy, Clone, Eq, PartialEq)]
struct State {
    cost: usize,
    position: usize,
}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        // Reverse ordering for min-heap
        other.cost.cmp(&self.cost)
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Dijkstra shortest path engine for surface code data qubit networks.
pub struct DijkstraSpacetime;

impl DijkstraSpacetime {
    /// Computes the shortest sequence of data qubits connecting two ancillae of the same type.
    pub fn shortest_path(lattice: &SurfaceLattice, a1: usize, a2: usize) -> Vec<usize> {
        lattice.shortest_data_qubit_path(a1, a2)
    }

    /// Computes the shortest sequence of data qubits from an ancilla to its relevant code boundary.
    pub fn shortest_boundary_path(lattice: &SurfaceLattice, a: usize) -> Vec<usize> {
        lattice.shortest_path_to_boundary(a)
    }

    /// Computes the all-pairs shortest data qubit path lengths across all data qubits.
    pub fn all_pairs_data_distances(lattice: &SurfaceLattice) -> Vec<Vec<usize>> {
        let n = lattice.num_data_qubits;
        let mut dist_matrix = vec![vec![usize::MAX; n]; n];

        for start in 0..n {
            let mut dist = vec![usize::MAX; n];
            let mut heap = BinaryHeap::new();

            dist[start] = 0;
            heap.push(State { cost: 0, position: start });

            while let Some(State { cost, position }) = heap.pop() {
                if cost > dist[position] {
                    continue;
                }

                let (r, c) = lattice.data_qubit_coord(position);
                let neighbors = [
                    (r.wrapping_sub(1), c),
                    (r + 1, c),
                    (r, c.wrapping_sub(1)),
                    (r, c + 1),
                ];

                for &(nr, nc) in &neighbors {
                    if nr < lattice.distance && nc < lattice.distance {
                        let next = lattice.data_qubit_index(nr, nc);
                        let next_cost = cost + 1;
                        if next_cost < dist[next] {
                            dist[next] = next_cost;
                            heap.push(State { cost: next_cost, position: next });
                        }
                    }
                }
            }

            dist_matrix[start] = dist;
        }

        dist_matrix
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_pairs_distances_d3() {
        let lattice = SurfaceLattice::new(3);
        let dists = DijkstraSpacetime::all_pairs_data_distances(&lattice);
        assert_eq!(dists.len(), 9);
        // Distance from qubit 0 (0,0) to qubit 8 (2,2) is |2-0| + |2-0| = 4
        assert_eq!(dists[0][8], 4);
        assert_eq!(dists[0][0], 0);
        assert_eq!(dists[0][1], 1);
    }
}
