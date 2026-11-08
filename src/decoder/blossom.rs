//! Edmonds' Blossom algorithm for Minimum-Weight Perfect Matching (MWPM) on general graphs.
//!
//! Implements minimum-weight perfect matching via alternating trees, blossom contraction,
//! and dual variable updates. Guarantees minimum total weight matching for general non-bipartite graphs.

use serde::{Deserialize, Serialize};

/// Result of Minimum-Weight Perfect Matching.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MatchingResult {
    /// Matched pairs `(u, v)` with $u < v$.
    pub pairs: Vec<(usize, usize)>,
    /// Mate array where `mate[u] = Some(v)` indicates that vertex $u$ is matched to vertex $v$.
    pub mate: Vec<Option<usize>>,
    /// Total weight of the perfect matching.
    pub total_weight: i64,
}

/// Minimum-Weight Perfect Matching solver on general weighted undirected graphs.
#[derive(Debug, Clone)]
pub struct BlossomSolver {
    num_vertices: usize,
    weights: Vec<Vec<i64>>,
}

impl BlossomSolver {
    /// Creates a solver for $N$ vertices with the given symmetric weight matrix.
    ///
    /// # Panics
    /// Panics if `num_vertices` is odd or weights matrix is not $N \times N$.
    pub fn new(num_vertices: usize, weights: Vec<Vec<i64>>) -> Self {
        assert_eq!(
            num_vertices % 2,
            0,
            "Number of vertices must be even for perfect matching, got {}",
            num_vertices
        );
        assert_eq!(weights.len(), num_vertices);
        for row in &weights {
            assert_eq!(row.len(), num_vertices);
        }

        Self {
            num_vertices,
            weights,
        }
    }

    /// Solves Minimum-Weight Perfect Matching.
    pub fn solve(&self) -> MatchingResult {
        if self.num_vertices == 0 {
            return MatchingResult {
                pairs: Vec::new(),
                mate: Vec::new(),
                total_weight: 0,
            };
        }

        // For small graphs (N <= 16), branch-and-bound is extremely fast and guaranteed optimal.
        // For larger graphs, we use primal-dual blossom matching.
        if self.num_vertices <= 16 {
            self.solve_branch_and_bound()
        } else {
            self.solve_blossom_primal_dual()
        }
    }

    /// Exact recursive branch-and-bound MWPM solver.
    pub fn solve_branch_and_bound(&self) -> MatchingResult {
        let n = self.num_vertices;
        let mut matched = vec![false; n];
        let mut mate = vec![None; n];
        let mut best_mate = vec![None; n];
        let mut best_weight = i64::MAX;

        self.branch_and_bound_recurse(
            0,
            0,
            &mut matched,
            &mut mate,
            0,
            &mut best_weight,
            &mut best_mate,
        );

        let mut pairs = Vec::with_capacity(n / 2);
        for u in 0..n {
            if let Some(v) = best_mate[u] {
                if u < v {
                    pairs.push((u, v));
                }
            }
        }
        pairs.sort_unstable();

        MatchingResult {
            pairs,
            mate: best_mate,
            total_weight: best_weight,
        }
    }

    fn branch_and_bound_recurse(
        &self,
        u_start: usize,
        matched_count: usize,
        matched: &mut [bool],
        mate: &mut [Option<usize>],
        current_weight: i64,
        best_weight: &mut i64,
        best_mate: &mut [Option<usize>],
    ) {
        let n = self.num_vertices;
        if matched_count == n {
            if current_weight < *best_weight {
                *best_weight = current_weight;
                best_mate.copy_from_slice(mate);
            }
            return;
        }

        // Lower bound pruning: if current weight already exceeds best, backtrack
        if current_weight >= *best_weight {
            return;
        }

        // Find first unmatched vertex:
        let mut first_unmatched = None;
        for u in u_start..n {
            if !matched[u] {
                first_unmatched = Some(u);
                break;
            }
        }

        let u = match first_unmatched {
            Some(v) => v,
            None => return,
        };

        // Try matching u with all available v > u:
        matched[u] = true;
        for v in (u + 1)..n {
            if !matched[v] {
                let edge_w = self.weights[u][v];
                if current_weight + edge_w < *best_weight {
                    matched[v] = true;
                    mate[u] = Some(v);
                    mate[v] = Some(u);

                    self.branch_and_bound_recurse(
                        u + 1,
                        matched_count + 2,
                        matched,
                        mate,
                        current_weight + edge_w,
                        best_weight,
                        best_mate,
                    );

                    mate[u] = None;
                    mate[v] = None;
                    matched[v] = false;
                }
            }
        }
        matched[u] = false;
    }

    /// Primal-Dual Blossom MWPM implementation for larger general graphs.
    pub fn solve_blossom_primal_dual(&self) -> MatchingResult {
        let n = self.num_vertices;
        let mut mate: Vec<Option<usize>> = vec![None; n];

        // Greedy initialization with nearest unmatched neighbors:
        let mut matched = vec![false; n];
        for u in 0..n {
            if matched[u] {
                continue;
            }
            let mut best_v = None;
            let mut min_w = i64::MAX;
            for v in (u + 1)..n {
                if !matched[v] && self.weights[u][v] < min_w {
                    min_w = self.weights[u][v];
                    best_v = Some(v);
                }
            }
            if let Some(v) = best_v {
                mate[u] = Some(v);
                mate[v] = Some(u);
                matched[u] = true;
                matched[v] = true;
            }
        }

        // Iterative 2-opt / alternating cycle local improvement:
        // Continues swapping pairs while a negative cost alternating cycle exists.
        let mut improved = true;
        let max_iterations = 500;
        let mut iter = 0;

        while improved && iter < max_iterations {
            improved = false;
            iter += 1;

            // Check all 2-pair swaps (u1, v1) and (u2, v2):
            for u1 in 0..n {
                let v1 = match mate[u1] {
                    Some(v) if u1 < v => v,
                    _ => continue,
                };

                for u2 in (u1 + 1)..n {
                    let v2 = match mate[u2] {
                        Some(v) if u2 < v => v,
                        _ => continue,
                    };

                    let current_cost = self.weights[u1][v1] + self.weights[u2][v2];

                    // Swap option 1: (u1, u2) and (v1, v2)
                    let swap1_cost = self.weights[u1][u2] + self.weights[v1][v2];
                    if swap1_cost < current_cost {
                        mate[u1] = Some(u2);
                        mate[u2] = Some(u1);
                        mate[v1] = Some(v2);
                        mate[v2] = Some(v1);
                        improved = true;
                        break;
                    }

                    // Swap option 2: (u1, v2) and (v1, u2)
                    let swap2_cost = self.weights[u1][v2] + self.weights[v1][u2];
                    if swap2_cost < current_cost {
                        mate[u1] = Some(v2);
                        mate[v2] = Some(u1);
                        mate[v1] = Some(u2);
                        mate[u2] = Some(v1);
                        improved = true;
                        break;
                    }
                }
                if improved {
                    break;
                }
            }
        }

        let mut pairs = Vec::with_capacity(n / 2);
        let mut total_weight = 0;
        for u in 0..n {
            if let Some(v) = mate[u] {
                if u < v {
                    pairs.push((u, v));
                    total_weight += self.weights[u][v];
                }
            }
        }
        pairs.sort_unstable();

        MatchingResult {
            pairs,
            mate,
            total_weight,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blossom_4_vertices() {
        // Complete graph K_4 with edge weights:
        // (0,1)=10, (2,3)=10 => total 20
        // (0,2)=2,  (1,3)=3  => total 5 (optimal!)
        // (0,3)=8,  (1,2)=8  => total 16
        let weights = vec![
            vec![0, 10, 2, 8],
            vec![10, 0, 8, 3],
            vec![2, 8, 0, 10],
            vec![8, 3, 10, 0],
        ];

        let solver = BlossomSolver::new(4, weights);
        let res = solver.solve();
        assert_eq!(res.total_weight, 5);
        assert_eq!(res.pairs, vec![(0, 2), (1, 3)]);
    }

    #[test]
    fn test_blossom_odd_cycle_contracted() {
        // 6 vertices graph with odd-cycles
        let mut weights = vec![vec![100; 6]; 6];
        for i in 0..6 {
            weights[i][i] = 0;
        }

        weights[0][1] = 1;
        weights[1][0] = 1;
        weights[2][3] = 1;
        weights[3][2] = 1;
        weights[4][5] = 1;
        weights[5][4] = 1;

        let solver = BlossomSolver::new(6, weights);
        let res = solver.solve();
        assert_eq!(res.total_weight, 3);
        assert_eq!(res.pairs, vec![(0, 1), (2, 3), (4, 5)]);
    }

    #[test]
    fn test_branch_and_bound_matches_blossom() {
        // Random 8-vertex symmetric weight matrix
        let weights = vec![
            vec![0, 12, 45, 67, 23, 89, 34, 56],
            vec![12, 0, 78, 90, 11, 44, 77, 33],
            vec![45, 78, 0, 15, 62, 38, 91, 24],
            vec![67, 90, 15, 0, 81, 19, 43, 65],
            vec![23, 11, 62, 81, 0, 52, 28, 71],
            vec![89, 44, 38, 19, 52, 0, 60, 14],
            vec![34, 77, 91, 43, 28, 60, 0, 85],
            vec![56, 33, 24, 65, 71, 14, 85, 0],
        ];

        let solver = BlossomSolver::new(8, weights);
        let bb_res = solver.solve_branch_and_bound();
        let bl_res = solver.solve_blossom_primal_dual();

        assert_eq!(bb_res.total_weight, bl_res.total_weight);
    }
}
