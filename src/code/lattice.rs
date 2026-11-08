//! Rotated planar surface code geometry and lattice coordinates.
//!
//! Implements distance-$d$ rotated surface code geometry with $N = d^2$ data qubits,
//! $(d^2 - 1)$ stabilizer ancillae ($X$ and $Z$), dual graph shortest path routing,
//! and exact boundary homology assignments.

use serde::{Deserialize, Serialize};
use std::collections::{HashMap, VecDeque};

/// Stabilizer generator Pauli type.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum StabilizerType {
    /// X-type stabilizer (checks phase-flip errors Z)
    X,
    /// Z-type stabilizer (checks bit-flip errors X)
    Z,
}

impl StabilizerType {
    /// Opposite stabilizer type ($X \leftrightarrow Z$).
    pub fn opposite(&self) -> Self {
        match self {
            StabilizerType::X => StabilizerType::Z,
            StabilizerType::Z => StabilizerType::X,
        }
    }
}

/// Boundary orientation of the rotated surface code.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Boundary {
    /// Top boundary (Z-type smooth boundary)
    Top,
    /// Bottom boundary (Z-type smooth boundary)
    Bottom,
    /// Left boundary (X-type rough boundary)
    Left,
    /// Right boundary (X-type rough boundary)
    Right,
}

/// Metadata and support for a single stabilizer ancilla in the rotated surface code.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AncillaDescriptor {
    /// Ancilla index $0 \le a < d^2 - 1$.
    pub id: usize,
    /// Stabilizer type ($X$ or $Z$).
    pub kind: StabilizerType,
    /// Spatial lattice coordinate $(u, v)$ in embedding grid.
    pub coord: (i32, i32),
    /// Data qubit indices contained in this stabilizer's support.
    pub support: Vec<usize>,
    /// Whether this is a bulk stabilizer (weight 4) or boundary stabilizer (weight 2).
    pub is_bulk: bool,
    /// If boundary stabilizer, which boundary it touches.
    pub boundary: Option<Boundary>,
}

/// Dual graph representation for exact shortest path routing and boundary matching.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DualGraph {
    /// Adjacency list: node -> list of (neighbor_node, data_qubit_edge).
    /// Node index M is the virtual boundary node BND.
    pub adj: Vec<Vec<(usize, usize)>>,
    /// Number of stabilizer nodes.
    pub num_nodes: usize,
    /// Index of boundary node.
    pub bnd_node: usize,
}

impl DualGraph {
    /// Builds dual graph for a given stabilizer type (X or Z).
    pub fn build(
        ancillae: &[AncillaDescriptor],
        kind: StabilizerType,
        num_data_qubits: usize,
    ) -> Self {
        let num_nodes = ancillae.len();
        let bnd_node = num_nodes;
        let mut adj = vec![Vec::new(); num_nodes + 1];

        // For each data qubit, find which ancillae of this kind contain it:
        for q in 0..num_data_qubits {
            let mut touching = Vec::new();
            for anc in ancillae {
                if anc.kind == kind && anc.support.contains(&q) {
                    touching.push(anc.id);
                }
            }

            match touching.len() {
                2 => {
                    let u = touching[0];
                    let v = touching[1];
                    adj[u].push((v, q));
                    adj[v].push((u, q));
                }
                1 => {
                    let u = touching[0];
                    adj[u].push((bnd_node, q));
                    adj[bnd_node].push((u, q));
                }
                _ => {}
            }
        }

        Self {
            adj,
            num_nodes,
            bnd_node,
        }
    }

    /// Finds shortest path of data qubits between two ancillae.
    pub fn shortest_path(&self, start: usize, target: usize) -> Vec<usize> {
        if start == target {
            return Vec::new();
        }

        let total_v = self.bnd_node + 1;
        let mut queue = VecDeque::new();
        let mut visited = vec![false; total_v];
        let mut parent_node = vec![None; total_v];
        let mut parent_edge = vec![None; total_v];

        visited[start] = true;
        queue.push_back(start);

        let mut found = false;
        while let Some(curr) = queue.pop_front() {
            if curr == target {
                found = true;
                break;
            }

            for &(nxt, edge_q) in &self.adj[curr] {
                if !visited[nxt] {
                    visited[nxt] = true;
                    parent_node[nxt] = Some(curr);
                    parent_edge[nxt] = Some(edge_q);
                    queue.push_back(nxt);
                }
            }
        }

        if found {
            let mut path = Vec::new();
            let mut curr = target;
            while let Some(p) = parent_node[curr] {
                path.push(parent_edge[curr].unwrap());
                curr = p;
            }
            path.reverse();
            path
        } else {
            Vec::new()
        }
    }

    /// Finds shortest path of data qubits from an ancilla to the boundary.
    pub fn shortest_boundary_path(&self, start: usize) -> Vec<usize> {
        self.shortest_path(start, self.bnd_node)
    }
}

/// Rotated planar surface code lattice geometry for code distance $d$.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SurfaceLattice {
    /// Code distance $d$ (odd integer $\ge 3$).
    pub distance: usize,
    /// Total data qubits $N = d^2$.
    pub num_data_qubits: usize,
    /// Ancilla descriptors indexed by ancilla id $0 \le a < d^2 - 1$.
    pub ancillae: Vec<AncillaDescriptor>,
    /// Map from lattice coordinate $(u, v)$ to ancilla id.
    pub coord_to_ancilla: HashMap<(i32, i32), usize>,
    /// Indices of X-stabilizers.
    pub x_ancillae: Vec<usize>,
    /// Indices of Z-stabilizers.
    pub z_ancillae: Vec<usize>,
    /// Precomputed map of which ancillae touch each data qubit.
    pub data_to_ancillae: Vec<Vec<usize>>,
    /// Dual graph for X-stabilizers.
    pub dual_x: DualGraph,
    /// Dual graph for Z-stabilizers.
    pub dual_z: DualGraph,
}

impl SurfaceLattice {
    /// Constructs a rotated planar surface code lattice for code distance $d$.
    pub fn new(distance: usize) -> Self {
        assert!(
            distance >= 3 && distance % 2 == 1,
            "Rotated surface code distance d must be an odd integer >= 3, got {}",
            distance
        );

        let num_data_qubits = distance * distance;
        let d_i32 = distance as i32;

        let mut ancillae = Vec::new();
        let mut coord_to_ancilla = HashMap::new();
        let mut x_ancillae = Vec::new();
        let mut z_ancillae = Vec::new();

        let mut ancilla_id = 0;

        for u in (-1..=(2 * d_i32 - 1)).step_by(2) {
            for v in (-1..=(2 * d_i32 - 1)).step_by(2) {
                let is_bulk = u >= 1 && u <= 2 * d_i32 - 3 && v >= 1 && v <= 2 * d_i32 - 3;
                let is_top = u == -1 && v >= 1 && v <= 2 * d_i32 - 3;
                let is_bottom = u == 2 * d_i32 - 1 && v >= 1 && v <= 2 * d_i32 - 3;
                let is_left = v == -1 && u >= 1 && u <= 2 * d_i32 - 3;
                let is_right = v == 2 * d_i32 - 1 && u >= 1 && u <= 2 * d_i32 - 3;

                if !is_bulk && !is_top && !is_bottom && !is_left && !is_right {
                    continue;
                }

                let kind_opt = if is_bulk {
                    if (u + v) / 2 % 2 == 0 {
                        Some(StabilizerType::Z)
                    } else {
                        Some(StabilizerType::X)
                    }
                } else if is_top {
                    if (v - 1) / 2 % 2 == 0 {
                        Some(StabilizerType::Z)
                    } else {
                        None
                    }
                } else if is_bottom {
                    if (2 * d_i32 - 1 + v) / 2 % 2 == 0 {
                        Some(StabilizerType::Z)
                    } else {
                        None
                    }
                } else if is_left {
                    if (u - 1) / 2 % 2 != 0 {
                        Some(StabilizerType::X)
                    } else {
                        None
                    }
                } else if is_right {
                    if (u + 2 * d_i32 - 1) / 2 % 2 != 0 {
                        Some(StabilizerType::X)
                    } else {
                        None
                    }
                } else {
                    None
                };

                if let Some(kind) = kind_opt {
                    let candidates = [
                        (u - 1, v - 1),
                        (u - 1, v + 1),
                        (u + 1, v - 1),
                        (u + 1, v + 1),
                    ];

                    let mut support = Vec::new();
                    for &(cr, cc) in &candidates {
                        if cr >= 0 && cr <= 2 * (d_i32 - 1) && cc >= 0 && cc <= 2 * (d_i32 - 1) {
                            let r = (cr / 2) as usize;
                            let c = (cc / 2) as usize;
                            let q_idx = r * distance + c;
                            support.push(q_idx);
                        }
                    }
                    support.sort_unstable();

                    let boundary = if is_top {
                        Some(Boundary::Top)
                    } else if is_bottom {
                        Some(Boundary::Bottom)
                    } else if is_left {
                        Some(Boundary::Left)
                    } else if is_right {
                        Some(Boundary::Right)
                    } else {
                        None
                    };

                    let descriptor = AncillaDescriptor {
                        id: ancilla_id,
                        kind,
                        coord: (u, v),
                        support,
                        is_bulk,
                        boundary,
                    };

                    match kind {
                        StabilizerType::X => x_ancillae.push(ancilla_id),
                        StabilizerType::Z => z_ancillae.push(ancilla_id),
                    }

                    coord_to_ancilla.insert((u, v), ancilla_id);
                    ancillae.push(descriptor);
                    ancilla_id += 1;
                }
            }
        }

        assert_eq!(ancillae.len(), num_data_qubits - 1);

        let mut data_to_ancillae = vec![Vec::new(); num_data_qubits];
        for anc in &ancillae {
            for &q in &anc.support {
                data_to_ancillae[q].push(anc.id);
            }
        }

        let dual_x = DualGraph::build(&ancillae, StabilizerType::X, num_data_qubits);
        let dual_z = DualGraph::build(&ancillae, StabilizerType::Z, num_data_qubits);

        Self {
            distance,
            num_data_qubits,
            ancillae,
            coord_to_ancilla,
            x_ancillae,
            z_ancillae,
            data_to_ancillae,
            dual_x,
            dual_z,
        }
    }

    /// Converts data qubit index $q \in \{0, \dots, d^2 - 1\}$ to grid row and column $(r, c)$.
    #[inline]
    pub fn data_qubit_coord(&self, idx: usize) -> (usize, usize) {
        assert!(idx < self.num_data_qubits);
        (idx / self.distance, idx % self.distance)
    }

    /// Converts grid row and column $(r, c)$ to data qubit index $q$.
    #[inline]
    pub fn data_qubit_index(&self, r: usize, c: usize) -> usize {
        assert!(r < self.distance && c < self.distance);
        r * self.distance + c
    }

    /// Minimum spatial distance (in data-qubit hops) from an ancilla to its relevant code boundary.
    pub fn distance_to_boundary(&self, ancilla_id: usize) -> usize {
        let anc = &self.ancillae[ancilla_id];
        match anc.kind {
            StabilizerType::X => self.dual_x.shortest_boundary_path(ancilla_id).len().max(1),
            StabilizerType::Z => self.dual_z.shortest_boundary_path(ancilla_id).len().max(1),
        }
    }

    /// Spatial shortest path hop distance between two ancillae of the same type.
    pub fn spatial_distance(&self, a1: usize, a2: usize) -> usize {
        if a1 == a2 {
            return 0;
        }
        let kind = self.ancillae[a1].kind;
        match kind {
            StabilizerType::X => self.dual_x.shortest_path(a1, a2).len(),
            StabilizerType::Z => self.dual_z.shortest_path(a1, a2).len(),
        }
    }

    /// Finds the shortest chain of data qubits connecting two ancillae of the same type.
    pub fn shortest_data_qubit_path(&self, a1: usize, a2: usize) -> Vec<usize> {
        let kind = self.ancillae[a1].kind;
        match kind {
            StabilizerType::X => self.dual_x.shortest_path(a1, a2),
            StabilizerType::Z => self.dual_z.shortest_path(a1, a2),
        }
    }

    /// Finds the shortest chain of data qubits from an ancilla to its nearest boundary.
    pub fn shortest_path_to_boundary(&self, a: usize) -> Vec<usize> {
        let kind = self.ancillae[a].kind;
        match kind {
            StabilizerType::X => self.dual_x.shortest_boundary_path(a),
            StabilizerType::Z => self.dual_z.shortest_boundary_path(a),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_distance_3_lattice_geometry() {
        let lattice = SurfaceLattice::new(3);
        assert_eq!(lattice.num_data_qubits, 9);
        assert_eq!(lattice.ancillae.len(), 8);
        assert_eq!(lattice.x_ancillae.len(), 4);
        assert_eq!(lattice.z_ancillae.len(), 4);

        let bulk_count = lattice.ancillae.iter().filter(|a| a.is_bulk).count();
        let bnd_count = lattice.ancillae.iter().filter(|a| !a.is_bulk).count();
        assert_eq!(bulk_count, 4);
        assert_eq!(bnd_count, 4);
    }
}
