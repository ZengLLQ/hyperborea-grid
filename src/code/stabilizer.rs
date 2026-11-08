//! Stabilizer generators and stabilizer group verification.
//!
//! Constructs the $d^2 - 1$ stabilizer generators ($X$ and $Z$), verifies their mutual
//! commutation relations $[S_i, S_j] = 0$, and computes syndrome outcomes for given data qubit errors.

use super::lattice::{StabilizerType, SurfaceLattice};
use super::pauli::{Pauli, PauliString};
use serde::{Deserialize, Serialize};

/// A stabilizer generator of the surface code.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StabilizerGenerator {
    /// Generator identifier $0 \le i < d^2 - 1$.
    pub id: usize,
    /// Stabilizer type ($X$ or $Z$).
    pub kind: StabilizerType,
    /// Data qubits participating in this stabilizer.
    pub support: Vec<usize>,
    /// Full $N$-qubit Pauli string representation.
    pub pauli_string: PauliString,
    /// Whether this is a weight-4 bulk stabilizer or weight-2 boundary stabilizer.
    pub is_bulk: bool,
    /// Weight (number of non-identity Pauli factors, either 2 or 4).
    pub weight: usize,
}

/// The rotated surface code stabilizer structure.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StabilizerCode {
    /// Associated surface lattice.
    pub lattice: SurfaceLattice,
    /// List of all $d^2 - 1$ stabilizer generators.
    pub generators: Vec<StabilizerGenerator>,
}

impl StabilizerCode {
    /// Constructs and validates the stabilizer generators for a surface code lattice.
    pub fn new(lattice: SurfaceLattice) -> Self {
        let n = lattice.num_data_qubits;
        let mut generators = Vec::with_capacity(lattice.ancillae.len());

        for anc in &lattice.ancillae {
            let pauli_type = match anc.kind {
                StabilizerType::X => Pauli::X,
                StabilizerType::Z => Pauli::Z,
            };

            let pauli_string = PauliString::from_support(n, &anc.support, pauli_type);
            let weight = anc.support.len();

            generators.push(StabilizerGenerator {
                id: anc.id,
                kind: anc.kind,
                support: anc.support.clone(),
                pauli_string,
                is_bulk: anc.is_bulk,
                weight,
            });
        }

        let code = Self {
            lattice,
            generators,
        };

        // Mathematically verify that all stabilizer generators mutually commute:
        code.verify_commutation_algebra();

        code
    }

    /// Verifies that all stabilizer generators mutually commute: $[S_i, S_j] = 0 \quad \forall i, j$.
    pub fn verify_commutation_algebra(&self) {
        let m = self.generators.len();
        for i in 0..m {
            for j in (i + 1)..m {
                let g1 = &self.generators[i];
                let g2 = &self.generators[j];
                assert!(
                    g1.pauli_string.commutes_with(&g2.pauli_string),
                    "Stabilizer violation: S_{} (kind {:?}, coord {:?}) and S_{} (kind {:?}, coord {:?}) anticommute!",
                    i, g1.kind, self.lattice.ancillae[i].coord,
                    j, g2.kind, self.lattice.ancillae[j].coord,
                );
            }
        }
    }

    /// Total number of stabilizer generators: $d^2 - 1$.
    #[inline]
    pub fn num_generators(&self) -> usize {
        self.generators.len()
    }

    /// Measures the ideal syndrome for an $N$-qubit data Pauli error $E$.
    ///
    /// For each generator $S_i$:
    /// - If $E$ and $S_i$ commute: outcome is $+1$ (syndrome bit 0, no defect).
    /// - If $E$ and $S_i$ anticommute: outcome is $-1$ (syndrome bit 1, defect detected).
    pub fn measure_syndrome(&self, error: &PauliString) -> Vec<bool> {
        self.generators
            .iter()
            .map(|gen| gen.pauli_string.anticommutes_with(error))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stabilizer_commutation_distance_3() {
        let lattice = SurfaceLattice::new(3);
        let code = StabilizerCode::new(lattice);
        assert_eq!(code.num_generators(), 8);
        code.verify_commutation_algebra();
    }

    #[test]
    fn test_stabilizer_commutation_distance_5() {
        let lattice = SurfaceLattice::new(5);
        let code = StabilizerCode::new(lattice);
        assert_eq!(code.num_generators(), 24);
        code.verify_commutation_algebra();
    }

    #[test]
    fn test_single_qubit_error_syndromes() {
        let lattice = SurfaceLattice::new(3);
        let code = StabilizerCode::new(lattice);

        // Single X error on data qubit 0:
        // Qubit 0 is at (0, 0), which is in support of Z-stabilizer (-1, 1) and X-stabilizer (1, 1).
        // Since error is X, it anticommutes with Z-stabilizers that touch qubit 0!
        let err_x0 = PauliString::single_qubit(9, 0, Pauli::X);
        let syndrome = code.measure_syndrome(&err_x0);

        // Only Z-stabilizers containing qubit 0 should be triggered (syndrome bit true):
        for (i, &triggered) in syndrome.iter().enumerate() {
            let gen = &code.generators[i];
            if gen.support.contains(&0) && gen.kind == StabilizerType::Z {
                assert!(triggered);
            } else {
                assert!(!triggered);
            }
        }
    }
}
