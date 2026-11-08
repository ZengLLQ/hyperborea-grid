//! Correction chain operator reconstruction and logical fidelity evaluation.
//!
//! Maps matched defect pairs and boundary terminations into Pauli correction strings,
//! applies corrections to data qubits, verifies defect clearance, and computes logical error status.

use crate::code::lattice::SurfaceLattice;
use crate::code::logical::{LogicalOperators, LogicalStatus};
use crate::code::pauli::{Pauli, PauliString};
use crate::code::stabilizer::StabilizerCode;
use crate::decoder::boundary::MatchedPair;
use crate::decoder::dijkstra::DijkstraSpacetime;
use serde::{Deserialize, Serialize};

/// Detailed outcome of the quantum error correction decoding cycle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CorrectionResult {
    /// Whether the logical quantum information was preserved without logical fault.
    pub success: bool,
    /// Logical Pauli status of the residual error ($I_L, X_L, Z_L, Y_L$).
    pub logical_status: LogicalStatus,
    /// The reconstructed Pauli correction operator applied to data qubits.
    pub correction_chain: PauliString,
    /// The net residual error operator $E_{\text{res}} = E_{\text{true}} \cdot C$.
    pub residual_error: PauliString,
    /// Number of physical defects corrected for X-type and Z-type stabilizers.
    pub num_x_defects: usize,
    pub num_z_defects: usize,
    /// Physical pairings resolved for X-stabilizer defects.
    pub x_pairs: Vec<MatchedPair>,
    /// Physical pairings resolved for Z-stabilizer defects.
    pub z_pairs: Vec<MatchedPair>,
    /// Total MWPM weight score.
    pub total_weight: f64,
}

/// Engine for reconstructing Pauli correction chains from matched defect pairs.
pub struct CorrectionReconstructor;

impl CorrectionReconstructor {
    /// Reconstructs full Pauli correction operators for both X and Z syndromes,
    /// applies them to data qubits, and checks logical fidelity.
    pub fn reconstruct_and_evaluate(
        lattice: &SurfaceLattice,
        code: &StabilizerCode,
        logical_ops: &LogicalOperators,
        x_pairs: Vec<MatchedPair>,
        z_pairs: Vec<MatchedPair>,
        num_x_defects: usize,
        num_z_defects: usize,
        true_data_error: &PauliString,
        total_weight: f64,
    ) -> CorrectionResult {
        let n = lattice.num_data_qubits;
        let mut correction = PauliString::identity(n);

        // 1. Z-stabilizer defects are caused by X errors on data qubits.
        // Therefore, matched pairs of Z-stabilizer defects yield Pauli X corrections!
        for pair in &z_pairs {
            match pair {
                MatchedPair::DefectPair { defect1, defect2 } => {
                    if defect1.ancilla_id != defect2.ancilla_id {
                        let path = DijkstraSpacetime::shortest_path(
                            lattice,
                            defect1.ancilla_id,
                            defect2.ancilla_id,
                        );
                        for &q in &path {
                            let curr = correction.get(q);
                            let new_p = match curr {
                                Pauli::I => Pauli::X,
                                Pauli::X => Pauli::I,
                                Pauli::Y => Pauli::Z,
                                Pauli::Z => Pauli::Y,
                            };
                            correction.set(q, new_p);
                        }
                    }
                }
                MatchedPair::BoundaryPair { defect, .. } => {
                    let path = DijkstraSpacetime::shortest_boundary_path(lattice, defect.ancilla_id);
                    for &q in &path {
                        let curr = correction.get(q);
                        let new_p = match curr {
                            Pauli::I => Pauli::X,
                            Pauli::X => Pauli::I,
                            Pauli::Y => Pauli::Z,
                            Pauli::Z => Pauli::Y,
                        };
                        correction.set(q, new_p);
                    }
                }
            }
        }

        // 2. X-stabilizer defects are caused by Z errors on data qubits.
        // Therefore, matched pairs of X-stabilizer defects yield Pauli Z corrections!
        for pair in &x_pairs {
            match pair {
                MatchedPair::DefectPair { defect1, defect2 } => {
                    if defect1.ancilla_id != defect2.ancilla_id {
                        let path = DijkstraSpacetime::shortest_path(
                            lattice,
                            defect1.ancilla_id,
                            defect2.ancilla_id,
                        );
                        for &q in &path {
                            let curr = correction.get(q);
                            let new_p = match curr {
                                Pauli::I => Pauli::Z,
                                Pauli::Z => Pauli::I,
                                Pauli::Y => Pauli::X,
                                Pauli::X => Pauli::Y,
                            };
                            correction.set(q, new_p);
                        }
                    }
                }
                MatchedPair::BoundaryPair { defect, .. } => {
                    let path = DijkstraSpacetime::shortest_boundary_path(lattice, defect.ancilla_id);
                    for &q in &path {
                        let curr = correction.get(q);
                        let new_p = match curr {
                            Pauli::I => Pauli::Z,
                            Pauli::Z => Pauli::I,
                            Pauli::Y => Pauli::X,
                            Pauli::X => Pauli::Y,
                        };
                        correction.set(q, new_p);
                    }
                }
            }
        }

        // 3. Compute net residual error operator: E_res = E_true * C
        let residual_error = true_data_error.multiply(&correction);

        // 4. Verify that residual error produces no active defects in the stabilizer group:
        let residual_syndrome = code.measure_syndrome(&residual_error);
        assert!(
            residual_syndrome.iter().all(|&defect| !defect),
            "Correction chain failed to clear syndrome! Residual syndrome: {:?}",
            residual_syndrome
        );

        // 5. Evaluate logical status against logical operators X_L and Z_L:
        let logical_status = logical_ops.check_residual_error(&residual_error);
        let success = logical_status.is_success();

        CorrectionResult {
            success,
            logical_status,
            correction_chain: correction,
            residual_error,
            num_x_defects,
            num_z_defects,
            x_pairs,
            z_pairs,
            total_weight,
        }
    }
}
