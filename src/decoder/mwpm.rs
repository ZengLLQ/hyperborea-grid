//! Minimum-Weight Perfect Matching (MWPM) Decoder for rotated surface codes.
//!
//! Orchestrates spacetime defect graph construction, Edmonds' Blossom matching on general graphs,
//! boundary pairing resolution, and correction chain operator reconstruction.

use crate::code::lattice::SurfaceLattice;
use crate::code::logical::LogicalOperators;
use crate::code::pauli::PauliString;
use crate::code::stabilizer::StabilizerCode;
use crate::decoder::blossom::BlossomSolver;
use crate::decoder::boundary::BoundaryResolver;
use crate::decoder::correction::{CorrectionReconstructor, CorrectionResult};
use crate::syndrome::extraction::ExtractedDefects;
use crate::syndrome::graph::SpacetimeDefectGraph;
use crate::syndrome::weights::MetricWeights;

/// Minimum-Weight Perfect Matching decoder for fault-tolerant rotated surface codes.
#[derive(Debug, Clone)]
pub struct MwpmDecoder {
    /// Surface code lattice geometry.
    pub lattice: SurfaceLattice,
    /// Stabilizer code structure.
    pub code: StabilizerCode,
    /// Logical operators and fidelity evaluator.
    pub logical_ops: LogicalOperators,
    /// Metric weights for spacetime defect edges.
    pub weights: MetricWeights,
}

impl MwpmDecoder {
    /// Constructs an MWPM decoder for code distance $d$ and physical/measurement error probabilities.
    pub fn new(distance: usize, p: f64, p_meas: f64) -> Self {
        let lattice = SurfaceLattice::new(distance);
        let code = StabilizerCode::new(lattice.clone());
        let logical_ops = LogicalOperators::new(&lattice, &code);
        let weights = MetricWeights::from_probabilities(p, p_meas);

        Self {
            lattice,
            code,
            logical_ops,
            weights,
        }
    }

    /// Decodes an extracted spacetime defect structure given the ground-truth physical error.
    pub fn decode_defects(
        &self,
        extracted: &ExtractedDefects,
        true_data_error: &PauliString,
    ) -> CorrectionResult {
        let num_x_defects = extracted.x_defects.len();
        let num_z_defects = extracted.z_defects.len();

        let mut total_weight = 0.0;

        // 1. Solve MWPM for X-stabilizer defects (caused by Z errors):
        let x_pairs = if num_x_defects > 0 {
            let x_graph = SpacetimeDefectGraph::build(
                &self.lattice,
                extracted.x_defects.clone(),
                self.weights.clone(),
            );
            let x_solver = BlossomSolver::new(x_graph.total_vertices, x_graph.weight_matrix);
            let x_res = x_solver.solve();
            total_weight += x_res.total_weight as f64 / self.weights.scale_factor;
            BoundaryResolver::resolve_pairs(&extracted.x_defects, &x_res)
        } else {
            Vec::new()
        };

        // 2. Solve MWPM for Z-stabilizer defects (caused by X errors):
        let z_pairs = if num_z_defects > 0 {
            let z_graph = SpacetimeDefectGraph::build(
                &self.lattice,
                extracted.z_defects.clone(),
                self.weights.clone(),
            );
            let z_solver = BlossomSolver::new(z_graph.total_vertices, z_graph.weight_matrix);
            let z_res = z_solver.solve();
            total_weight += z_res.total_weight as f64 / self.weights.scale_factor;
            BoundaryResolver::resolve_pairs(&extracted.z_defects, &z_res)
        } else {
            Vec::new()
        };

        // 3. Reconstruct correction chain and evaluate logical status:
        CorrectionReconstructor::reconstruct_and_evaluate(
            &self.lattice,
            &self.code,
            &self.logical_ops,
            x_pairs,
            z_pairs,
            num_x_defects,
            num_z_defects,
            true_data_error,
            total_weight,
        )
    }

    /// Decodes a multi-round syndrome history.
    pub fn decode_syndromes(
        &self,
        measured_syndromes: &[Vec<bool>],
        true_data_error: &PauliString,
    ) -> CorrectionResult {
        let extracted = ExtractedDefects::from_syndromes(&self.lattice, measured_syndromes);
        self.decode_defects(&extracted, true_data_error)
    }

    /// Decodes single-round code-capacity syndrome.
    pub fn decode_single_round(
        &self,
        syndrome: &[bool],
        true_data_error: &PauliString,
    ) -> CorrectionResult {
        let round_0 = vec![false; syndrome.len()];
        let round_1 = syndrome.to_vec();
        let syndromes = vec![round_0, round_1];
        self.decode_syndromes(&syndromes, true_data_error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::pauli::Pauli;

    #[test]
    fn test_decoder_single_qubit_x_error_corrected() {
        let decoder = MwpmDecoder::new(3, 0.01, 0.01);
        let n = decoder.lattice.num_data_qubits;

        // Inject single X error on qubit 4 (the central data qubit of d=3)
        let true_err = PauliString::single_qubit(n, 4, Pauli::X);
        let syndrome = decoder.code.measure_syndrome(&true_err);

        let result = decoder.decode_single_round(&syndrome, &true_err);
        assert!(result.success, "Single qubit X error must be corrected");
        assert_eq!(result.logical_status, crate::code::logical::LogicalStatus::Success);
    }

    #[test]
    fn test_decoder_single_qubit_z_error_corrected() {
        let decoder = MwpmDecoder::new(3, 0.01, 0.01);
        let n = decoder.lattice.num_data_qubits;

        // Inject single Z error on qubit 4
        let true_err = PauliString::single_qubit(n, 4, Pauli::Z);
        let syndrome = decoder.code.measure_syndrome(&true_err);

        let result = decoder.decode_single_round(&syndrome, &true_err);
        assert!(result.success, "Single qubit Z error must be corrected");
    }

    #[test]
    fn test_decoder_d5_two_qubit_errors_corrected() {
        let decoder = MwpmDecoder::new(5, 0.01, 0.01);
        let n = decoder.lattice.num_data_qubits;

        // In d=5, code can correct any (d-1)/2 = 2 errors!
        let mut true_err = PauliString::single_qubit(n, 12, Pauli::X);
        let err2 = PauliString::single_qubit(n, 13, Pauli::X);
        true_err.multiply_in_place(&err2);

        let syndrome = decoder.code.measure_syndrome(&true_err);
        let result = decoder.decode_single_round(&syndrome, &true_err);
        assert!(result.success, "Two-qubit errors on d=5 must be corrected");
    }
}
