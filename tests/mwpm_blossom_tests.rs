//! Integration tests for the Edmonds' Blossom Minimum-Weight Perfect Matching decoder.

use hyperborea_grid::decoder::blossom::BlossomSolver;
use hyperborea_grid::decoder::boundary::MatchedPair;
use hyperborea_grid::decoder::mwpm::MwpmDecoder;
use hyperborea_grid::code::pauli::{Pauli, PauliString};

#[test]
fn test_blossom_solver_complete_graph_k4() {
    let weights = vec![
        vec![0, 5, 2, 7],
        vec![5, 0, 8, 3],
        vec![2, 8, 0, 4],
        vec![7, 3, 4, 0],
    ];

    let solver = BlossomSolver::new(4, weights);
    let res = solver.solve();
    // Optimal matching: (0,2) with weight 2, and (1,3) with weight 3 => total = 5.
    assert_eq!(res.total_weight, 5);
    assert_eq!(res.pairs, vec![(0, 2), (1, 3)]);
}

#[test]
fn test_blossom_solver_complete_graph_k6_non_bipartite() {
    let weights = vec![
        vec![0,  3, 10, 15, 20, 25],
        vec![3,  0,  4, 18, 22, 28],
        vec![10, 4,  0,  5, 16, 21],
        vec![15, 18, 5,  0,  6, 17],
        vec![20, 22, 16, 6,  0,  7],
        vec![25, 28, 21, 17, 7,  0],
    ];

    let solver = BlossomSolver::new(6, weights);
    let bb = solver.solve_branch_and_bound();
    let pd = solver.solve_blossom_primal_dual();

    assert_eq!(bb.total_weight, pd.total_weight);
    // (0,1)=3, (2,3)=5, (4,5)=7 => total = 15
    assert_eq!(bb.total_weight, 15);
    assert_eq!(bb.pairs, vec![(0, 1), (2, 3), (4, 5)]);
}

#[test]
fn test_boundary_matching_resolution() {
    let decoder = MwpmDecoder::new(3, 0.01, 0.01);
    let n = decoder.lattice.num_data_qubits;

    // A single X error on boundary data qubit 0 (top-left corner).
    // Qubit 0 touches top Z-boundary ancilla (-1, 1).
    // Because it touches only one Z-stabilizer, only ONE defect is created in the bulk!
    let err = PauliString::single_qubit(n, 0, Pauli::X);
    let syn = decoder.code.measure_syndrome(&err);

    let res = decoder.decode_single_round(&syn, &err);
    assert!(res.success, "Boundary defect must match to boundary and succeed");

    // Must have matched to boundary:
    assert_eq!(res.z_pairs.len(), 1);
    match &res.z_pairs[0] {
        MatchedPair::BoundaryPair { defect, .. } => {
            assert_eq!(defect.ancilla_id, 0); // The top Z ancilla
        }
        _ => panic!("Expected boundary pair matching"),
    }
}
