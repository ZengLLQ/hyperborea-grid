//! Integration tests for rotated planar surface code geometry, stabilizers, and logical operators.

use hyperborea_grid::code::lattice::{Boundary, StabilizerType, SurfaceLattice};
use hyperborea_grid::code::logical::{LogicalOperators, LogicalStatus};
use hyperborea_grid::code::pauli::PauliString;
use hyperborea_grid::code::stabilizer::StabilizerCode;

#[test]
fn test_rotated_surface_code_scaling_d3_d5_d7() {
    for &d in &[3, 5, 7] {
        let lattice = SurfaceLattice::new(d);
        let expected_data = d * d;
        let expected_stabilizers = d * d - 1;
        let expected_half = (d * d - 1) / 2;

        assert_eq!(lattice.distance, d);
        assert_eq!(lattice.num_data_qubits, expected_data);
        assert_eq!(lattice.ancillae.len(), expected_stabilizers);
        assert_eq!(lattice.x_ancillae.len(), expected_half);
        assert_eq!(lattice.z_ancillae.len(), expected_half);

        let code = StabilizerCode::new(lattice.clone());
        assert_eq!(code.num_generators(), expected_stabilizers);

        // Verify mutual commutation of all stabilizers
        code.verify_commutation_algebra();

        // Verify logical operators
        let ops = LogicalOperators::new(&lattice, &code);
        assert_eq!(ops.x_logical.weight(), d);
        assert_eq!(ops.z_logical.weight(), d);
        assert!(ops.x_logical.anticommutes_with(&ops.z_logical));
    }
}

#[test]
fn test_stabilizer_weights_and_boundaries() {
    let lattice = SurfaceLattice::new(3);

    let mut weight2_count = 0;
    let mut weight4_count = 0;

    for anc in &lattice.ancillae {
        match anc.support.len() {
            2 => {
                weight2_count += 1;
                assert!(!anc.is_bulk);
                assert!(anc.boundary.is_some());
                match anc.boundary.unwrap() {
                    Boundary::Top | Boundary::Bottom => {
                        assert_eq!(anc.kind, StabilizerType::Z);
                    }
                    Boundary::Left | Boundary::Right => {
                        assert_eq!(anc.kind, StabilizerType::X);
                    }
                }
            }
            4 => {
                weight4_count += 1;
                assert!(anc.is_bulk);
                assert!(anc.boundary.is_none());
            }
            other => panic!("Unexpected stabilizer weight: {}", other),
        }
    }

    assert_eq!(weight2_count, 4);
    assert_eq!(weight4_count, 4);
}

#[test]
fn test_logical_operator_commutation_with_stabilizers() {
    let lattice = SurfaceLattice::new(5);
    let code = StabilizerCode::new(lattice.clone());
    let ops = LogicalOperators::new(&lattice, &code);

    // Every stabilizer must commute with X_L and Z_L
    for gen in &code.generators {
        assert!(gen.pauli_string.commutes_with(&ops.x_logical));
        assert!(gen.pauli_string.commutes_with(&ops.z_logical));
    }
}

#[test]
fn test_logical_error_classification() {
    let lattice = SurfaceLattice::new(3);
    let code = StabilizerCode::new(lattice.clone());
    let ops = LogicalOperators::new(&lattice, &code);

    let id = PauliString::identity(9);
    assert_eq!(ops.check_residual_error(&id), LogicalStatus::Success);

    // Entire X_L chain applied as error:
    let x_err = ops.x_logical.clone();
    assert_eq!(ops.check_residual_error(&x_err), LogicalStatus::LogicalX);

    // Entire Z_L chain applied as error:
    let z_err = ops.z_logical.clone();
    assert_eq!(ops.check_residual_error(&z_err), LogicalStatus::LogicalZ);

    // Product of X_L and Z_L:
    let y_err = x_err.multiply(&z_err);
    assert_eq!(ops.check_residual_error(&y_err), LogicalStatus::LogicalY);
}
