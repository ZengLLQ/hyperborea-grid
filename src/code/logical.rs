//! Logical operators ($X_L$ and $Z_L$) and logical error classification.
//!
//! Provides the weight-$d$ logical operator representatives for the rotated surface code,
//! verifies that they commute with all stabilizer generators and anticommute with each other,
//! and evaluates whether a residual error causes a logical fault.

use super::lattice::SurfaceLattice;
use super::pauli::{Pauli, PauliString};
use super::stabilizer::StabilizerCode;
use serde::{Deserialize, Serialize};

/// Type of logical Pauli error occurring on the encoded qubit.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum LogicalStatus {
    /// Preserved logical state ($I_L$). Successful error correction.
    Success,
    /// Logical bit-flip fault ($X_L$).
    LogicalX,
    /// Logical phase-flip fault ($Z_L$).
    LogicalZ,
    /// Logical bit-and-phase fault ($Y_L$).
    LogicalY,
}

impl LogicalStatus {
    /// Returns true if no logical error occurred ($I_L$).
    #[inline]
    pub fn is_success(&self) -> bool {
        matches!(self, LogicalStatus::Success)
    }
}

/// Logical operators ($X_L, Z_L$) for the encoded logical qubit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LogicalOperators {
    /// Code distance $d$.
    pub distance: usize,
    /// Total data qubits $N = d^2$.
    pub num_data_qubits: usize,
    /// Logical $X$ operator ($X_L$), weight $d$.
    pub x_logical: PauliString,
    /// Logical $Z$ operator ($Z_L$), weight $d$.
    pub z_logical: PauliString,
    /// Qubits participating in $X_L$.
    pub x_support: Vec<usize>,
    /// Qubits participating in $Z_L$.
    pub z_support: Vec<usize>,
}

impl LogicalOperators {
    /// Constructs and validates logical operators for a given surface lattice and stabilizer code.
    pub fn new(lattice: &SurfaceLattice, code: &StabilizerCode) -> Self {
        let d = lattice.distance;
        let n = lattice.num_data_qubits;

        // X_L: row 0 of data qubits, connecting left to right rough boundaries.
        // Support: data qubits (0, c) for c in 0..d.
        let mut x_support = Vec::with_capacity(d);
        for c in 0..d {
            x_support.push(lattice.data_qubit_index(0, c));
        }
        let x_logical = PauliString::from_support(n, &x_support, Pauli::X);

        // Z_L: col 0 of data qubits, connecting top to bottom smooth boundaries.
        // Support: data qubits (r, 0) for r in 0..d.
        let mut z_support = Vec::with_capacity(d);
        for r in 0..d {
            z_support.push(lattice.data_qubit_index(r, 0));
        }
        let z_logical = PauliString::from_support(n, &z_support, Pauli::Z);

        let ops = Self {
            distance: d,
            num_data_qubits: n,
            x_logical,
            z_logical,
            x_support,
            z_support,
        };

        // Mathematically verify that X_L and Z_L anticommute:
        assert!(
            ops.x_logical.anticommutes_with(&ops.z_logical),
            "Logical operator error: X_L and Z_L must anticommute! {{X_L, Z_L}} != 0"
        );

        // Verify that X_L and Z_L commute with every stabilizer generator:
        for (i, gen) in code.generators.iter().enumerate() {
            assert!(
                ops.x_logical.commutes_with(&gen.pauli_string),
                "Logical operator error: X_L does not commute with stabilizer S_{}",
                i
            );
            assert!(
                ops.z_logical.commutes_with(&gen.pauli_string),
                "Logical operator error: Z_L does not commute with stabilizer S_{}",
                i
            );
        }

        ops
    }

    /// Evaluates the logical status of a residual Pauli error operator $E_{\text{res}}$.
    ///
    /// - If $E_{\text{res}}$ anticommutes with $Z_L$, an odd number of $X$ errors traversed the lattice $\implies$ Logical $X$.
    /// - If $E_{\text{res}}$ anticommutes with $X_L$, an odd number of $Z$ errors traversed the lattice $\implies$ Logical $Z$.
    /// - If it anticommutes with both $\implies$ Logical $Y$.
    /// - If it commutes with both $\implies$ Success ($I_L$).
    pub fn check_residual_error(&self, residual_error: &PauliString) -> LogicalStatus {
        let flip_x = residual_error.anticommutes_with(&self.z_logical);
        let flip_z = residual_error.anticommutes_with(&self.x_logical);

        match (flip_x, flip_z) {
            (false, false) => LogicalStatus::Success,
            (true, false) => LogicalStatus::LogicalX,
            (false, true) => LogicalStatus::LogicalZ,
            (true, true) => LogicalStatus::LogicalY,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_logical_operators_distance_3() {
        let lattice = SurfaceLattice::new(3);
        let code = StabilizerCode::new(lattice.clone());
        let ops = LogicalOperators::new(&lattice, &code);

        assert_eq!(ops.x_logical.weight(), 3);
        assert_eq!(ops.z_logical.weight(), 3);
        assert!(ops.x_logical.anticommutes_with(&ops.z_logical));

        // Identity error is success:
        let id_err = PauliString::identity(9);
        assert_eq!(ops.check_residual_error(&id_err), LogicalStatus::Success);

        // Single stabilizer as error should also be logical success (stabilizers act trivially on codespace):
        for gen in &code.generators {
            assert_eq!(
                ops.check_residual_error(&gen.pauli_string),
                LogicalStatus::Success
            );
        }

        // Applying full X_L as error should be detected as LogicalX:
        assert_eq!(
            ops.check_residual_error(&ops.x_logical),
            LogicalStatus::LogicalX
        );

        // Applying full Z_L as error should be detected as LogicalZ:
        assert_eq!(
            ops.check_residual_error(&ops.z_logical),
            LogicalStatus::LogicalZ
        );
    }

    #[test]
    fn test_logical_operators_distance_5() {
        let lattice = SurfaceLattice::new(5);
        let code = StabilizerCode::new(lattice.clone());
        let ops = LogicalOperators::new(&lattice, &code);

        assert_eq!(ops.x_logical.weight(), 5);
        assert_eq!(ops.z_logical.weight(), 5);
        assert!(ops.x_logical.anticommutes_with(&ops.z_logical));
    }
}
