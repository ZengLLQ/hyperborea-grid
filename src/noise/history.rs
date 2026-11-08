//! Spacetime error history tracking for fault-tolerant quantum error correction.
//!
//! Records physical Pauli faults on data qubits and measurement readout flips on ancillae
//! across consecutive syndrome extraction rounds.

use crate::code::pauli::{Pauli, PauliString};
use serde::{Deserialize, Serialize};

/// A physical quantum error event occurring at a specific space-time coordinate.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ErrorEvent {
    /// Syndrome extraction round $0 \le t < T$.
    pub round: usize,
    /// Target data qubit index (if data error).
    pub data_qubit: Option<usize>,
    /// Target ancilla id (if measurement readout error).
    pub ancilla_id: Option<usize>,
    /// Pauli operator applied (for data errors).
    pub pauli: Pauli,
    /// Whether this is a measurement readout flip.
    pub is_measurement_error: bool,
}

/// Complete spacetime history of injected errors throughout a multi-round simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ErrorHistory {
    /// Total data qubits in the code.
    pub num_data_qubits: usize,
    /// Total syndrome rounds simulated.
    pub num_rounds: usize,
    /// Log of all discrete error occurrences.
    pub events: Vec<ErrorEvent>,
    /// Net cumulative Pauli error string on all data qubits at the end of the simulation.
    pub cumulative_data_error: PauliString,
}

impl ErrorHistory {
    /// Creates an empty error history tracker for $N$ data qubits over $T$ rounds.
    pub fn new(num_data_qubits: usize, num_rounds: usize) -> Self {
        Self {
            num_data_qubits,
            num_rounds,
            events: Vec::new(),
            cumulative_data_error: PauliString::identity(num_data_qubits),
        }
    }

    /// Records a Pauli fault on a data qubit at round $t$.
    pub fn record_data_error(&mut self, round: usize, qubit: usize, pauli: Pauli) {
        if pauli.is_identity() {
            return;
        }

        self.events.push(ErrorEvent {
            round,
            data_qubit: Some(qubit),
            ancilla_id: None,
            pauli,
            is_measurement_error: false,
        });

        // Multiply into cumulative data error:
        let single = PauliString::single_qubit(self.num_data_qubits, qubit, pauli);
        self.cumulative_data_error.multiply_in_place(&single);
    }

    /// Records a measurement readout error on an ancilla at round $t$.
    pub fn record_measurement_error(&mut self, round: usize, ancilla_id: usize) {
        self.events.push(ErrorEvent {
            round,
            data_qubit: None,
            ancilla_id: Some(ancilla_id),
            pauli: Pauli::X, // Readout bit-flip
            is_measurement_error: true,
        });
    }

    /// Total number of physical data qubit Pauli errors injected.
    pub fn data_error_count(&self) -> usize {
        self.events
            .iter()
            .filter(|e| !e.is_measurement_error)
            .count()
    }

    /// Total number of measurement readout flips injected.
    pub fn measurement_error_count(&self) -> usize {
        self.events.iter().filter(|e| e.is_measurement_error).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error_history_tracking() {
        let mut history = ErrorHistory::new(9, 3);
        assert_eq!(history.data_error_count(), 0);
        assert_eq!(history.measurement_error_count(), 0);

        history.record_data_error(0, 2, Pauli::X);
        history.record_data_error(1, 2, Pauli::X); // X * X = I on qubit 2
        history.record_measurement_error(0, 1);

        assert_eq!(history.data_error_count(), 2);
        assert_eq!(history.measurement_error_count(), 1);
        assert_eq!(history.cumulative_data_error.get(2), Pauli::I);
    }
}
