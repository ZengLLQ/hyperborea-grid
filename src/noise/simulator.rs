//! Quantum noise simulator across multi-round stabilizer cycles.
//!
//! Applies physical Pauli error processes to data qubits and injects measurement readout flips
//! into stabilizer ancillae, producing recorded syndrome measurements and spacetime error histories.

use super::channel::NoiseModel;
use super::history::ErrorHistory;
use crate::code::stabilizer::StabilizerCode;
use rand::Rng;

/// Quantum error simulation engine.
#[derive(Debug, Clone)]
pub struct NoiseSimulator {
    /// Noise model governing physical and measurement error rates.
    pub model: NoiseModel,
}

impl NoiseSimulator {
    /// Creates a new noise simulator with the specified noise channel.
    pub fn new(model: NoiseModel) -> Self {
        Self { model }
    }

    /// Simulates $T$ rounds of noisy syndrome extraction on the rotated surface code.
    ///
    /// At each round $t \in \{0, \dots, T - 1\}$:
    /// 1. Data qubits accumulate stochastic Pauli errors according to the noise model.
    /// 2. Stabilizer generators are measured against the current cumulative data error.
    /// 3. Measurement outcomes are subject to readout bit-flip errors.
    ///
    /// At round $T$ (terminal round), an ideal noiseless readout is performed to ensure
    /// spacetime defect graph closure.
    ///
    /// Returns:
    /// - `ErrorHistory`: full log of physical errors and the final cumulative data Pauli error.
    /// - `Vec<Vec<bool>>`: noisy syndrome bitvectors of length $d^2 - 1$ for rounds $0, \dots, T$.
    pub fn simulate_multi_round<R: Rng + ?Sized>(
        &self,
        code: &StabilizerCode,
        num_rounds: usize,
        rng: &mut R,
    ) -> (ErrorHistory, Vec<Vec<bool>>) {
        let n = code.lattice.num_data_qubits;
        let m = code.num_generators();
        let mut history = ErrorHistory::new(n, num_rounds);
        let mut measured_syndromes = Vec::with_capacity(num_rounds + 1);

        for t in 0..num_rounds {
            // 1. Inject data qubit errors:
            for q in 0..n {
                let pauli = self.model.sample_data_qubit_error(rng);
                history.record_data_error(t, q, pauli);
            }

            // 2. Measure ideal syndrome against cumulative data error:
            let ideal_syndrome = code.measure_syndrome(&history.cumulative_data_error);

            // 3. Inject measurement readout flip:
            let mut noisy_syndrome = ideal_syndrome;
            for a in 0..m {
                if self.model.sample_measurement_flip(rng) {
                    noisy_syndrome[a] = !noisy_syndrome[a];
                    history.record_measurement_error(t, a);
                }
            }

            measured_syndromes.push(noisy_syndrome);
        }

        // Terminal round T: noiseless final syndrome (representing data qubit basis readout):
        let final_noiseless_syndrome = code.measure_syndrome(&history.cumulative_data_error);
        measured_syndromes.push(final_noiseless_syndrome);

        (history, measured_syndromes)
    }

    /// Single-round code-capacity noise simulation (noiseless measurements, data errors only).
    pub fn simulate_code_capacity<R: Rng + ?Sized>(
        &self,
        code: &StabilizerCode,
        rng: &mut R,
    ) -> (ErrorHistory, Vec<bool>) {
        let n = code.lattice.num_data_qubits;
        let mut history = ErrorHistory::new(n, 1);

        for q in 0..n {
            let pauli = self.model.sample_data_qubit_error(rng);
            history.record_data_error(0, q, pauli);
        }

        let syndrome = code.measure_syndrome(&history.cumulative_data_error);
        (history, syndrome)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::code::lattice::SurfaceLattice;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_noiseless_simulation_gives_zero_defects() {
        let lattice = SurfaceLattice::new(3);
        let code = StabilizerCode::new(lattice);
        let sim = NoiseSimulator::new(NoiseModel::depolarizing(0.0));
        let mut rng = ChaCha8Rng::seed_from_u64(123);

        let (history, syndromes) = sim.simulate_multi_round(&code, 3, &mut rng);
        assert_eq!(history.data_error_count(), 0);
        assert_eq!(history.measurement_error_count(), 0);

        for round_syn in syndromes {
            assert!(round_syn.iter().all(|&b| !b));
        }
    }
}
