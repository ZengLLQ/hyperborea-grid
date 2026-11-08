//! Phenomenological and circuit-level quantum Pauli noise models.
//!
//! Provides depolarizing, independent bit-flip / phase-flip, measurement readout error,
//! and circuit-level noise channels.

use crate::code::pauli::Pauli;
use rand::Rng;
use serde::{Deserialize, Serialize};

/// Type of quantum noise model applied during fault-tolerant syndrome cycles.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum NoiseModel {
    /// Standard single-parameter depolarizing noise:
    /// with probability $1 - p$: identity $I$,
    /// with probability $p/3$: $X$,
    /// with probability $p/3$: $Y$,
    /// with probability $p/3$: $Z$.
    Depolarizing {
        /// Physical error probability $p$.
        p: f64,
        /// Measurement readout flip probability $p_{\text{meas}}$.
        p_meas: f64,
    },

    /// Independent $X$ (bit-flip) and $Z$ (phase-flip) noise:
    /// each data qubit independently undergoes $X$ error with probability $p_X$,
    /// and $Z$ error with probability $p_Z$.
    IndependentXZ {
        /// Bit-flip probability $p_X$.
        px: f64,
        /// Phase-flip probability $p_Z$.
        pz: f64,
        /// Measurement readout flip probability $p_{\text{meas}}$.
        p_meas: f64,
    },

    /// Pure bit-flip noise ($X$ only).
    BitFlip {
        /// Bit-flip error probability $p$.
        p: f64,
        /// Measurement readout flip probability $p_{\text{meas}}$.
        p_meas: f64,
    },

    /// Pure phase-flip noise ($Z$ only).
    PhaseFlip {
        /// Phase-flip error probability $p$.
        p: f64,
        /// Measurement readout flip probability $p_{\text{meas}}$.
        p_meas: f64,
    },

    /// Circuit-level noise parameterization for gate and readout operations.
    CircuitLevel {
        /// Single-qubit gate error rate.
        p_single: f64,
        /// Two-qubit entangling gate error rate (15 Pauli combinations).
        p_two: f64,
        /// Ancilla measurement readout error rate.
        p_meas: f64,
        /// Ancilla state reset error rate.
        p_reset: f64,
    },
}

impl NoiseModel {
    /// Creates a standard phenomenological depolarizing noise channel with equal data and measurement error rates $p$.
    pub fn depolarizing(p: f64) -> Self {
        Self::Depolarizing { p, p_meas: p }
    }

    /// Creates a depolarizing noise channel with distinct data and measurement error rates.
    pub fn depolarizing_with_meas(p: f64, p_meas: f64) -> Self {
        Self::Depolarizing { p, p_meas }
    }

    /// Creates an independent $X$ and $Z$ noise channel.
    pub fn independent_xz(px: f64, pz: f64, p_meas: f64) -> Self {
        Self::IndependentXZ { px, pz, p_meas }
    }

    /// Returns the data qubit physical error rate $p$.
    pub fn physical_error_rate(&self) -> f64 {
        match *self {
            Self::Depolarizing { p, .. } => p,
            Self::IndependentXZ { px, pz, .. } => (px + pz) / 2.0,
            Self::BitFlip { p, .. } => p,
            Self::PhaseFlip { p, .. } => p,
            Self::CircuitLevel { p_two, .. } => p_two,
        }
    }

    /// Returns the measurement readout error rate $p_{\text{meas}}$.
    pub fn measurement_error_rate(&self) -> f64 {
        match *self {
            Self::Depolarizing { p_meas, .. } => p_meas,
            Self::IndependentXZ { p_meas, .. } => p_meas,
            Self::BitFlip { p_meas, .. } => p_meas,
            Self::PhaseFlip { p_meas, .. } => p_meas,
            Self::CircuitLevel { p_meas, .. } => p_meas,
        }
    }

    /// Samples a single-qubit Pauli error for a data qubit in one round according to this model.
    pub fn sample_data_qubit_error<R: Rng + ?Sized>(&self, rng: &mut R) -> Pauli {
        match *self {
            Self::Depolarizing { p, .. } => {
                let r: f64 = rng.gen();
                if r < p / 3.0 {
                    Pauli::X
                } else if r < 2.0 * p / 3.0 {
                    Pauli::Y
                } else if r < p {
                    Pauli::Z
                } else {
                    Pauli::I
                }
            }
            Self::IndependentXZ { px, pz, .. } => {
                let rx: f64 = rng.gen();
                let rz: f64 = rng.gen();
                let has_x = rx < px;
                let has_z = rz < pz;
                match (has_x, has_z) {
                    (false, false) => Pauli::I,
                    (true, false) => Pauli::X,
                    (true, true) => Pauli::Y,
                    (false, true) => Pauli::Z,
                }
            }
            Self::BitFlip { p, .. } => {
                if rng.gen::<f64>() < p {
                    Pauli::X
                } else {
                    Pauli::I
                }
            }
            Self::PhaseFlip { p, .. } => {
                if rng.gen::<f64>() < p {
                    Pauli::Z
                } else {
                    Pauli::I
                }
            }
            Self::CircuitLevel { p_single, .. } => {
                let r: f64 = rng.gen();
                if r < p_single / 3.0 {
                    Pauli::X
                } else if r < 2.0 * p_single / 3.0 {
                    Pauli::Y
                } else if r < p_single {
                    Pauli::Z
                } else {
                    Pauli::I
                }
            }
        }
    }

    /// Samples whether a measurement readout outcome is flipped at round $t$.
    pub fn sample_measurement_flip<R: Rng + ?Sized>(&self, rng: &mut R) -> bool {
        rng.gen::<f64>() < self.measurement_error_rate()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_chacha::ChaCha8Rng;

    #[test]
    fn test_depolarizing_sampling_frequencies() {
        let mut rng = ChaCha8Rng::seed_from_u64(42);
        let p = 0.3;
        let model = NoiseModel::depolarizing(p);

        let n = 100_000;
        let mut count_x = 0;
        let mut count_y = 0;
        let mut count_z = 0;
        let mut count_i = 0;

        for _ in 0..n {
            match model.sample_data_qubit_error(&mut rng) {
                Pauli::I => count_i += 1,
                Pauli::X => count_x += 1,
                Pauli::Y => count_y += 1,
                Pauli::Z => count_z += 1,
            }
        }

        // Each of X, Y, Z should be close to p/3 = 10%
        let fx = count_x as f64 / n as f64;
        let fy = count_y as f64 / n as f64;
        let fz = count_z as f64 / n as f64;
        let fi = count_i as f64 / n as f64;

        assert!((fx - 0.1).abs() < 0.01);
        assert!((fy - 0.1).abs() < 0.01);
        assert!((fz - 0.1).abs() < 0.01);
        assert!((fi - 0.7).abs() < 0.01);
    }
}
