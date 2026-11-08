//! Integration tests for fault-tolerant error correction and threshold scaling.

use hyperborea_grid::code::pauli::{Pauli, PauliString};
use hyperborea_grid::decoder::mwpm::MwpmDecoder;
use hyperborea_grid::noise::channel::NoiseModel;
use hyperborea_grid::noise::simulator::NoiseSimulator;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

#[test]
fn test_exhaustive_single_qubit_correction_distance_3() {
    let decoder = MwpmDecoder::new(3, 0.005, 0.005);
    let n = decoder.lattice.num_data_qubits;

    // For code distance 3, every single-qubit Pauli error (X, Y, Z) on ANY of the 9 data qubits
    // MUST be 100% corrected without logical fault.
    for q in 0..n {
        for &p in &[Pauli::X, Pauli::Y, Pauli::Z] {
            let err = PauliString::single_qubit(n, q, p);
            let syn = decoder.code.measure_syndrome(&err);
            let res = decoder.decode_single_round(&syn, &err);

            assert!(
                res.success,
                "Failed to correct single-qubit Pauli {:?} on qubit {}",
                p, q
            );
        }
    }
}

#[test]
fn test_distance_5_two_qubit_arbitrary_correction() {
    let decoder = MwpmDecoder::new(5, 0.005, 0.005);
    let n = decoder.lattice.num_data_qubits;

    // Distance 5 can correct any t = (5-1)/2 = 2 errors.
    // Test pairs of X errors across various qubit combinations:
    let test_pairs = [
        (0, 1),
        (0, 12),
        (6, 18),
        (12, 13),
        (10, 14),
        (23, 24),
    ];

    for &(q1, q2) in &test_pairs {
        let mut err = PauliString::single_qubit(n, q1, Pauli::X);
        let err2 = PauliString::single_qubit(n, q2, Pauli::X);
        err.multiply_in_place(&err2);

        let syn = decoder.code.measure_syndrome(&err);
        let res = decoder.decode_single_round(&syn, &err);
        assert!(
            res.success,
            "Failed to correct two-qubit error on qubits {} and {}",
            q1, q2
        );
    }
}

#[test]
fn test_fault_tolerant_multi_round_decoding_with_measurement_noise() {
    let decoder = MwpmDecoder::new(3, 0.005, 0.005);
    let model = NoiseModel::depolarizing_with_meas(0.005, 0.005);
    let sim = NoiseSimulator::new(model);
    let mut rng = ChaCha8Rng::seed_from_u64(12345);

    let trials = 100;
    let mut successes = 0;

    for _ in 0..trials {
        let (history, syndromes) = sim.simulate_multi_round(&decoder.code, 3, &mut rng);
        let res = decoder.decode_syndromes(&syndromes, &history.cumulative_data_error);
        if res.success {
            successes += 1;
        }
    }

    // At low error rate p = 0.005 on d=3, success rate should be >= 95%
    let success_rate = successes as f64 / trials as f64;
    assert!(
        success_rate >= 0.95,
        "Success rate {:.2} is below 95% threshold for p=0.005",
        success_rate
    );
}
