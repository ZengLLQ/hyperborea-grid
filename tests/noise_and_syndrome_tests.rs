//! Integration tests for Pauli noise channels, error tracking, and syndrome defect extraction.

use hyperborea_grid::code::lattice::{StabilizerType, SurfaceLattice};
use hyperborea_grid::code::pauli::{Pauli, PauliString};
use hyperborea_grid::code::stabilizer::StabilizerCode;
use hyperborea_grid::noise::channel::NoiseModel;
use hyperborea_grid::noise::simulator::NoiseSimulator;
use hyperborea_grid::syndrome::extraction::ExtractedDefects;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

#[test]
fn test_pauli_error_algebra_and_syndrome_excitation() {
    let lattice = SurfaceLattice::new(3);
    let code = StabilizerCode::new(lattice);

    // Test X error on each data qubit:
    // Every X error must anticommute with at least one Z stabilizer and commute with all X stabilizers.
    for q in 0..9 {
        let err_x = PauliString::single_qubit(9, q, Pauli::X);
        let syn = code.measure_syndrome(&err_x);

        let mut triggered_z = 0;
        let mut triggered_x = 0;

        for (i, &trig) in syn.iter().enumerate() {
            if trig {
                match code.generators[i].kind {
                    StabilizerType::X => triggered_x += 1,
                    StabilizerType::Z => triggered_z += 1,
                }
            }
        }

        assert_eq!(triggered_x, 0, "X error must not trigger X stabilizers");
        assert!(triggered_z >= 1 && triggered_z <= 2, "X error in d=3 must trigger 1 or 2 Z stabilizers");
    }

    // Test Z error on each data qubit:
    // Every Z error must anticommute with at least one X stabilizer and commute with all Z stabilizers.
    for q in 0..9 {
        let err_z = PauliString::single_qubit(9, q, Pauli::Z);
        let syn = code.measure_syndrome(&err_z);

        let mut triggered_z = 0;
        let mut triggered_x = 0;

        for (i, &trig) in syn.iter().enumerate() {
            if trig {
                match code.generators[i].kind {
                    StabilizerType::X => triggered_x += 1,
                    StabilizerType::Z => triggered_z += 1,
                }
            }
        }

        assert_eq!(triggered_z, 0, "Z error must not trigger Z stabilizers");
        assert!(triggered_x >= 1 && triggered_x <= 2, "Z error in d=3 must trigger 1 or 2 X stabilizers");
    }
}

#[test]
fn test_measurement_error_creates_timelike_defect_pairs() {
    let lattice = SurfaceLattice::new(3);
    let _code = StabilizerCode::new(lattice.clone());

    // Construct 3 rounds with zero data error, but a single measurement readout flip on ancilla 2 at round 1:
    let syn_round_0 = vec![false; 8];
    let mut syn_round_1 = vec![false; 8];
    syn_round_1[2] = true; // Readout flip
    let syn_round_2 = vec![false; 8]; // Ideal terminal readout

    let syndromes = vec![syn_round_0, syn_round_1, syn_round_2];
    let extracted = ExtractedDefects::from_syndromes(&lattice, &syndromes);

    // The single measurement error at round 1 causes two defect events in time:
    // one at round 1 and one at round 2, both on ancilla 2!
    assert_eq!(extracted.all_defects.len(), 2);
    assert_eq!(extracted.all_defects[0].ancilla_id, 2);
    assert_eq!(extracted.all_defects[0].round, 1);
    assert_eq!(extracted.all_defects[1].ancilla_id, 2);
    assert_eq!(extracted.all_defects[1].round, 2);
}

#[test]
fn test_multi_round_noise_simulation_history() {
    let lattice = SurfaceLattice::new(3);
    let code = StabilizerCode::new(lattice);
    let model = NoiseModel::depolarizing_with_meas(0.02, 0.02);
    let sim = NoiseSimulator::new(model);
    let mut rng = ChaCha8Rng::seed_from_u64(777);

    let (history, syndromes) = sim.simulate_multi_round(&code, 5, &mut rng);
    assert_eq!(history.num_rounds, 5);
    assert_eq!(syndromes.len(), 6); // 5 noisy rounds + 1 final noiseless readout
    assert_eq!(history.cumulative_data_error.num_qubits(), 9);
}
