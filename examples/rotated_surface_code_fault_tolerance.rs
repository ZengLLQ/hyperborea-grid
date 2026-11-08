//! Realistic demonstration of Distance-3 and Distance-5 Rotated Surface Code Fault-Tolerance.
//!
//! Demonstrates multi-round syndrome extraction ($T = d$), spacetime defect event identification,
//! Edmonds' Blossom MWPM decoding, correction chain reconstruction, and logical fidelity verification.

use hyperborea_grid::decoder::mwpm::MwpmDecoder;
use hyperborea_grid::noise::channel::NoiseModel;
use hyperborea_grid::noise::simulator::NoiseSimulator;
use hyperborea_grid::syndrome::extraction::ExtractedDefects;
use hyperborea_grid::visualizer::AsciiVisualizer;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use std::time::Instant;

fn main() {
    println!("================================================================================");
    println!("   HYPERBOREA-GRID: ROTATED PLANAR SURFACE CODE FAULT-TOLERANCE DEMONSTRATION   ");
    println!("================================================================================\n");

    run_distance_experiment(3, 0.006, 0.006, 101);
    println!("\n--------------------------------------------------------------------------------\n");
    run_distance_experiment(5, 0.004, 0.004, 202);
}

fn run_distance_experiment(distance: usize, p: f64, p_meas: f64, seed: u64) {
    println!(">>> RUNNING FAULT-TOLERANT EXTRACTION & MWPM DECODING (DISTANCE d = {}) <<<", distance);
    println!("Physical Parameters:");
    println!("  Code distance d:           {}", distance);
    println!("  Data qubits (d^2):         {}", distance * distance);
    println!("  Stabilizers (d^2 - 1):     {}", distance * distance - 1);
    println!("  Syndrome cycles T:         {}", distance);
    println!("  Physical error rate p:     {:.4} ({:.2}%)", p, p * 100.0);
    println!("  Measurement flip p_meas:   {:.4} ({:.2}%)", p_meas, p_meas * 100.0);
    println!("  PRNG Seed:                 {}", seed);

    let decoder = MwpmDecoder::new(distance, p, p_meas);
    let model = NoiseModel::depolarizing_with_meas(p, p_meas);
    let simulator = NoiseSimulator::new(model);
    let mut rng = ChaCha8Rng::seed_from_u64(seed);

    // Multi-round fault-tolerant simulation:
    let sim_start = Instant::now();
    let (history, syndromes) = simulator.simulate_multi_round(&decoder.code, distance, &mut rng);
    let sim_duration = sim_start.elapsed();

    let extracted = ExtractedDefects::from_syndromes(&decoder.lattice, &syndromes);

    // MWPM Blossom Decoding:
    let dec_start = Instant::now();
    let result = decoder.decode_defects(&extracted, &history.cumulative_data_error);
    let dec_duration = dec_start.elapsed();

    println!("\nSpacetime Syndrome Defect Statistics:");
    println!("  Injected Data Qubit Errors:     {}", history.data_error_count());
    println!("  Injected Readout Flips:         {}", history.measurement_error_count());
    println!("  Active X-Stabilizer Defects:    {}", extracted.x_defects.len());
    println!("  Active Z-Stabilizer Defects:    {}", extracted.z_defects.len());
    println!("  Total Spacetime Defect Nodes:   {}", extracted.all_defects.len());
    println!("  MWPM Solution Weight:           {:.4}", result.total_weight);

    println!("\nDecoded Pairings:");
    println!("  Z-Defect Pairs (bit-flips):     {}", result.z_pairs.len());
    println!("  X-Defect Pairs (phase-flips):   {}", result.x_pairs.len());

    println!("\nLogical Fidelity Verdict:");
    println!("  Logical State Preserved:        {}", result.success);
    println!("  Residual Logical Operator:      {:?}", result.logical_status);
    println!("  Correction Operator Weight:     {}", result.correction_chain.weight());
    println!("  Correction Chain String:        {}", result.correction_chain);

    println!("\nTiming Benchmarks:");
    println!("  Multi-round simulation time:    {:.3} ms", sim_duration.as_secs_f64() * 1000.0);
    println!("  MWPM Blossom decode time:       {:.3} ms", dec_duration.as_secs_f64() * 1000.0);

    println!("\n2D Rotated Surface Code Lattice Visualization:");
    println!(
        "{}",
        AsciiVisualizer::render_lattice(
            &decoder.lattice,
            &extracted.all_defects,
            Some(&result.correction_chain)
        )
    );
}
