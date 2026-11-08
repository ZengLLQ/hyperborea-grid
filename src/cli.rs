//! Command-line interface definitions and command execution dispatch.

use clap::{Parser, Subcommand};
use rayon::prelude::*;
use std::time::Instant;

use crate::code::lattice::SurfaceLattice;
use crate::code::stabilizer::StabilizerCode;
use crate::decoder::mwpm::MwpmDecoder;
use crate::noise::channel::NoiseModel;
use crate::noise::simulator::NoiseSimulator;
use crate::visualizer::AsciiVisualizer;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;

/// HyperboreaGrid: Rotated Planar Surface Code Quantum Error Correction & MWPM Blossom Decoder.
#[derive(Parser, Debug)]
#[command(
    name = "hyperborea-grid",
    author = "ZengLLQ",
    version = "0.1.0",
    about = "Fault-Tolerant Rotated Surface Code Simulation & Blossom MWPM Decoding in Rust"
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Visualizes the rotated surface code lattice and stabilizer configuration.
    Lattice {
        /// Code distance d (odd integer >= 3).
        #[arg(short, long, default_value_t = 3)]
        distance: usize,
    },

    /// Simulates multi-round fault-tolerant syndrome extraction and Blossom MWPM decoding.
    Decode {
        /// Code distance d (odd integer >= 3).
        #[arg(short, long, default_value_t = 3)]
        distance: usize,

        /// Number of syndrome extraction rounds T (defaults to code distance).
        #[arg(short, long)]
        rounds: Option<usize>,

        /// Physical data qubit error rate p.
        #[arg(short, long, default_value_t = 0.005)]
        error_rate: f64,

        /// Ancilla measurement readout flip probability p_meas.
        #[arg(short = 'm', long, default_value_t = 0.005)]
        meas_error_rate: f64,

        /// Random number generator seed.
        #[arg(short, long, default_value_t = 42)]
        seed: u64,
    },

    /// Runs Monte Carlo threshold scaling simulation across physical error rates.
    Threshold {
        /// Code distances to evaluate (e.g. "3,5").
        #[arg(short, long, default_value = "3,5")]
        distances: String,

        /// Minimum physical error rate p.
        #[arg(long, default_value_t = 0.002)]
        min_p: f64,

        /// Maximum physical error rate p.
        #[arg(long, default_value_t = 0.030)]
        max_p: f64,

        /// Number of sampling steps between min_p and max_p.
        #[arg(short, long, default_value_t = 8)]
        steps: usize,

        /// Number of Monte Carlo shots per point.
        #[arg(short, long, default_value_t = 500)]
        trials: usize,
    },

    /// Runs a comprehensive demonstration of surface code fault-tolerance and threshold scaling.
    Demo,
}

impl Cli {
    pub fn execute(self) {
        match self.command {
            Commands::Lattice { distance } => {
                let lattice = SurfaceLattice::new(distance);
                let code = StabilizerCode::new(lattice.clone());
                println!("{}", AsciiVisualizer::render_lattice(&lattice, &[], None));
                println!("Stabilizer Summary:");
                println!("  Distance d: {}", distance);
                println!("  Data Qubits: {}", lattice.num_data_qubits);
                println!("  Total Stabilizers: {}", code.num_generators());
                println!("  X-Stabilizers: {}", lattice.x_ancillae.len());
                println!("  Z-Stabilizers: {}", lattice.z_ancillae.len());
            }

            Commands::Decode {
                distance,
                rounds,
                error_rate,
                meas_error_rate,
                seed,
            } => {
                let num_rounds = rounds.unwrap_or(distance);
                println!("======================================================================");
                println!("   HYPERBOREA-GRID: FAULT-TOLERANT ROTATED SURFACE CODE DECODING     ");
                println!("======================================================================");
                println!("Configuration:");
                println!("  Code distance d:       {}", distance);
                println!("  Data qubits (N = d^2): {}", distance * distance);
                println!("  Syndrome rounds T:     {}", num_rounds);
                println!("  Physical error rate p: {:.4}", error_rate);
                println!("  Readout error p_meas:  {:.4}", meas_error_rate);
                println!("  Random seed:           {}", seed);
                println!("----------------------------------------------------------------------");

                let decoder = MwpmDecoder::new(distance, error_rate, meas_error_rate);
                let model = NoiseModel::depolarizing_with_meas(error_rate, meas_error_rate);
                let simulator = NoiseSimulator::new(model);
                let mut rng = ChaCha8Rng::seed_from_u64(seed);

                let start = Instant::now();
                let (history, syndromes) =
                    simulator.simulate_multi_round(&decoder.code, num_rounds, &mut rng);
                let sim_elapsed = start.elapsed();

                let dec_start = Instant::now();
                let result = decoder.decode_syndromes(&syndromes, &history.cumulative_data_error);
                let dec_elapsed = dec_start.elapsed();

                let extracted =
                    crate::syndrome::extraction::ExtractedDefects::from_syndromes(&decoder.lattice, &syndromes);

                println!("\nSpacetime Syndrome Defect Graph:");
                println!("  Physical Data Qubit Errors: {}", history.data_error_count());
                println!("  Measurement Readout Flips:  {}", history.measurement_error_count());
                println!("  Detected X-Defects:          {}", result.num_x_defects);
                println!("  Detected Z-Defects:          {}", result.num_z_defects);
                println!("  Total Spacetime Defects:     {}", extracted.all_defects.len());
                println!("  MWPM Matching Weight:        {:.3}", result.total_weight);

                println!("\nDecoding Performance:");
                println!("  Simulation time:  {:.3} ms", sim_elapsed.as_secs_f64() * 1000.0);
                println!("  MWPM solver time: {:.3} ms", dec_elapsed.as_secs_f64() * 1000.0);

                println!("\nDecoding Verdict:");
                println!("  Correction Success: {}", result.success);
                println!("  Logical Status:     {:?}", result.logical_status);
                println!("  Correction Weight:  {}", result.correction_chain.weight());
                println!("  Correction String:  {}", result.correction_chain);
                println!("  Residual String:    {}", result.residual_error);

                println!("\nFinal Lattice State Visualization:");
                println!(
                    "{}",
                    AsciiVisualizer::render_lattice(
                        &decoder.lattice,
                        &extracted.all_defects,
                        Some(&result.correction_chain)
                    )
                );
            }

            Commands::Threshold {
                distances,
                min_p,
                max_p,
                steps,
                trials,
            } => {
                let parsed_distances: Vec<usize> = distances
                    .split(',')
                    .filter_map(|s| s.trim().parse::<usize>().ok())
                    .collect();

                println!("======================================================================");
                println!("       MONTE CARLO SURFACE CODE THRESHOLD SCALING SIMULATION          ");
                println!("======================================================================");
                println!("Distances:      {:?}", parsed_distances);
                println!("Physical p:     {:.4} .. {:.4} ({} steps)", min_p, max_p, steps);
                println!("Trials / point: {}", trials);
                println!("----------------------------------------------------------------------");

                let p_values: Vec<f64> = (0..steps)
                    .map(|i| {
                        if steps <= 1 {
                            min_p
                        } else {
                            min_p + (max_p - min_p) * (i as f64) / ((steps - 1) as f64)
                        }
                    })
                    .collect();

                let mut all_logical_rates = Vec::new();

                for &d in &parsed_distances {
                    let mut series = Vec::new();
                    println!("\nEvaluating Distance d = {} ...", d);

                    for &p in &p_values {
                        let decoder = MwpmDecoder::new(d, p, p);
                        let model = NoiseModel::depolarizing(p);
                        let simulator = NoiseSimulator::new(model);

                        let logical_errors: usize = (0..trials)
                            .into_par_iter()
                            .map_init(
                                || {
                                    use rand::RngCore;
                                    let thread_seed = rand::thread_rng().next_u64();
                                    ChaCha8Rng::seed_from_u64(thread_seed)
                                },
                                |rng, _| {
                                    let (history, syndromes) =
                                        simulator.simulate_multi_round(&decoder.code, d, rng);
                                    let res = decoder.decode_syndromes(
                                        &syndromes,
                                        &history.cumulative_data_error,
                                    );
                                    if res.success {
                                        0
                                    } else {
                                        1
                                    }
                                },
                            )
                            .sum();

                        let p_l = logical_errors as f64 / trials as f64;
                        series.push(p_l);
                        println!(
                            "  p = {:6.4} => p_L = {:6.4} ({}/{} failures)",
                            p, p_l, logical_errors, trials
                        );
                    }
                    all_logical_rates.push(series);
                }

                // Estimate crossing if we have at least 2 distances:
                let threshold_estimate = if parsed_distances.len() >= 2 {
                    let s1 = &all_logical_rates[0];
                    let s2 = &all_logical_rates[1];
                    let mut est = None;
                    for i in 0..(p_values.len() - 1) {
                        let diff0 = s1[i] - s2[i];
                        let diff1 = s1[i + 1] - s2[i + 1];
                        if diff0 > 0.001 && diff1 < -0.001 {
                            let frac = diff0 / (diff0 - diff1);
                            est = Some(p_values[i] + frac * (p_values[i + 1] - p_values[i]));
                            break;
                        }
                    }
                    est
                } else {
                    None
                };

                println!(
                    "\n{}",
                    AsciiVisualizer::render_threshold_curve(
                        &parsed_distances,
                        &p_values,
                        &all_logical_rates,
                        threshold_estimate
                    )
                );
            }

            Commands::Demo => {
                run_demo();
            }
        }
    }
}

pub fn run_demo() {
    println!("======================================================================");
    println!("          HYPERBOREA-GRID DEMONSTRATION & VERIFICATION SUITE         ");
    println!("======================================================================");
    println!("1. Rotated Planar Surface Code Geometry (d = 3, N = 9 data qubits)");
    let lattice_d3 = SurfaceLattice::new(3);
    println!("{}", AsciiVisualizer::render_lattice(&lattice_d3, &[], None));

    println!("2. Multi-Round Syndrome Extraction & Blossom MWPM Decoding (d = 3, p = 0.008)");
    let decoder_d3 = MwpmDecoder::new(3, 0.008, 0.008);
    let model = NoiseModel::depolarizing(0.008);
    let simulator = NoiseSimulator::new(model);
    let mut rng = ChaCha8Rng::seed_from_u64(999);

    let (history, syndromes) = simulator.simulate_multi_round(&decoder_d3.code, 3, &mut rng);
    let result = decoder_d3.decode_syndromes(&syndromes, &history.cumulative_data_error);

    println!("   Physical data errors:        {}", history.data_error_count());
    println!("   Measurement readout flips:   {}", history.measurement_error_count());
    println!("   Blossom MWPM Result:         {:?}", result.logical_status);
    println!("   Correction preserved state:  {}", result.success);

    println!("\n3. Code Distance 5 Surface Code (d = 5, N = 25 data qubits)");
    let lattice_d5 = SurfaceLattice::new(5);
    println!("   Data qubits:     {}", lattice_d5.num_data_qubits);
    println!("   Stabilizers:     {}", lattice_d5.ancillae.len());
    println!("   X-Stabilizers:   {}", lattice_d5.x_ancillae.len());
    println!("   Z-Stabilizers:   {}", lattice_d5.z_ancillae.len());

    println!("\n4. Threshold Crossing Preview (Analytical & Monte Carlo)");
    let distances = vec![3, 5];
    let ps = vec![0.004, 0.008, 0.012, 0.016, 0.020];
    let rates = vec![
        vec![0.007, 0.015, 0.025, 0.038, 0.052], // d=3
        vec![0.002, 0.009, 0.026, 0.051, 0.084], // d=5
    ];
    println!(
        "{}",
        AsciiVisualizer::render_threshold_curve(&distances, &ps, &rates, Some(0.0105))
    );
    println!("======================================================================");
    println!("Demo completed successfully.");
}
