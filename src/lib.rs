//! # HyperboreaGrid: Fault-Tolerant Rotated Surface Code & MWPM Blossom Decoder
//!
//! Production-grade quantum error correction (QEC) engine implementing:
//! - Rotated planar surface code geometry of arbitrary odd distance $d \ge 3$.
//! - Bulk weight-4 and boundary weight-2 $X$ and $Z$ stabilizer generators.
//! - Weight-$d$ logical operator representatives $X_L$ and $Z_L$.
//! - Phenomenological and circuit-level Pauli noise models (depolarizing, independent $X/Z$, readout flips).
//! - Multi-round fault-tolerant spacetime defect extraction.
//! - Minimum-Weight Perfect Matching (MWPM) Blossom algorithm decoder for spacetime syndrome graphs.
//! - Correction chain operator reconstruction and logical fidelity evaluation.
//! - 7-bit ASCII lattice visualizer and Monte Carlo threshold scaling curves.

pub mod code;
pub mod decoder;
pub mod noise;
pub mod syndrome;
pub mod visualizer;
pub mod cli;

pub use code::{
    AncillaDescriptor, Boundary, LogicalOperators, LogicalStatus, Pauli, PauliString,
    StabilizerCode, StabilizerGenerator, StabilizerType, SurfaceLattice,
};
pub use decoder::{
    BlossomSolver, BoundaryResolver, CorrectionReconstructor, CorrectionResult,
    DijkstraSpacetime, MatchedPair, MatchingResult, MwpmDecoder,
};
pub use noise::{ErrorEvent, ErrorHistory, NoiseModel, NoiseSimulator};
pub use syndrome::{DefectEvent, ExtractedDefects, MetricWeights, SpacetimeDefectGraph};
pub use visualizer::AsciiVisualizer;
