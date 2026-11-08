//! Minimum-Weight Perfect Matching (MWPM) Blossom decoders and correction reconstructors.

pub mod blossom;
pub mod boundary;
pub mod correction;
pub mod dijkstra;
pub mod mwpm;

pub use blossom::{BlossomSolver, MatchingResult};
pub use boundary::{BoundaryResolver, MatchedPair};
pub use correction::{CorrectionReconstructor, CorrectionResult};
pub use dijkstra::DijkstraSpacetime;
pub use mwpm::MwpmDecoder;
