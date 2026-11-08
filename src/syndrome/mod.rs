//! Syndrome extraction, spacetime defect event identification, and matching graph construction.

pub mod extraction;
pub mod graph;
pub mod weights;

pub use extraction::{DefectEvent, ExtractedDefects};
pub use graph::SpacetimeDefectGraph;
pub use weights::MetricWeights;
