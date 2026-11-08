//! Boundary matching resolution for virtual augmented spacetime syndrome graphs.
//!
//! Maps Blossom perfect matching outcomes on the $2m$-vertex graph back to physical
//! defect-to-defect pairs and defect-to-boundary terminations.

use crate::code::lattice::Boundary;
use crate::syndrome::extraction::DefectEvent;
use crate::decoder::blossom::MatchingResult;
use serde::{Deserialize, Serialize};

/// A physical pairing outcome determined by the decoder.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum MatchedPair {
    /// Two physical defect events matched to each other.
    DefectPair {
        /// First defect event.
        defect1: DefectEvent,
        /// Second defect event.
        defect2: DefectEvent,
    },
    /// A physical defect matched to its code boundary.
    BoundaryPair {
        /// The defect event matched to boundary.
        defect: DefectEvent,
        /// Which boundary it is terminated at.
        boundary: Option<Boundary>,
    },
}

/// Helper for resolving augmented matching into physical pairings.
pub struct BoundaryResolver;

impl BoundaryResolver {
    /// Resolves Blossom matching pairs on the $2m$-vertex graph.
    ///
    /// - $0 \le i < m$: Real defect events.
    /// - $m \le i < 2m$: Virtual boundary companion vertices.
    pub fn resolve_pairs(
        defects: &[DefectEvent],
        matching: &MatchingResult,
    ) -> Vec<MatchedPair> {
        let m = defects.len();
        let mut physical_pairs = Vec::new();

        for &(u, v) in &matching.pairs {
            if u < m && v < m {
                // Defect-to-defect pairing:
                physical_pairs.push(MatchedPair::DefectPair {
                    defect1: defects[u].clone(),
                    defect2: defects[v].clone(),
                });
            } else if u < m && v >= m {
                // Defect u matched to virtual boundary v:
                physical_pairs.push(MatchedPair::BoundaryPair {
                    defect: defects[u].clone(),
                    boundary: None,
                });
            } else if u >= m && v < m {
                // Defect v matched to virtual boundary u:
                physical_pairs.push(MatchedPair::BoundaryPair {
                    defect: defects[v].clone(),
                    boundary: None,
                });
            } else {
                // Virtual boundary matched to virtual boundary (B_i, B_j): ignore (cost 0)
            }
        }

        physical_pairs
    }
}
