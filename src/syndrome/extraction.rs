//! Multi-round syndrome extraction and spacetime defect event identification.
//!
//! Computes defect events $\Delta s(a, t) = s(a, t) \oplus s(a, t-1)$ across consecutive measurement rounds
//! and separates them into independent $X$-type and $Z$-type syndrome sets for decoding.

use crate::code::lattice::{StabilizerType, SurfaceLattice};
use serde::{Deserialize, Serialize};

/// A detected defect event in the spacetime syndrome lattice.
///
/// A defect occurs when a stabilizer's measurement eigenvalue flips relative to the previous round:
/// $\Delta s(a, t) = s(a, t) \oplus s(a, t-1) = 1$.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct DefectEvent {
    /// Defect identifier within the decoded set.
    pub id: usize,
    /// Ancilla index $0 \le a < d^2 - 1$.
    pub ancilla_id: usize,
    /// Temporal round $0 \le t \le T$.
    pub round: usize,
    /// Stabilizer generator type ($X$ or $Z$).
    pub kind: StabilizerType,
    /// Spatial lattice coordinate $(u, v)$.
    pub coord: (i32, i32),
}

/// Extracted spacetime syndrome defect sets.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExtractedDefects {
    /// Total measurement rounds $T$ plus terminal round ($T+1$).
    pub num_rounds: usize,
    /// Defect events detected on X-stabilizers (arising from Z-errors).
    pub x_defects: Vec<DefectEvent>,
    /// Defect events detected on Z-stabilizers (arising from X-errors).
    pub z_defects: Vec<DefectEvent>,
    /// Full record of raw defect events across all stabilizers.
    pub all_defects: Vec<DefectEvent>,
}

impl ExtractedDefects {
    /// Computes defect events from multi-round syndrome measurement records.
    ///
    /// # Arguments
    /// - `lattice`: The surface code lattice structure.
    /// - `measured_syndromes`: Syndrome bitvectors for rounds $t = 0, \dots, T$.
    pub fn from_syndromes(lattice: &SurfaceLattice, measured_syndromes: &[Vec<bool>]) -> Self {
        let num_rounds = measured_syndromes.len();
        let num_ancillae = lattice.ancillae.len();

        let mut x_defects = Vec::new();
        let mut z_defects = Vec::new();
        let mut all_defects = Vec::new();

        let mut all_id = 0;
        let mut x_id = 0;
        let mut z_id = 0;

        for t in 0..num_rounds {
            let curr = &measured_syndromes[t];
            let prev = if t == 0 {
                None
            } else {
                Some(&measured_syndromes[t - 1])
            };

            for a in 0..num_ancillae {
                let s_curr = curr[a];
                let s_prev = prev.map(|p| p[a]).unwrap_or(false);
                let is_defect = s_curr ^ s_prev;

                if is_defect {
                    let anc = &lattice.ancillae[a];
                    let defect = DefectEvent {
                        id: all_id,
                        ancilla_id: a,
                        round: t,
                        kind: anc.kind,
                        coord: anc.coord,
                    };

                    match anc.kind {
                        StabilizerType::X => {
                            let mut d = defect.clone();
                            d.id = x_id;
                            x_defects.push(d);
                            x_id += 1;
                        }
                        StabilizerType::Z => {
                            let mut d = defect.clone();
                            d.id = z_id;
                            z_defects.push(d);
                            z_id += 1;
                        }
                    }

                    all_defects.push(defect);
                    all_id += 1;
                }
            }
        }

        Self {
            num_rounds,
            x_defects,
            z_defects,
            all_defects,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_round_defect_extraction() {
        let lattice = SurfaceLattice::new(3);
        let syn_round_0 = vec![false; 8];
        let mut syn_round_1 = vec![false; 8];
        syn_round_1[0] = true; // One defect on ancilla 0

        let syndromes = vec![syn_round_0, syn_round_1];
        let extracted = ExtractedDefects::from_syndromes(&lattice, &syndromes);

        assert_eq!(extracted.all_defects.len(), 1);
        assert_eq!(extracted.all_defects[0].ancilla_id, 0);
        assert_eq!(extracted.all_defects[0].round, 1);
    }
}
