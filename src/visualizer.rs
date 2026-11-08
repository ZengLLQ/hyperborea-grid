//! 7-bit ASCII visualization for rotated surface code lattices, defects, and threshold curves.
//!
//! Provides clean, zero-dependency text visualization of the 2D rotated lattice geometry,
//! defect events, correction chains, and Monte Carlo threshold scaling curves using 7-bit ASCII only.

use crate::code::lattice::{StabilizerType, SurfaceLattice};
use crate::code::pauli::PauliString;
use crate::syndrome::extraction::DefectEvent;

/// ASCII Visualizer for rotated surface codes and threshold scaling curves.
pub struct AsciiVisualizer;

impl AsciiVisualizer {
    /// Renders a 7-bit ASCII map of the rotated surface code lattice showing data qubits,
    /// stabilizer generators, detected defects, and applied correction chains.
    pub fn render_lattice(
        lattice: &SurfaceLattice,
        defects: &[DefectEvent],
        corrections: Option<&PauliString>,
    ) -> String {
        let d = lattice.distance;
        let d_i32 = d as i32;
        let mut out = String::new();

        out.push_str(&format!(
            "=== Rotated Surface Code Lattice (Distance d = {}, Data Qubits = {}, Stabilizers = {}) ===\n\n",
            d, lattice.num_data_qubits, lattice.ancillae.len()
        ));

        let defect_ancillae: Vec<usize> = defects.iter().map(|d| d.ancilla_id).collect();

        // Embedding grid dimensions:
        // u runs from -1 to 2d-1
        // v runs from -1 to 2d-1
        for u in -1..=(2 * d_i32 - 1) {
            let mut row_str = String::new();

            for v in -1..=(2 * d_i32 - 1) {
                if u % 2 == 0 && v % 2 == 0 {
                    // Data qubit location (2r, 2c)
                    let r = u / 2;
                    let c = v / 2;
                    if r >= 0 && r < d_i32 && c >= 0 && c < d_i32 {
                        let q_idx = (r as usize) * d + (c as usize);
                        let corr = corrections.map(|c| c.get(q_idx)).unwrap_or(crate::code::pauli::Pauli::I);
                        match corr {
                            crate::code::pauli::Pauli::I => {
                                row_str.push_str(&format!("(d{:02})", q_idx));
                            }
                            crate::code::pauli::Pauli::X => {
                                row_str.push_str(&format!("[X{:02}]", q_idx));
                            }
                            crate::code::pauli::Pauli::Z => {
                                row_str.push_str(&format!("[Z{:02}]", q_idx));
                            }
                            crate::code::pauli::Pauli::Y => {
                                row_str.push_str(&format!("[Y{:02}]", q_idx));
                            }
                        }
                    } else {
                        row_str.push_str("     ");
                    }
                } else if let Some(&anc_idx) = lattice.coord_to_ancilla.get(&(u, v)) {
                    // Ancilla / Stabilizer location
                    let anc = &lattice.ancillae[anc_idx];
                    let is_defect = defect_ancillae.contains(&anc_idx);

                    let symbol = match (anc.kind, is_defect) {
                        (StabilizerType::X, false) => "  X  ",
                        (StabilizerType::Z, false) => "  Z  ",
                        (StabilizerType::X, true) => " *X* ",
                        (StabilizerType::Z, true) => " !Z! ",
                    };
                    row_str.push_str(symbol);
                } else {
                    row_str.push_str("     ");
                }
            }

            // Only print row if it contains non-whitespace characters
            if row_str.trim().is_empty() {
                continue;
            }

            out.push_str(&row_str);
            out.push('\n');
        }

        out.push_str("\nLegend:\n");
        out.push_str("  (d##)  : Data qubit (identity, no correction)\n");
        out.push_str("  [X##]  : Data qubit with Pauli X correction applied\n");
        out.push_str("  [Z##]  : Data qubit with Pauli Z correction applied\n");
        out.push_str("  [Y##]  : Data qubit with Pauli Y correction applied\n");
        out.push_str("    X    : Inactive X-stabilizer (weight-2 or weight-4)\n");
        out.push_str("    Z    : Inactive Z-stabilizer (weight-2 or weight-4)\n");
        out.push_str("   *X*   : Active X-defect event (triggered by phase flip)\n");
        out.push_str("   !Z!   : Active Z-defect event (triggered by bit flip)\n");

        out
    }

    /// Renders a 7-bit ASCII Monte Carlo threshold scaling plot of $p_L$ vs $p$
    /// for code distances $d=3, 5, 7$, displaying the crossing threshold near ~1%.
    pub fn render_threshold_curve(
        distances: &[usize],
        physical_ps: &[f64],
        logical_rates: &[Vec<f64>], // [distance_idx][p_idx]
        threshold_estimate: Option<f64>,
    ) -> String {
        let width = 60;
        let height = 18;
        let mut out = String::new();

        out.push_str("======================================================================\n");
        out.push_str("       MONTE CARLO SURFACE CODE THRESHOLD SCALING (7-BIT ASCII)       \n");
        out.push_str("======================================================================\n\n");

        // Find min and max values:
        let min_p = physical_ps.first().copied().unwrap_or(0.001);
        let max_p = physical_ps.last().copied().unwrap_or(0.03);

        let mut max_pl = 0.05f64;
        for series in logical_rates {
            for &val in series {
                if val > max_pl {
                    max_pl = val;
                }
            }
        }
        max_pl = (max_pl * 1.15).min(1.0);

        // Marker characters for different distances:
        let markers = ['3', '5', '7', '9'];

        let mut canvas = vec![vec![' '; width]; height];

        // Draw plot grid:
        for y in 0..height {
            canvas[y][0] = '|';
        }
        for x in 0..width {
            canvas[height - 1][x] = '-';
        }
        canvas[height - 1][0] = '+';

        // Draw vertical threshold line if available:
        if let Some(p_th) = threshold_estimate {
            if p_th >= min_p && p_th <= max_p {
                let x_th = ((p_th - min_p) / (max_p - min_p) * (width as f64 - 2.0)).round() as usize + 1;
                if x_th < width {
                    for y in 0..(height - 1) {
                        if canvas[y][x_th] == ' ' {
                            canvas[y][x_th] = ':';
                        }
                    }
                }
            }
        }

        // Plot data series:
        for (d_idx, series) in logical_rates.iter().enumerate() {
            let marker = markers[d_idx % markers.len()];

            for (&p, &pl) in physical_ps.iter().zip(series.iter()) {
                let x_frac = ((p - min_p) / (max_p - min_p)).clamp(0.0, 1.0);
                let y_frac = (pl / max_pl).clamp(0.0, 1.0);

                let x = (x_frac * (width as f64 - 2.0)).round() as usize + 1;
                let y = ((1.0 - y_frac) * (height as f64 - 2.0)).round() as usize;

                if x < width && y < height {
                    canvas[y][x] = marker;
                }
            }
        }

        // Output lines with Y-axis labels:
        for y in 0..height {
            let label = if y == 0 {
                format!("{:5.3} ", max_pl)
            } else if y == height / 2 {
                format!("{:5.3} ", max_pl / 2.0)
            } else if y == height - 1 {
                "0.000 ".to_string()
            } else {
                "      ".to_string()
            };

            let row_str: String = canvas[y].iter().collect();
            out.push_str(&label);
            out.push_str(&row_str);
            out.push('\n');
        }

        // X-axis label:
        out.push_str("       +");
        for _ in 1..width {
            out.push('-');
        }
        out.push('\n');
        out.push_str(&format!(
            "       {:5.3}{:>width$}\n",
            min_p,
            format!("{:5.3}", max_p),
            width = width - 6
        ));
        out.push_str("                     Physical Error Rate p (per operation)\n\n");

        out.push_str("Series Legend:\n");
        for (d_idx, &dist) in distances.iter().enumerate() {
            let marker = markers[d_idx % markers.len()];
            out.push_str(&format!("  '{}' : Distance d = {}\n", marker, dist));
        }
        if let Some(p_th) = threshold_estimate {
            out.push_str(&format!(
                "  ':' : Estimated Threshold Crossing p_th ~ {:.3} ({:.2}%)\n",
                p_th,
                p_th * 100.0
            ));
        }

        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_render_lattice_d3() {
        let lattice = SurfaceLattice::new(3);
        let rendered = AsciiVisualizer::render_lattice(&lattice, &[], None);
        assert!(rendered.contains("Distance d = 3"));
        assert!(rendered.contains("(d00)"));
        assert!(rendered.is_ascii());
    }

    #[test]
    fn test_render_threshold_curve_ascii() {
        let distances = vec![3, 5];
        let ps = vec![0.005, 0.010, 0.015, 0.020];
        let rates = vec![
            vec![0.008, 0.015, 0.025, 0.038], // d=3
            vec![0.002, 0.010, 0.028, 0.055], // d=5 (crosses near 0.010)
        ];

        let rendered = AsciiVisualizer::render_threshold_curve(&distances, &ps, &rates, Some(0.010));
        assert!(rendered.contains("MONTE CARLO SURFACE CODE THRESHOLD SCALING"));
        assert!(rendered.contains("Estimated Threshold Crossing"));
        assert!(rendered.is_ascii());
    }
}
