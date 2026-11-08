//! Pauli group operators and multi-qubit Pauli strings.
//!
//! Provides single-qubit Pauli operators $\{I, X, Y, Z\}$ and $N$-qubit Pauli strings
//! with symplectic bitvector representation, commutation testing, operator multiplication,
//! and phase tracking.

use serde::{Deserialize, Serialize};
use std::fmt;

/// Single-qubit Pauli operator.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Pauli {
    /// Identity operator $I = \begin{pmatrix} 1 & 0 \\ 0 & 1 \end{pmatrix}$
    I,
    /// Pauli-X bit-flip operator $X = \begin{pmatrix} 0 & 1 \\ 1 & 0 \end{pmatrix}$
    X,
    /// Pauli-Y bit-and-phase-flip operator $Y = \begin{pmatrix} 0 & -i \\ i & 0 \end{pmatrix}$
    Y,
    /// Pauli-Z phase-flip operator $Z = \begin{pmatrix} 1 & 0 \\ 0 & -1 \end{pmatrix}$
    Z,
}

impl Pauli {
    /// Checks if this operator is identity $I$.
    #[inline]
    pub fn is_identity(&self) -> bool {
        matches!(self, Pauli::I)
    }

    /// Checks if this operator has an X component ($X$ or $Y$).
    #[inline]
    pub fn has_x(&self) -> bool {
        matches!(self, Pauli::X | Pauli::Y)
    }

    /// Checks if this operator has a Z component ($Z$ or $Y$).
    #[inline]
    pub fn has_z(&self) -> bool {
        matches!(self, Pauli::Z | Pauli::Y)
    }

    /// Checks if two single-qubit Pauli operators commute.
    #[inline]
    pub fn commutes_with(&self, other: Pauli) -> bool {
        if self.is_identity() || other.is_identity() || *self == other {
            true
        } else {
            false
        }
    }

    /// Multiplies two Pauli operators, returning `(phase_power_of_i, result_pauli)`.
    /// Phase is represented as an integer $k \in \{0, 1, 2, 3\}$ where phase $= i^k$.
    pub fn multiply(&self, other: Pauli) -> (u8, Pauli) {
        match (*self, other) {
            (Pauli::I, p) | (p, Pauli::I) => (0, p),
            (Pauli::X, Pauli::X) => (0, Pauli::I),
            (Pauli::Y, Pauli::Y) => (0, Pauli::I),
            (Pauli::Z, Pauli::Z) => (0, Pauli::I),
            (Pauli::X, Pauli::Y) => (1, Pauli::Z), // X Y = i Z
            (Pauli::Y, Pauli::X) => (3, Pauli::Z), // Y X = -i Z = i^3 Z
            (Pauli::Y, Pauli::Z) => (1, Pauli::X), // Y Z = i X
            (Pauli::Z, Pauli::Y) => (3, Pauli::X), // Z Y = -i X = i^3 X
            (Pauli::Z, Pauli::X) => (1, Pauli::Y), // Z X = i Y
            (Pauli::X, Pauli::Z) => (3, Pauli::Y), // X Z = -i Y = i^3 Y
        }
    }

    /// Single character representation: 'I', 'X', 'Y', 'Z'.
    pub fn as_char(&self) -> char {
        match self {
            Pauli::I => 'I',
            Pauli::X => 'X',
            Pauli::Y => 'Y',
            Pauli::Z => 'Z',
        }
    }
}

impl fmt::Display for Pauli {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.as_char())
    }
}

/// An $N$-qubit Pauli string representing an element of the $N$-qubit Pauli group $\mathcal{P}_N$.
///
/// Internally stored in symplectic representation $(\mathbf{v}_X, \mathbf{v}_Z) \in \{0, 1\}^N \times \{0, 1\}^N$
/// with global phase $i^k, k \in \{0, 1, 2, 3\}$.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct PauliString {
    num_qubits: usize,
    x_bits: Vec<bool>,
    z_bits: Vec<bool>,
    phase_power: u8,
}

impl PauliString {
    /// Creates an identity Pauli string on $N$ qubits.
    pub fn identity(num_qubits: usize) -> Self {
        Self {
            num_qubits,
            x_bits: vec![false; num_qubits],
            z_bits: vec![false; num_qubits],
            phase_power: 0,
        }
    }

    /// Creates a single-qubit Pauli operator embedded on qubit `qubit_idx` out of $N$ qubits.
    pub fn single_qubit(num_qubits: usize, qubit_idx: usize, pauli: Pauli) -> Self {
        assert!(qubit_idx < num_qubits, "Qubit index out of range");
        let mut x_bits = vec![false; num_qubits];
        let mut z_bits = vec![false; num_qubits];
        if pauli.has_x() {
            x_bits[qubit_idx] = true;
        }
        if pauli.has_z() {
            z_bits[qubit_idx] = true;
        }
        Self {
            num_qubits,
            x_bits,
            z_bits,
            phase_power: 0,
        }
    }

    /// Creates a Pauli string from an array of single-qubit Pauli operators.
    pub fn from_paulis(paulis: &[Pauli]) -> Self {
        let num_qubits = paulis.len();
        let mut x_bits = vec![false; num_qubits];
        let mut z_bits = vec![false; num_qubits];
        for (i, &p) in paulis.iter().enumerate() {
            if p.has_x() {
                x_bits[i] = true;
            }
            if p.has_z() {
                z_bits[i] = true;
            }
        }
        Self {
            num_qubits,
            x_bits,
            z_bits,
            phase_power: 0,
        }
    }

    /// Creates a multi-qubit Pauli operator with given Pauli type on a subset of qubits.
    pub fn from_support(num_qubits: usize, support: &[usize], pauli: Pauli) -> Self {
        let mut x_bits = vec![false; num_qubits];
        let mut z_bits = vec![false; num_qubits];
        for &q in support {
            assert!(q < num_qubits, "Qubit index {} out of range", q);
            if pauli.has_x() {
                x_bits[q] = true;
            }
            if pauli.has_z() {
                z_bits[q] = true;
            }
        }
        Self {
            num_qubits,
            x_bits,
            z_bits,
            phase_power: 0,
        }
    }

    /// Number of qubits $N$.
    #[inline]
    pub fn num_qubits(&self) -> usize {
        self.num_qubits
    }

    /// Gets the Pauli operator acting on qubit `idx`.
    pub fn get(&self, idx: usize) -> Pauli {
        assert!(idx < self.num_qubits, "Qubit index out of range");
        let x = self.x_bits[idx];
        let z = self.z_bits[idx];
        match (x, z) {
            (false, false) => Pauli::I,
            (true, false) => Pauli::X,
            (true, true) => Pauli::Y,
            (false, true) => Pauli::Z,
        }
    }

    /// Sets the Pauli operator acting on qubit `idx`.
    pub fn set(&mut self, idx: usize, pauli: Pauli) {
        assert!(idx < self.num_qubits, "Qubit index out of range");
        self.x_bits[idx] = pauli.has_x();
        self.z_bits[idx] = pauli.has_z();
    }

    /// Computes the support of this Pauli string (indices where $P_j \ne I$).
    pub fn support(&self) -> Vec<usize> {
        (0..self.num_qubits)
            .filter(|&i| self.x_bits[i] || self.z_bits[i])
            .collect()
    }

    /// Computes the weight (number of non-identity single-qubit terms).
    pub fn weight(&self) -> usize {
        (0..self.num_qubits)
            .filter(|&i| self.x_bits[i] || self.z_bits[i])
            .count()
    }

    /// Checks whether this operator commutes with another Pauli string.
    ///
    /// Two Pauli operators commute if and only if their symplectic inner product is zero:
    /// $\langle P_1, P_2 \rangle_S = \mathbf{v}_{X,1} \cdot \mathbf{v}_{Z,2} + \mathbf{v}_{Z,1} \cdot \mathbf{v}_{X,2} \pmod 2 = 0$.
    pub fn commutes_with(&self, other: &PauliString) -> bool {
        assert_eq!(self.num_qubits, other.num_qubits, "Qubit count mismatch");
        let mut inner_product: usize = 0;
        for i in 0..self.num_qubits {
            if (self.x_bits[i] && other.z_bits[i]) ^ (self.z_bits[i] && other.x_bits[i]) {
                inner_product ^= 1;
            }
        }
        inner_product == 0
    }

    /// Checks whether this operator anticommutes with another Pauli string.
    pub fn anticommutes_with(&self, other: &PauliString) -> bool {
        !self.commutes_with(other)
    }

    /// Multiplies this Pauli string with another Pauli string: $P_{\text{new}} = P_{\text{self}} \cdot P_{\text{other}}$.
    pub fn multiply(&self, other: &PauliString) -> PauliString {
        assert_eq!(self.num_qubits, other.num_qubits, "Qubit count mismatch");
        let mut res_x = vec![false; self.num_qubits];
        let mut res_z = vec![false; self.num_qubits];
        let mut total_phase = (self.phase_power + other.phase_power) % 4;

        for i in 0..self.num_qubits {
            let p1 = self.get(i);
            let p2 = other.get(i);
            let (phase, p_res) = p1.multiply(p2);
            total_phase = (total_phase + phase) % 4;
            res_x[i] = p_res.has_x();
            res_z[i] = p_res.has_z();
        }

        PauliString {
            num_qubits: self.num_qubits,
            x_bits: res_x,
            z_bits: res_z,
            phase_power: total_phase,
        }
    }

    /// Multiplies in-place by another Pauli string: $\text{self} \leftarrow \text{self} \cdot P_{\text{other}}$.
    pub fn multiply_in_place(&mut self, other: &PauliString) {
        *self = self.multiply(other);
    }

    /// Returns a compact string representation like `"X..Z."` where '.' denotes identity.
    pub fn to_compact_string(&self) -> String {
        self.x_bits
            .iter()
            .zip(self.z_bits.iter())
            .map(|(&x, &z)| match (x, z) {
                (false, false) => '.',
                (true, false) => 'X',
                (true, true) => 'Y',
                (false, true) => 'Z',
            })
            .collect()
    }
}

impl fmt::Display for PauliString {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_compact_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_qubit_pauli_algebra() {
        assert_eq!(Pauli::I.multiply(Pauli::X), (0, Pauli::X));
        assert_eq!(Pauli::X.multiply(Pauli::X), (0, Pauli::I));
        assert_eq!(Pauli::X.multiply(Pauli::Y), (1, Pauli::Z)); // X Y = i Z
        assert_eq!(Pauli::Y.multiply(Pauli::X), (3, Pauli::Z)); // Y X = -i Z = i^3 Z
        assert_eq!(Pauli::Y.multiply(Pauli::Z), (1, Pauli::X)); // Y Z = i X
        assert_eq!(Pauli::Z.multiply(Pauli::X), (1, Pauli::Y)); // Z X = i Y

        assert!(Pauli::X.commutes_with(Pauli::X));
        assert!(!Pauli::X.commutes_with(Pauli::Z));
        assert!(Pauli::I.commutes_with(Pauli::Z));
    }

    #[test]
    fn test_pauli_string_commutation() {
        // P1 = X_0 X_1
        let p1 = PauliString::from_support(4, &[0, 1], Pauli::X);
        // P2 = Z_0 Z_1
        let p2 = PauliString::from_support(4, &[0, 1], Pauli::Z);
        // P1 and P2 share two qubits (0 and 1) where X and Z anticommute.
        // (-1)^2 = +1 => they should commute!
        assert!(p1.commutes_with(&p2));

        // P3 = Z_0 (shares 1 qubit) => should anticommute
        let p3 = PauliString::single_qubit(4, 0, Pauli::Z);
        assert!(p1.anticommutes_with(&p3));

        // P4 = Z_2 (shares 0 qubits) => should commute
        let p4 = PauliString::single_qubit(4, 2, Pauli::Z);
        assert!(p1.commutes_with(&p4));
    }

    #[test]
    fn test_pauli_string_multiplication() {
        let p1 = PauliString::single_qubit(3, 0, Pauli::X);
        let p2 = PauliString::single_qubit(3, 0, Pauli::X);
        let prod = p1.multiply(&p2);
        assert_eq!(prod.get(0), Pauli::I);
        assert_eq!(prod.weight(), 0);
    }
}
