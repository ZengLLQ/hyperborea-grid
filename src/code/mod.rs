//! Rotated planar surface code geometry, stabilizer algebra, and logical operators.

pub mod lattice;
pub mod logical;
pub mod pauli;
pub mod stabilizer;

pub use lattice::{AncillaDescriptor, Boundary, StabilizerType, SurfaceLattice};
pub use logical::{LogicalOperators, LogicalStatus};
pub use pauli::{Pauli, PauliString};
pub use stabilizer::{StabilizerCode, StabilizerGenerator};
