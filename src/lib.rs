//! # qsimulator
//!
//! A minimal state-vector quantum circuit simulator.
//!
//! The crate exposes these main pieces:
//! - [`state::State`] — the `2^n` complex amplitude vector of an `n`-qubit register.
//! - [`gates`] — unitary matrices for the standard gate set.
//! - [`circuit::Circuit`] — a builder that sequences gates and runs them.
//! - [`rng::Rng`] — a seedable RNG used for reproducible measurement sampling.

pub mod circuit;
pub mod gates;
pub mod program;
pub mod qasm;
pub mod rng;
pub mod state;

pub use circuit::Circuit;
pub use rng::Rng;
pub use state::State;
