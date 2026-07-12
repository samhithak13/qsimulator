//! # qsimulator
//!
//! A minimal state-vector quantum circuit simulator.
//!
//! The crate exposes three main pieces:
//! - [`state::State`] — the `2^n` complex amplitude vector of an `n`-qubit register.
//! - [`gates`] — unitary matrices for the standard gate set.
//! - [`circuit::Circuit`] — a builder that sequences gates and runs them.

pub mod circuit;
pub mod gates;
pub mod state;

pub use circuit::Circuit;
pub use state::State;
