//! # qsimulator
//!
//! A noiseless state-vector quantum circuit simulator.
//!
//! The crate exposes these main pieces:
//! - [`state::State`] — the `2^n` complex amplitude vector of an `n`-qubit register.
//! - [`gates`] — unitary matrices for the standard gate set.
//! - [`circuit::Circuit`] — a builder that sequences gates and runs them.
//! - [`rng::Rng`] — a seedable RNG used for reproducible measurement sampling.
//! - [`program`] and [`qasm`] — text and OpenQASM 2.0 front ends.
//!
//! ```
//! use qsimulator::Circuit;
//!
//! let mut c = Circuit::new(2);
//! c.h(0).cnot(0, 1); // Bell state
//! let state = c.run();
//! assert!((state.probability(0b00) - 0.5).abs() < 1e-12);
//! assert!((state.probability(0b11) - 0.5).abs() < 1e-12);
//! ```
#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod circuit;
pub mod density;
pub mod error;
mod expr;
pub mod gates;
pub mod noise;
pub mod program;
pub mod qasm;
pub mod qasm3;
pub mod rng;
pub mod state;

pub use circuit::{Circuit, DensityError, ExportError};
pub use error::ParseError;
pub use rng::Rng;
pub use state::State;
