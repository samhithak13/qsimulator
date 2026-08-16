//! # qsimulator
//!
//! A quantum circuit simulator: gates, mid-circuit measurement and reset,
//! classical feed-forward, and noise.
//!
//! There are two backends, and they check each other.
//! [`state::State`] carries a `2^n` state vector, reaching about 30 qubits, and
//! simulates a noise channel by sampling one Kraus operator per shot — so noisy
//! results converge as `1/sqrt(shots)`. [`density::DensityMatrix`] carries the
//! `2^n x 2^n` mixture instead, which makes channels, measurement and
//! feed-forward *exact* at the cost of reaching about half as far.
//!
//! The crate exposes these main pieces:
//! - [`circuit::Circuit`] — a builder that sequences operations and runs them,
//!   via [`run`](circuit::Circuit::run), [`sample`](circuit::Circuit::sample)
//!   or [`run_density`](circuit::Circuit::run_density).
//! - [`state::State`] — the `2^n` complex amplitude vector of an `n`-qubit
//!   register.
//! - [`density::DensityMatrix`] — the density-matrix backend, for exact noise.
//! - [`gates`] — unitary matrices for the standard gate set.
//! - [`noise`] — the standard single-qubit channels, as Kraus operators.
//! - [`rng::Rng`] — a seedable RNG, so sampling is reproducible.
//! - [`program`], [`qasm`] and [`qasm3`] — the text, OpenQASM 2 and OpenQASM 3
//!   front ends.
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
//!
//! The same circuit under noise, exactly rather than sampled:
//!
//! ```
//! use qsimulator::Circuit;
//!
//! let mut c = Circuit::new(2);
//! c.h(0).cnot(0, 1).depolarizing(0.1, 0);
//! let rho = c.run_density().unwrap();
//!
//! // Still normalized, but no longer a pure state.
//! assert!((rho.trace() - 1.0).abs() < 1e-12);
//! assert!(rho.purity() < 1.0);
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
