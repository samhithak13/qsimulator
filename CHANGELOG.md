# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims
to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- State-vector core: register construction, gate application (single,
  controlled, and multi-controlled), and qubit swap.
- Gate set: X, Y, Z, H, S, T, S†, T†, the phase gate P(λ), the general
  single-qubit U2/U3, Rx/Ry/Rz, SWAP, CNOT, CZ, controlled-Rz, controlled
  phase, Toffoli, and arbitrary multi-controlled unitaries.
- Measurement: per-qubit and full-register collapse, plus Born-rule sampling
  driven by a seedable, dependency-free RNG.
- Circuit builder with ASCII diagram rendering (`Circuit::diagram` / `Display`).
- Front ends: a line-based text program format and an OpenQASM 2.0 subset
  importer and exporter, with a CLI over both.
- Worked examples (GHZ, Grover, teleportation) and a broad oracle-based test
  suite.
- `--statevector` CLI flag that prints the final amplitudes as JSON.
- Cross-validation against Qiskit (`crossval/`): random circuits over the
  shared gate set match Qiskit's state vectors to floating-point precision.
- CI on Linux and macOS: formatting, clippy (`-D warnings`), build, test,
  documentation (`-D warnings`), example execution, and the Qiskit
  cross-validation.
- `cargo bench` throughput harness (`benches/throughput.rs`).
- Optional `parallel` feature: a rayon-backed multithreaded `apply_1q`
  (~1.5–2.4x on an Apple M1), off by default so the core stays dependency-free.

### Fixed
- The OpenQASM importer no longer panics on malformed bracket order (a `qreg`
  or qubit reference like `q]0[`); it now returns an error.
- Both parsers reject a register larger than 30 qubits instead of aborting the
  process on allocation.

### Changed
- The parsers and exporter now return typed errors implementing
  `std::error::Error` (`ParseError`, `ExportError`) instead of `String`, so the
  library composes with `?` and `Box<dyn Error>`.
- Restructured `apply_1q` into a bounds-check-free walk over `chunks_exact_mut`
  and hoisted gate-matrix entries into locals across the gate kernels; roughly
  1.3–1.6× higher single-gate throughput on an Apple M1 (see benches/README.md).

[Unreleased]: https://github.com/samhithak13/qsimulator/commits/main
