# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims
to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Controlled-U3 gate: `Circuit::cu3`, plus OpenQASM `cu3` import and export.
- Controlled-Y (`cy`), controlled-Hadamard (`ch`), and Fredkin / controlled-SWAP
  (`cswap`) gates, with OpenQASM import and export.

### Changed
- The `parallel` feature now also threads the controlled and multi-controlled
  kernels (~1.2x on an Apple M1), via a shared `apply_masked` routine.
- `Circuit::to_qasm` now exports an arbitrary controlled-U (`cu`) by
  decomposing it into a control phase (`u1`) and a `cu3`, instead of returning
  an error. `ExportError::ControlledU` is removed; only `mcu` and a
  multi-controlled-X with a control count other than two still error.

## [0.1.0] - 2026-08-06

Initial release.

### Added
- State-vector core: register construction and in-place gate application
  (single-, controlled-, and multi-controlled), with a bounds-check-free
  single-qubit kernel.
- Gate set: X, Y, Z, H, S, T, S†, T†, the phase gate P(λ), the general
  single-qubit U2/U3, Rx/Ry/Rz, SWAP, CNOT, CZ, controlled-Rz, controlled
  phase, Toffoli, and arbitrary multi-controlled unitaries.
- Measurement: per-qubit and full-register collapse, plus Born-rule sampling
  driven by a seedable, dependency-free RNG.
- Circuit builder with ASCII diagram rendering (`Circuit::diagram` / `Display`).
- Front ends: a line-based text program format and an OpenQASM 2.0 subset
  importer and exporter, with a CLI over both (`--emit-qasm`, `--statevector`).
- Typed error types (`ParseError`, `ExportError`) implementing
  `std::error::Error`.
- Optional `parallel` feature: a rayon-backed multithreaded `apply_1q`
  (~1.5–2.4x on an Apple M1), off by default so the core depends only on
  `num-complex`.
- Worked examples: GHZ, Grover, and quantum teleportation.
- Verification and tooling: a broad oracle-based test suite, Qiskit
  cross-validation (`crossval/`), parser fuzzing (`fuzz/` and a stable
  robustness test), a `cargo bench` throughput harness, and CI on Linux and
  macOS covering formatting, clippy (`-D warnings`), a warning-clean
  `cargo doc`, tests, examples, cross-validation, fuzzing, and `cargo audit`.

[Unreleased]: https://github.com/samhithak13/qsimulator/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/samhithak13/qsimulator/releases/tag/v0.1.0
