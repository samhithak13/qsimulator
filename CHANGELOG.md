# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims
to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- OpenQASM import split statements on `;`, which glued a `gate` body's closing
  `}` onto the statement after it. Since `gate` blocks conventionally precede
  `qreg` — as Qiskit emits them — the register declaration was swallowed and
  the error read `no \`qreg\` declared` for a file that plainly had one. The
  splitter is now brace-aware, and unbalanced braces are reported as such.

## [0.2.0] - 2026-08-09

Every gate in the builder now round-trips through OpenQASM 2.0, which was the
last interop gap. Removing the export-error variants that can no longer occur
makes this a breaking change.

### Added
- Controlled-U3 gate: `Circuit::cu3`, plus OpenQASM `cu3` import and export.
- Controlled-Y (`cy`), controlled-Hadamard (`ch`), and Fredkin / controlled-SWAP
  (`cswap`) gates, with OpenQASM import and export.
- A quantum Fourier transform example (`examples/qft.rs`).
- Property tests (`tests/properties.rs`): random circuits preserve the norm,
  self-inverse gates cancel, and inverse rotations cancel.
- CLI integration tests (`tests/cli.rs`) driving the built binary, and broader
  parser coverage — raising line coverage to ~97%.
- A CI coverage job (`cargo llvm-cov`) that fails if line coverage drops below
  90%.
- An `id` (identity) gate: `Circuit::id`, `gates::id`, the `id Q` program
  instruction, and OpenQASM `id` import and export.
- CLI `--shots N` and `--seed S` flags, so any program can be sampled --
  including an OpenQASM file, which has no way to carry a `sample` directive.
  They override the program's own directive when it has one.
- Text program instructions `mcx C... T` and `mcu3 THETA PHI LAMBDA C... T`,
  which take any number of controls, so multi-controlled gates are reachable
  from the CLI.

### Changed
- The `parallel` feature now also threads the controlled and multi-controlled
  kernels (~1.2x on an Apple M1), via a shared `apply_masked` routine.
- `Circuit::to_qasm` now exports an arbitrary controlled-U (`cu`) by
  decomposing it into a control phase (`u1`) and a `cu3`, instead of returning
  an error.
- `Circuit::to_qasm` now also exports multi-controlled gates (`mcx`, `mcu`) at
  any width, via a Barenco-style decomposition into Toffolis and single-qubit
  rotations: a borrowed-qubit Toffoli ladder when the register has a spare
  qubit, and the square-root recursion when it does not. Every built-in gate
  now round-trips through OpenQASM.

### Removed
- **Breaking:** `ExportError::ControlledU`, `ExportError::MultiControlledU`,
  and `ExportError::MultiControlledX`, none of which can occur now that every
  controlled and multi-controlled gate exports. `ExportError::SingleGate` remains but is
  unreachable through the builder API. Code matching on `ExportError` needs
  those arms dropped.

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

[Unreleased]: https://github.com/samhithak13/qsimulator/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/samhithak13/qsimulator/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/samhithak13/qsimulator/releases/tag/v0.1.0
