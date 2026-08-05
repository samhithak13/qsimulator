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
- CI on Linux and macOS: formatting, clippy (`-D warnings`), build, test,
  documentation (`-D warnings`), and example execution.

[Unreleased]: https://github.com/samhithak13/qsimulator/commits/main
