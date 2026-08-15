# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims
to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- A density-matrix backend (`density::DensityMatrix`, `Circuit::run_density`,
  `qsimulator --density`). Where the state vector samples one trajectory
  through noise, this carries the whole mixture, so channels, measurement and
  reset are applied *exactly* — no shots, no sampling error — and `purity`
  measures directly how much noise has cost. The trade is memory: `4^n` entries
  against `2^n`, so the ceiling is 12 qubits rather than 30.
  `DensityError::ClassicalFeedForward` refuses `if`, which needs a classical
  outcome the matrix does not carry.
- A fifth cross-validation phase comparing density matrices against Qiskit
  entrywise. Both sides are exact, so it checks coherences rather than a
  distribution, and agrees to about 1e-16.
- A noise example (`examples/noise.rs`), checking each channel against its
  analytic rate: depolarizing breaking a Bell pair's correlation at `2p/3`,
  amplitude damping decaying as `(1-gamma)^k` over rounds, and phase damping
  destroying interference without moving populations. CI runs it.

### Fixed
- `noise::depolarizing` documented the wrong convention: it said `p` was the
  probability of replacing the qubit with the maximally mixed state, when the
  Kraus operators apply X, Y or Z each with `p/3`. Those differ by a factor of
  4/3 — the state is maximally mixed at `p = 3/4`, not `p = 1`. Only the doc
  was wrong; the operators were always the ones cross-validated against Qiskit.
- `docs/design.md` listed neither `noise.rs` nor `expr.rs` in the module
  layout, and its status section predated mid-circuit measurement.

## [0.5.0] - 2026-08-13

Noise. A channel would normally need a density matrix, squaring memory from
`2^n` amplitudes to `2^n x 2^n`; simulating it by quantum trajectories instead
keeps the state vector and reuses the per-shot re-run that measurement already
required. The trade is statistical: one run is a single trajectory, and results
converge as `1/sqrt(shots)`.

### Added
- Noise channels, simulated by the quantum-trajectory method: each shot samples
  one Kraus operator rather than carrying a density matrix, which fits the
  existing per-shot re-run and keeps the state vector. New `noise` module with
  `depolarizing`, `bit_flip`, `phase_flip`, `amplitude_damping` and
  `phase_damping`, plus `is_trace_preserving`; `State::apply_kraus`;
  `Circuit::channel` and a builder per channel; the same names as program
  instructions. A single `run` is one trajectory, so noisy results need shots.
- A fourth cross-validation phase comparing noisy circuits against Qiskit Aer.
  A channel has no OpenQASM form, so it builds the Qiskit circuit from the same
  op list rather than through a shared source file.

### Changed
- **Breaking:** `ExportError` gains a `Noise` variant, since OpenQASM 2 has no
  syntax for a channel and a noisy circuit cannot be written out without
  silently dropping the noise. The enum is now `#[non_exhaustive]`, so later
  variants will be additive — but code matching it exhaustively needs a
  catch-all arm.
- The teleportation example is built with `Circuit` rather than the low-level
  `State` API: the protocol needed classical feed-forward, which now exists.
  It walks seeds until all four Bell outcomes have come up and checks each, so
  a correction that is wrong for one branch cannot hide behind the others.
- A bundled `programs/teleport.qsim` shows the same protocol in the text
  format, so feed-forward is reachable from the CLI.

## [0.4.0] - 2026-08-13

A circuit is no longer a sequence of unitaries. Mid-circuit measurement,
`reset` and classical feed-forward (`if`) all execute, which completes the
OpenQASM 2 surface apart from `opaque` — a declaration with no body to
simulate. The headline fix is that `measure` used to be accepted and silently
ignored, so a measurement with gates after it returned a confidently wrong
answer.

### Added
- Classical feed-forward: `Circuit::if_classical_eq`, the `if VALUE
  INSTRUCTION` program instruction, and OpenQASM `if (c == value) <gate>;`
  on import and export, plus `Circuit::measure_into`
  and `Circuit::n_qubits`. The circuit has one classical register as wide as
  the quantum one; `measure` writes qubit `i` to bit `i`, and declared `creg`s
  share that bit space in declaration order. `if` was the last unsupported
  OpenQASM 2 feature apart from `opaque`, which declares no body to simulate.
- `reset`: `Circuit::reset`, the `reset Q` program instruction, OpenQASM
  `reset q[i];` and `reset q;` on import, and `reset` lines on export. It was
  previously rejected as unsupported; implementing measurement made it a
  collapse plus a conditional flip. Unlike a trailing measurement, a reset
  always applies — it changes the state that gets sampled even as the last
  operation.
- A third cross-validation phase: circuits with a mid-circuit measurement,
  sampled by qsimulator and by Qiskit Aer over the exported program and
  compared by total variation distance. Such a circuit has no single state
  vector, so the two existing phases could not reach it.
- Mid-circuit measurement: `Circuit::measure`, the `measure Q` program
  instruction, OpenQASM `measure q[i] -> c[j]` and `measure q -> c` on import,
  and a `creg` plus `measure` lines on export. `Circuit::run_seeded` chooses
  the collapse stream; `Circuit::sample` re-runs the circuit per shot when a
  measurement can affect what follows it. Validated against Qiskit Aer, since a
  branching circuit cannot be checked by comparing state vectors.

### Fixed
- CI actions moved off the deprecated Node 20 runtime: `actions/checkout` v4
  to v7, `actions/cache` v4 to v6, `actions/setup-python` v5 to v7.
- OpenQASM `measure` was accepted and silently ignored, so a measurement with
  gates after it gave the wrong answer rather than an error — `h; measure; h`
  reported |0> with certainty where the real result is a coin flip, because
  ignoring the collapse turns the two Hadamards into an identity.

### Changed
- Measurements at the *end* of a circuit are treated as readout and are not
  applied by `run`, which keeps the prepared state visible and the result
  deterministic — nearly every written-out program ends in a measurement, and
  collapsing it would report one arbitrary branch instead. Sampling draws the
  same distribution either way, so this changes no histogram.

## [0.3.0] - 2026-08-09

The OpenQASM bridge now works in both directions. 0.2.0 made every builder gate
exportable; this release makes every gate Qiskit emits importable — `gate`
declarations, the whole of `qelib1`, and the arithmetic angle expressions both
of those need. Purely additive: nothing removed or changed in meaning.

### Added
- The rest of `qelib1.inc` on import: `u`, `u0`, `sx`, `sxdg`, `crx`, `cry`,
  `csx`, `cu`, `rxx`, `rzz`, `rccx`, `rc3x`, `c3x`, `c3sqrtx`, `c4x`. Together
  with `gate` declarations this means a file Qiskit exports imports and agrees
  with it. Thirteen match Qiskit's unitary exactly; `rxx` and `rzz` follow
  qelib1's own decomposition and differ from Qiskit's gate object by a global
  phase of `theta/2`, which is unobservable and inexpressible in OpenQASM 2.
- `sx` and `sxdg` gates (`Circuit::sx`/`sxdg`, `gates::sx`/`sxdg`), with the
  `sx`/`sxdg` program instructions and OpenQASM import and export. The matrix
  already existed internally as the `V` of the multi-controlled decomposition.
- OpenQASM `gate name(params) qargs { ... }` declarations are now imported,
  expanded at each call site with the actual angles and qubits substituted in.
  Bodies may call other declarations and use expressions over the formal
  parameters. This is what Qiskit emits for any composite gate, so files it
  exports now import. A file's own declaration shadows the built-in of the same
  name. Expansion is bounded by a nesting depth and a total-gate budget, so a
  short file cannot describe an unbounded one.
- The OpenQASM 2 primitives `U(theta,phi,lambda)` and `CX`.
- Angles in both front ends are now arithmetic expressions rather than a
  literal or a multiple of pi: `+ - * / ^`, parentheses, unary sign, and
  `sin`/`cos`/`tan`/`exp`/`ln`/`sqrt`, over numbers, `pi`, and (inside a
  `gate` body) the declaration's formal parameters. Every previously accepted
  form still parses.

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

[Unreleased]: https://github.com/samhithak13/qsimulator/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/samhithak13/qsimulator/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/samhithak13/qsimulator/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/samhithak13/qsimulator/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/samhithak13/qsimulator/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/samhithak13/qsimulator/releases/tag/v0.1.0
