# Contributing

Thanks for your interest in qsimulator. This is a small, focused project; the
bar for changes is correctness and clarity.

## Development

```bash
cargo build
cargo test
cargo run                 # the Bell-state demo
cargo run -- --help       # CLI usage
```

Before opening a pull request, run what CI runs:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets -- -D warnings
cargo test
cargo test --features parallel
cargo doc --no-deps                       # with RUSTDOCFLAGS="-D warnings"
```

## Expectations

- **Every change lands with its test in the same commit.** New gates and
  primitives need both a unit test (matrix or truth table) and a circuit-level
  integration test; every bug fix needs a regression test.
- **Public items are documented.** `missing_docs` is a warning and CI denies
  warnings, so an undocumented public item fails the build.
- **No `unsafe`.** The crate is `#![forbid(unsafe_code)]`.
- **The core stays dependency-free** apart from `num-complex`. Optional
  functionality goes behind a feature (see `parallel`).

## Adding a gate

1. Add the matrix to `src/gates.rs` with a unit test.
2. Add a builder to `src/circuit.rs` (and an `Op` arm if it needs new state
   machinery).
3. If it should be importable, wire it into `src/program.rs` and `src/qasm.rs`,
   and into `Circuit::to_qasm` for export.
4. Add integration tests, and — for anything in the OpenQASM subset — confirm
   it round-trips and still matches Qiskit (`crossval/`).

## Frozen conventions

Little-endian qubit order, row-major gate storage, rotations as
`exp(-i·θ/2·P)`, and deterministic RNG output. See `docs/design.md` for the
full list; changing any of these silently breaks results.

## Verification

Correctness rests on independent oracles, not hand-copied expected values.
See `docs/design.md` (testing strategy), `crossval/` (Qiskit cross-validation),
and `fuzz/` (parser fuzzing).
