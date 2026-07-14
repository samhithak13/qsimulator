# qsimulator

A **quantum circuit simulator** written in Rust.

`qsimulator` models an ideal (noiseless) quantum computer by maintaining
the full state vector of an *n*-qubit register and applying unitary gate
operations to it. It is intended as a clear, well-tested foundation for
experimenting with quantum algorithms and for learning how state-vector
simulation works under the hood.

## Scope & goals

- **State-vector simulation** of up to ~20–25 qubits on a typical laptop
  (memory grows as `2^n` complex amplitudes).
- A small, readable **gate set**: X, Y, Z, H, S, T, the rotations Rx/Ry/Rz,
  SWAP, plus controlled gates (CNOT, Toffoli) and arbitrary single-qubit
  unitaries.
- A **circuit builder** API for composing gates into programs.
- **Measurement** in the computational basis with correct Born-rule
  probabilities, post-measurement collapse, and seedable sampling.
- Correctness first: every gate and circuit primitive is covered by tests.

### Non-goals (for now)

- Noise / density-matrix simulation.
- GPU acceleration or distributed simulation.
- A hardware backend or transpiler.

These may be revisited once the core is stable — see the tracking issue.

## Project layout

```
qsimulator/
├── src/            # Library + CLI source
│   ├── lib.rs      # Crate root, re-exports
│   ├── state.rs    # State vector representation
│   ├── gates.rs    # Gate definitions (unitary matrices)
│   ├── circuit.rs  # Circuit builder & execution
│   ├── program.rs  # Text program parser for the CLI
│   └── main.rs     # CLI: runs program files, stdin, or a Bell demo
├── examples/       # Sample .qsim program files
├── tests/          # Integration tests
├── docs/           # Design notes & documentation
└── .github/        # CI workflows
```

## Quick start

```bash
# Run the built-in demo (prepares and measures a Bell state)
cargo run

# Run a program file …
cargo run -- examples/ghz.qsim

# … or pipe a program in on standard input
echo 'qubits 1
h 0
sample 1000' | cargo run -- -

# See the program grammar and all directives
cargo run -- --help

# Run the test suite
cargo test
```

### Program files

A program is one directive per line (`#` starts a comment). The first
directive must declare the register size:

```text
# A 3-qubit GHZ state, sampled 1000 times.
qubits 3
h 0
cnot 0 1
cnot 0 2
sample 1000 42
```

Supported instructions: `h|x|z <t>`, `rx|ry|rz <angle> <t>` (angle is a
float or a `pi`-expression like `pi/2` or `3pi/4`), `cnot <c> <t>`,
`swap <a> <b>`, `toffoli <c1> <c2> <t>`, and `sample <shots> [seed]`.
Running a program prints the final amplitudes and, if a `sample` directive
is present, a measurement histogram.

## Roadmap

Milestones are tracked in the repository issues. High level:

1. **v0.1 — Core** (this scaffold): state vector, single-qubit gates,
   CNOT, measurement, circuit builder.  ✅ measurement done (seedable
   sampling, single-qubit + full-register collapse).
2. **v0.2 — Ergonomics**: more gates (rotations, SWAP, Toffoli ✅), circuit
   diagram printing, richer CLI (program files + stdin ✅).
3. **v0.3 — Performance**: in-place gate application, sparse fast paths,
   benchmarks.

## Contributing

Contributions welcome. Please keep new gates/primitives covered by tests
and run `cargo fmt` + `cargo clippy` before opening a PR.

## License

Licensed under the [MIT License](./LICENSE).
