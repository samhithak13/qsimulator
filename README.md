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
- A small, readable **gate set**: X, Y, Z, H, S, T, S†, T†, the phase gate
  P(λ), the rotations Rx/Ry/Rz, SWAP, plus controlled gates (CNOT, CZ,
  Toffoli) and arbitrary controlled unitaries with any number of controls
  (`cu`, `mcx`, `mcu`).
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
│   └── main.rs     # Small demo CLI (Bell state)
├── tests/          # Integration tests
├── docs/           # Design notes & documentation
└── .github/        # CI workflows
```

## Quick start

```bash
# Run the built-in demo (prepares and samples a Bell state)
cargo run

# Run a circuit described in a text program (see programs/ghz.qsim)
cargo run -- programs/ghz.qsim

# Import and run an OpenQASM 2.0 file (see programs/bell.qasm)
cargo run -- programs/bell.qasm

# Read a program from stdin, or show the program format
printf 'qubits 1\nx 0\n' | cargo run -- -
cargo run -- --help

# Run the test suite
cargo test
```

### Program format

Circuits can be written as a small line-based text format instead of Rust —
one instruction per line, `#` for comments:

```text
qubits 3
h 0
cnot 0 1
cnot 1 2        # GHZ state
sample 1000 42
```

Supported: `qubits N`; `h|x|y|z|s|t|sdg|tdg Q`; `rx|ry|rz|p THETA Q` (THETA is
a float or a symbolic multiple of pi like `pi/2`, `-pi/4`, `2pi`);
`cnot|cz C T`; `swap A B`; `toffoli C1 C2 T`; and an optional `sample SHOTS
SEED`.

### OpenQASM 2.0 import

A file ending in `.qasm` (or starting with an `OPENQASM` header) is parsed as
an OpenQASM 2.0 subset — enough for hand-written textbook circuits:

```text
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
h q[0];
cx q[0],q[1];
```

Supported gates: `x y z h s t sdg tdg`, `rx/ry/rz(theta)`, `u1(lambda)`
(alias `p`), `cx`, `cz`, `swap`, `ccx`.
Multiple `qreg`s map into a single flat qubit space in declaration order;
`creg`, `barrier`, and `measure` are accepted and ignored; unsupported
features (custom `gate`, `if`, `reset`, …) are reported rather than
mis-simulated.

Export goes the other way — `Circuit::to_qasm()` (and `qsimulator --emit-qasm
<FILE>`) writes a circuit back out as OpenQASM 2.0. Any circuit using only the
supported gates round-trips exactly through import; gates outside the subset
(arbitrary controlled-U, C³X, …) are reported as an export error.

## Roadmap

Milestones are tracked in the repository issues. High level:

1. **v0.1 — Core** (this scaffold): state vector, single-qubit gates,
   CNOT, measurement, circuit builder.  ✅ measurement done (seedable
   sampling, single-qubit + full-register collapse).
2. **v0.2 — Ergonomics**: more gates (rotations, SWAP, Toffoli ✅), circuit
   diagram printing ✅, richer CLI ✅ (text program format).
3. **v0.3 — Interop & performance**: OpenQASM 2.0 import ✅; benchmarks,
   sparse fast paths (planned).

## Contributing

Contributions welcome. Please keep new gates/primitives covered by tests
and run `cargo fmt` + `cargo clippy` before opening a PR.

## License

Licensed under the [MIT License](./LICENSE).
