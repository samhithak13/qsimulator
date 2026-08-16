# qsimulator

[![CI](https://github.com/samhithak13/qsimulator/actions/workflows/ci.yml/badge.svg)](https://github.com/samhithak13/qsimulator/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/samhithak13/qsimulator/blob/main/LICENSE)

A noiseless quantum circuit simulator in Rust. It holds the full `2^n` state
vector of an *n*-qubit register, applies unitary gates, measures in the
computational basis, and reads and writes OpenQASM 2.0. The priority is being
correct and readable, not fast.

## Example

```rust
use qsimulator::Circuit;

let mut c = Circuit::new(2);
c.h(0).cnot(0, 1); // prepare a Bell state

print!("{}", c.diagram());
// q0: -H-*-
// q1: ---X-

let counts = c.sample(1000, 42); // seeded; only 00 and 11 ever appear
```

The same circuits can be written as a small text program or as OpenQASM and
run from the command line:

```bash
cargo run -- programs/ghz.qsim    # text program
cargo run -- programs/bell.qasm   # OpenQASM 2.0
cargo run -- --help               # program format and options
```

Worked examples live in [`examples/`](https://github.com/samhithak13/qsimulator/tree/main/examples) and each verify their own
output:

```bash
cargo run --example ghz            # N-qubit GHZ state
cargo run --example grover         # Grover search (found with certainty)
cargo run --example qft            # quantum Fourier transform
cargo run --example teleportation  # measurement + classical feed-forward
cargo run --example noise          # channels against their analytic rates
cargo run -- programs/teleport.qsim # the same protocol as a program file
```

## What's here

- **State-vector simulation** up to ~25 qubits on a laptop; amplitudes take
  `16·2^n` bytes, which is the practical limit.
- **Gates** — X, Y, Z, H, S, T and their daggers, the phase gate P(λ), the
  general single-qubit U2/U3, Rx/Ry/Rz, SWAP, CNOT, CY, CZ, CH, controlled
  rotations, controlled-U3, Toffoli, Fredkin (CSWAP), and arbitrary
  (multi-)controlled unitaries.
- **Measurement** — single-qubit and full-register collapse plus Born-rule
  sampling, driven by a seedable, dependency-free RNG so runs reproduce.
- **Three front ends** — a Rust builder API, a line-based text format, and
  OpenQASM 2.0 import and export.
- **ASCII circuit diagrams** via `Circuit::diagram()` / `Display`.

The core depends only on `num-complex`.

## Testing

Every quantity has an independent oracle rather than a hand-copied expected
value:

- Gate matrices are checked against their algebraic identities — `Rx(π)`
  equals X up to global phase, `S·S† = I`, `p(π/2) = S`, and so on.
- The OpenQASM importer is validated by running an imported circuit and
  comparing its amplitudes to the equivalent builder circuit; export is
  checked by round-tripping back through the importer.
- Measurement is checked statistically: a Bell state only ever collapses to
  00 or 11, measuring one half of a pair fixes the other, and the same seed
  reproduces the same histogram.

Beyond the in-tree tests, [`crossval/`](https://github.com/samhithak13/qsimulator/tree/main/crossval) checks qsimulator against
Qiskit: it runs hundreds of random circuits over the shared gate set through
both engines (via OpenQASM) and confirms the state vectors match to
floating-point precision. See [crossval/README.md](https://github.com/samhithak13/qsimulator/blob/main/crossval/README.md).

The parsers are fuzzed with [`cargo fuzz`](https://github.com/samhithak13/qsimulator/tree/main/fuzz) (targets under `fuzz/`), with
a stable-toolchain random-input smoke test in `tests/robustness.rs` for
continuous coverage.

`cargo test` runs 120+ unit and integration tests (including a CLI driven
end-to-end). CI runs the full suite on Linux and macOS, enforces
`cargo fmt --check`, `cargo clippy -D warnings`, and a warning-clean
`cargo doc`, executes the examples, runs the Qiskit cross-validation, fuzzes
the parsers, audits dependencies with `cargo audit`, and gates line coverage
(~97%) at 90% via `cargo llvm-cov`.

## Benchmarks

`cargo bench` measures gate-application throughput (the simulator's inner
loop). The harness and a note on the kernel's structure — a bounds-check-free
walk over the target bit's `|0>`/`|1>` halves — are in
[benches/README.md](https://github.com/samhithak13/qsimulator/blob/main/benches/README.md). Numbers are machine-dependent, so
regenerate them on the target machine rather than trusting a checked-in figure.

An optional `parallel` feature runs the single-qubit kernel across threads
(via rayon); it is off by default so the core has no dependency beyond
`num-complex`. Enable it with `cargo build --features parallel`.

## Project layout

```
qsimulator/
├── src/
│   ├── lib.rs       # crate root, re-exports
│   ├── state.rs     # state vector and gate application
│   ├── gates.rs     # gate matrices
│   ├── circuit.rs   # circuit builder, execution, diagrams, QASM export
│   ├── program.rs   # text program parser
│   ├── qasm.rs      # OpenQASM 2.0 importer
│   ├── rng.rs       # seedable xorshift64 RNG
│   └── main.rs      # CLI
├── programs/        # example programs (.qsim and .qasm)
├── tests/           # integration tests
├── docs/design.md   # design notes and status
└── .github/         # CI
```

## Text program format

One instruction per line; `#` starts a comment. The first line sets the
register size.

```text
qubits 3
h 0
cnot 0 1
cnot 1 2        # GHZ state
sample 1000 42
```

Instructions: `qubits N`; `id|h|x|y|z|s|t|sdg|tdg|sx|sxdg Q`; `rx|ry|rz|p THETA Q`
(THETA is an arithmetic expression over numbers and `pi`, like `pi/2`, `-pi/4`,
`2pi`, or `(pi/4+0.1)*2` — one unspaced token, since instructions split on
whitespace);
`u2 PHI LAMBDA Q`; `u3 THETA PHI LAMBDA Q`; `cnot|cy|cz|ch C T`;
`crz|cp THETA C T`; `cu3 THETA PHI LAMBDA C T`; `swap A B`; `cswap C A B`;
`toffoli C1 C2 T`; `measure Q`; `reset Q`; `if VALUE INSTRUCTION`;
`depolarizing|bit_flip|phase_flip|amplitude_damping|phase_damping P Q`; the open-ended `mcx C... T` and
`mcu3 THETA PHI LAMBDA C... T`; and an optional `sample SHOTS SEED`.

Sampling flags work on any program, including OpenQASM files, which have no
way to carry a `sample` directive:

```bash
qsimulator --shots 1000 --seed 42 programs/bell.qasm
```

`--shots` overrides the program's own `sample` directive; `--seed` on its own
re-seeds it. Sampling is a pure function of the seed, so a run repeats exactly.

## OpenQASM 2.0 and 3

A file ending in `.qasm`, or one starting with an `OPENQASM` header, is parsed
as OpenQASM — the version in the header picks 2 or 3. The supported subset is
enough for hand-written textbook circuits and for what Qiskit exports:

```text
OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
h q[0];
cx q[0],q[1];
```

Supported gates: `id x y z h s t sdg tdg sx sxdg`, `rx/ry/rz(theta)`, `u1(lambda)`
(alias `p`), `u2(phi,lambda)`, `u3(theta,phi,lambda)`, `cx`, `cy`, `cz`, `ch`,
`crz(theta)`, `cu1(lambda)` (alias `cp`), `cu3(theta,phi,lambda)`, `swap`,
`ccx`, `cswap`, and the primitives `U`, `CX`; plus `gate` declarations, which
are expanded at each call site (so files Qiskit exports for composite gates
import), and the rest of `qelib1.inc` (`u`, `u0`, `crx`, `cry`, `csx`, `cu`,
`rxx`, `rzz`, `rccx`, `rc3x`, `c3x`, `c3sqrtx`, `c4x`) implemented natively,
since `include` is ignored rather than honoured.
Multiple `qreg`s map into one flat qubit space in declaration order;
`creg` and `barrier` are accepted and ignored; anything outside
the subset is reported rather than silently mis-simulated.

`Circuit::to_qasm()` (and `qsimulator --emit-qasm <FILE>`) writes a circuit
back out as OpenQASM 2.0, and every gate in the builder round-trips exactly
through the importer. Gates with no OpenQASM 2 equivalent are decomposed into
ones that have one: an arbitrary controlled-U into a control phase plus `cu3`,
and a multi-controlled gate of any width into Toffolis and single-qubit
rotations.

## Status

The core is complete and covered by tests: simulation, the gate set above,
measurement and sampling, the three front ends, ASCII diagrams,
cross-validation against Qiskit, a throughput benchmark harness, and an
optional parallel gate kernel.

Planned:

- `cu3` import and export, to round-trip arbitrary controlled-U gates.

Design notes and a detailed status table live in [docs/design.md](https://github.com/samhithak13/qsimulator/blob/main/docs/design.md).

## Contributing

See [CONTRIBUTING.md](https://github.com/samhithak13/qsimulator/blob/main/CONTRIBUTING.md). In short: keep new gates and primitives
covered by tests, and run `cargo fmt` and `cargo clippy` before opening a pull
request.

## License

[MIT](https://github.com/samhithak13/qsimulator/blob/main/LICENSE).
