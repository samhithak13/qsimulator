# qsimulator — Creation Methodology & Design Rationale

This document is the "why" behind `qsimulator`. Where `docs/design.md` is a
living status/handoff log and the code comments explain *what* each piece
does, this file records the **reasoning**: which decisions were made, what
alternatives were rejected, and how a future contributor (human or agent)
should proceed and why. It is meant to be read top-to-bottom once, then used
as a reference.

---

## 1. What we are building, and the one idea it rests on

`qsimulator` is an **ideal (noiseless) state-vector quantum circuit
simulator** written in Rust. It models an `n`-qubit register as the full
vector of `2^n` complex amplitudes and evolves it by applying unitary gates,
then measures using the Born rule.

The entire project rests on a single modeling choice: **represent the state
explicitly and exactly, and pay the `2^n` memory cost for it.** Everything
else — the gate kernels, the builder API, the CLI — is a consequence of that
choice. We accepted it because the goal is *clarity and correctness first*, as
a foundation for learning and experimentation, not maximal qubit count. State
vector is the model you can read, test against textbook truth tables, and
reason about line by line. The stated ceiling (~20–25 qubits on a laptop) is a
feature of the honesty of the model, not a bug to engineer around yet.

### Why Rust

- **Correctness with zero GC pauses** and predictable memory: a `2^n` vector
  is the dominant cost, and Rust lets us keep it as one flat `Vec<Complex64>`
  with in-place kernels and no hidden allocation per gate.
- **A test culture that's cheap to honor**: `cargo test`, `cargo clippy`, and
  `cargo fmt` are one command each and gate every change (see §6).
- **`num-complex`** gives us a well-tested `Complex64` so we never hand-roll
  complex arithmetic — the one external dependency, and a deliberate one.

### Non-goals (and why they're deferred, not rejected)

Noise/density-matrix simulation, GPU/distributed execution, and a hardware
backend are all **out of scope until the core is stable**. Each of them
multiplies the state representation or the execution model; adding any before
the noiseless core is rock-solid would mean debugging two hard things at once.
They are revisitable — the roadmap leaves the door open.

---

## 2. The frozen conventions (decide once, never drift)

A simulator is only trustworthy if its conventions are *uniform*. Three were
fixed at the very first commit and are treated as immutable, because changing
any of them silently would corrupt every downstream test and golden value.

1. **Little-endian qubit order.** Bit `q` of a basis-state index is qubit `q`.
   So the 2-qubit index `0b10` means qubit 1 is |1⟩ and qubit 0 is |0⟩.
   *Why:* it makes a single-qubit gate a stride-`2^q` butterfly over amplitude
   pairs — the index arithmetic falls out directly from the bit position, with
   no remapping table. Big-endian would have forced an `n - 1 - q` flip into
   every kernel and every test.

2. **Row-major gate storage**, `[[m00, m01], [m10, m11]]` of `Complex64`.
   *Why:* it reads exactly like the matrix on the page, so translating a
   textbook unitary into code is transcription, not transposition — the single
   biggest source of silent gate bugs is eliminated by matching the notation.

3. **Rotations are `exp(-i·θ/2·P)`** for Pauli `P`, so `R_axis(π)` equals its
   Pauli *up to the global phase `-i`*. *Why:* this is the standard physics
   convention; picking it (and asserting the `-i` relationship exactly in
   tests) means our rotations compose correctly with any external reference
   and the global-phase bookkeeping is explicit rather than accidental.

A fourth invariant is enforced structurally, not just documented: **controls
must differ from their target** (asserted), and every gate kernel processes
each amplitude pair **exactly once**, always from the target-bit-0 side. This
is what guarantees a gate is applied once and only once — see §4.

---

## 3. Module architecture — and why the boundaries fall where they do

```
src/
├── lib.rs      # crate root, module wiring, re-exports (Circuit, State, Rng)
├── state.rs    # the 2^n amplitude vector + all apply_*/measure_*/swap kernels
├── gates.rs    # 2x2 unitary matrices (pure data), + rotation unit tests
├── circuit.rs  # Op enum + Circuit builder + run() + sample()
├── program.rs  # text program parser -> Program { circuit, shots, seed }
├── rng.rs      # seedable xorshift64 RNG (SplitMix64 seeding)
└── main.rs     # CLI: program file / stdin / built-in demo
```

The layering is deliberate and one-directional:

- **`gates.rs` is pure data.** A gate is just a `[[Complex64; 2]; 2]`; the
  functions there return matrices and know nothing about state, circuits, or
  qubits. *Why:* gates are the most-copied-from-a-textbook part of the code,
  so isolating them means each can be unit-tested as a matrix (unitarity,
  `R(π) = -i·P`) with zero machinery around it.

- **`state.rs` owns all the physics.** Every way the amplitude vector can
  change — single-qubit, controlled, multi-controlled, swap, and all three
  measurement variants — lives here as a method on `State`. *Why:* the
  amplitude vector is the one piece of mutable, correctness-critical state; if
  mutations were scattered across modules, the "touch each pair once"
  invariant would be impossible to audit. Keeping it in one file makes the
  invariant a local property you can verify by reading one screen.

- **`circuit.rs` is orchestration only.** The `Op` enum records *what* to do
  (`Single`, `Controlled`, `Swap`, `MultiControlled`); `run()` replays those
  ops against a fresh `State`. The builder never does amplitude math — it only
  pushes `Op`s. *Why:* this separation means adding a gate to the *builder*
  never risks the *kernels*, and a `Circuit` is a cheap, inspectable,
  re-runnable value (essential for `sample()`, which re-runs per shot).

- **`program.rs` and `main.rs`** are the outermost layer: text in, `Circuit`
  out, results printed. They depend on everything below and nothing depends on
  them. *Why:* the library is fully usable without the CLI; the CLI is a thin
  convenience that can grow (or be replaced) without touching the core.

The rule of thumb this encodes: **dependencies point inward toward the
amplitude vector, and the physics never leaks outward.**

---

## 4. The gate kernels — one pattern, generalized

Every gate mutation in `state.rs` is the same butterfly, specialized by a
predicate on the basis index `j`:

- **`apply_1q`** walks `j` in blocks of `2^(target+1)`, pairing `j` with
  `j + 2^target`, and applies the 2×2 to *every* such pair.
- **`apply_controlled_1q`** applies the same 2×2 only where the control bit is
  set (`(j & cmask) != 0`) and the target bit is 0.
- **`apply_multi_controlled_1q`** generalizes the control test to
  `(j & cmask) == cmask` — *all* controls set. Zero controls degenerate to an
  unconditional gate, one control reproduces the single-control case, and two
  controls with X is Toffoli.
- **`swap_qubits`** exchanges the amplitudes of the two indices that differ
  only in the two swapped bits, touching each pair once from the
  (bit-a-set, bit-b-clear) side.

**Why build it this way:** the multi-controlled kernel is the *general* form;
CNOT and Toffoli are not special cases in the code so much as thin builders on
top of it (`cnot` → controlled-X, `toffoli` → `MultiControlled` with two
controls). This is why `docs/design.md` can say the state layer "already
supports any control count" — the hard part was written once. When the
in-flight builder PRs add `mcx`/`mcu`, they add **no new state machinery**;
they're pure builder + `Op` pass-through. That was the payoff we were buying
by generalizing early.

Everything is **in place**: one `Vec<Complex64>` for the whole run, no
per-gate allocation. That's the §1 memory choice honored all the way down.

---

## 5. Measurement and randomness — reproducibility as a first principle

Measurement uses the Born rule (`p(i) = |amplitude(i)|²`) with post-measurement
collapse and renormalization. Three entry points cover the real use cases:
`prob_qubit_one` (read-only), `measure_qubit` (single-qubit sample + collapse),
and `measure_all` (full-register sample + collapse). `Circuit::sample` runs the
circuit once, then measures independent **clones** of the final state — so the
shots are genuinely independent and the expensive `run()` happens once.

The RNG choice is the load-bearing decision here. It is a **seedable,
dependency-free `xorshift64`**, seeded through the **SplitMix64** finalizer:

- **Dependency-free** because a simulator's randomness must be *understood and
  stable*, not inherited from a crate that might change its stream between
  versions.
- **Seedable and fully deterministic in the seed** because reproducibility is
  a correctness property, not a convenience: tests assert exact histograms,
  and a user re-running `sample(shots, seed)` must get identical results. This
  is exactly why the design doc marks the RNG recurrence and seeding as frozen
  — changing them breaks every golden value.
- **SplitMix64 seeding (and forcing the low bit set)** because raw xorshift is
  fatal at a zero state and weak on structured seeds; running the seed through
  a strong finalizer means even `seed = 0` yields a healthy, non-degenerate
  stream (`seed_zero_is_not_degenerate` guards this).

Two numerical guards in the kernels are deliberate, not incidental:
`measure_qubit` avoids dividing by a zero-probability branch, and `measure_all`
defaults its outcome to the last nonzero-amplitude state so floating-point
roundoff in the cumulative sum can never fall through to a wrong index.

---

## 6. Testing philosophy — every gate has a known truth, so assert it exactly

The strategy is simple and strictly enforced: **every gate and circuit
primitive has a known matrix or truth table, so tests assert exact
probabilities within a tight epsilon (`1e-12`).** There is no "looks about
right" tolerance for algebraic facts.

Two layers:

- **Unit tests** (in `gates.rs`, `rng.rs`, `program.rs`) check the smallest
  claim in isolation: a rotation is unitary, `R(π) = -i·P`, the RNG is
  reproducible, an angle string parses.
- **Integration tests** (`tests/`) check *behavior* of canonical circuits:
  Bell and GHZ probability splits, single-gate flips, SWAP, Toffoli truth
  tables, and program parsing/sampling end-to-end.

The **process rule** is the important part: *new functionality lands with its
test in the same commit*, and all four checks —
`cargo fmt --all -- --check`, `cargo clippy --all-targets -- -D warnings`,
`cargo build`, `cargo test` — must pass before committing. CI
(`.github/workflows/ci.yml`) runs exactly those four on every push and PR, so
"green locally" and "green in CI" are the same bar. This is why the test count
is a tracked number in the handoff log (currently **42**): it's a proxy for
how much of the surface is pinned down.

**Why clippy at `-D warnings`:** it turns lint into a hard gate, which caught a
real subtlety during the CLI work — a test used `1.5708` as a plain float and
clippy flagged it as an approximation of `π/2`. The fix (use `0.75`) kept the
test honest about what it was checking. A soft warning would have been ignored.

---

## 7. How the project was grown — incremental, each commit shippable

The git history is the methodology in miniature. Each commit is one coherent,
tested capability, in dependency order:

1. **Scaffold** — state vector, `apply_1q`, the Pauli/H/S/T matrices, the
   Bell-state demo, CI. The skeleton that everything hangs on.
2. **Measurement** — the RNG, Born-rule sampling, collapse. Added early
   because *observing* the state is how you test everything after it.
3. **Rotations** — `rx/ry/rz`, with the `R(π) = -i·P` invariant asserted.
4. **SWAP** — the first two-qubit permutation kernel.
5. **Toffoli** — introduced `apply_multi_controlled_1q`, the *general*
   controlled kernel, and made Toffoli a thin builder on it (§4).
6. **Handoff log** — `docs/design.md`, so any session (including a fresh
   agent picking up from the web app) can continue without re-deriving context.
7. **Text program format + CLI** (this branch) — the most recent step.

The through-line: **never add a feature whose foundation isn't already tested.**
Measurement before rotations (so rotations can be verified by sampling); the
general multi-controlled kernel before exposing `mcx`/`mcu` builders (so the
builders are trivial and safe). Each commit leaves `main` in a shippable,
green state — there is no "work in progress" commit that doesn't build.

---

## 8. The most recent decision: the text program format & CLI

This is the clearest worked example of the methodology applied to a *choice
under constraints*, so it's worth spelling out.

**The situation.** The `docs/design.md` handoff listed six next steps. But
inspection of the open PRs showed that steps #1 (more builders: `y/s/t/cz/cu`)
and #2 (`mcx`/`mcu`) were **already in flight in PRs #2 and #3**, and #3
(circuit `Display`) was included in PR #2.

**The decision.** Rather than duplicate in-flight work — which would create
merge conflicts and waste effort — this branch took the next *non-overlapping*
item: **#4, a richer CLI**, which also incidentally delivers #5's GHZ fixture.

**Why the scope is deliberately narrow.** The parser (`src/program.rs`)
supports **only the builders already on `main`** (`h, x, z, rx, ry, rz, cnot,
swap, toffoli`). It does *not* add `y/s/t/cz/cu/mcx/mcu`, precisely because
those builders are what PRs #2/#3 add. Using only what's already merged means
this branch **cannot conflict** with them, and extending the parser once they
land is a trivial follow-up (noted in the design doc). This is branch
coordination as a design constraint: *pick work that composes with, rather
than races against, parallel branches.*

**What was built, and the small choices inside it:**

- A **line-based format** — `qubits N`, one gate per line, `sample <shots>
  [seed]`, `#` comments, blank lines ignored. Line-based (not JSON/TOML)
  because a quantum circuit *is* a sequence of instructions; the format should
  read like the circuit, matching the same "code mirrors the notation"
  principle as row-major gates (§2).
- **`pi`-expression angles** (`pi`, `-pi`, `pi/2`, `2pi`, `2*pi`, `3pi/4`)
  because rotation angles are almost always fractions of π; forcing users to
  precompute `1.5708` would be hostile and error-prone.
- **Errors carry a 1-based line number** — a parser that says *where* is worth
  far more than one that only says *what*.
- **Validation mirrors the kernels' asserts** at parse time (control ≠ target,
  qubit in range) so a bad program fails with a friendly message instead of a
  panic deep in `state.rs`.
- **`main.rs` returns a real `ExitCode`** and supports file / stdin (`-`) /
  no-args-demo / `--help`, so the tool is scriptable and pipeable.

It shipped with `examples/ghz.qsim`, unit tests in `program.rs`, integration
tests in `tests/program.rs`, README and design-doc updates, and all four
checks green — 28 → 42 tests. Exactly the process rule from §6.

---

## 9. How to proceed from here — and why, in this order

The priority order below reflects the same principle throughout: **do the
smallest thing that unlocks the most, and never build on an untested base.**

1. **Merge/land the in-flight builder PRs (#2/#3) first.** They add
   `y/s/t/cz/cu/mcx/mcu` and a circuit `Display`. Everything else is smoother
   once the builder surface is complete. *Why first:* they're already written
   and reviewed; leaving them open forces every other branch to route around
   them (as this one did).

2. **Then extend `program.rs` to the new builders.** Once `y/s/t/cz/cu/mcx/mcu`
   exist on `Circuit`, adding them to the parser is a few match arms plus tests
   — trivial *because* the CLI was deliberately built on the stable subset.
   *Why:* it closes the gap between what the library can do and what a program
   file can express, at near-zero cost.

3. **Round out inspection.** With `Display` merged, print the circuit diagram
   in the CLI output and from the demo. *Why:* cheap, high-signal ergonomics;
   makes the tool self-documenting.

4. **Richer program semantics only if needed** — e.g. named registers,
   parameters, or barriers. *Why the caution:* every syntax addition is a
   permanent compatibility surface; add it when a real use case demands it, not
   speculatively.

5. **Performance (v0.3) last.** In-place kernels already exist, so the honest
   next step is *measurement*: add benchmarks first, then consider a sparse
   fast path only where the benchmark shows it pays. *Why last:* optimizing
   before measuring, on a correctness-first project, risks trading the thing we
   actually value (clarity, testability) for speed we can't prove we need.

**The rule for any new gate**, restated because it's the whole method: add the
matrix in `gates.rs`, a builder in `circuit.rs` (plus an `Op` variant and
`run` arm if it needs new state machinery), a parser arm in `program.rs`, and
*both* a unit test (matrix/truth-table) and an integration test (circuit
behavior) — in the same commit, with all four checks green.

---

## 10. Where to look

| You want to… | Read |
|---|---|
| Understand the model & conventions | `docs/design.md` §State-vector, this file §2 |
| See current status / what's done / test count | `docs/design.md` "Current status & handoff" |
| Understand *why* a decision was made | this file |
| Change a gate kernel | `src/state.rs` (mind the "touch once" invariant) |
| Add a gate matrix | `src/gates.rs` |
| Add a builder | `src/circuit.rs` |
| Extend the CLI language | `src/program.rs` + `src/main.rs` |
| Know what must never change | this file §2, `docs/design.md` "Frozen conventions" |

`design.md` is the *dashboard*; this file is the *rationale*. Keep them in
sync: when a decision here is superseded, record the new reasoning rather than
deleting the old — the trail of *why* is the point.
