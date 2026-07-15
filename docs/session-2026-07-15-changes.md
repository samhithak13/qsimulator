# Session handoff — builder expansion (2026-07-15)

Hand this file to the next Claude session to continue without re-deriving
context. It records **what changed, why, and what's left**. The living
design/status doc is `docs/design.md`; this file is the changelog for one
session's work.

## TL;DR

Resumed from the handoff log in `docs/design.md` and completed **steps 1 & 2**
of its "Suggested next steps" list: expanded the `Circuit` builder API with the
remaining single-qubit gates and generalized the multi-control machinery. No
changes to state-vector internals or the project's frozen conventions.

- Branch: `claude/resume-previous-session-8ga5ls`
- PR: **#5** (draft) → `main`, CI green ("Build, test, lint")
- Tests: **28 → 38** (added `tests/builders.rs`, 10 tests)

## What changed

All changes are in `src/circuit.rs` (new builder methods) plus a new test file.
The underlying `State` gate-application methods and the `gates::*` matrices
already existed — this session only added builder surface on top of them.

| New builder | Signature | What it does |
|---|---|---|
| `y` | `y(target)` | Pauli-Y single-qubit gate |
| `s` | `s(target)` | Phase gate S = diag(1, i) |
| `t` | `t(target)` | T gate = diag(1, e^{iπ/4}) |
| `cz` | `cz(control, target)` | Controlled-Z (phase-flip on \|11>) |
| `cu` | `cu(gate, control, target)` | Controlled arbitrary 2×2 unitary |
| `mcx` | `mcx(controls: &[usize], target)` | Multi-controlled X (any # of controls) |
| `mcu` | `mcu(gate, controls: &[usize], target)` | Multi-controlled arbitrary 2×2 unitary |

`toffoli(c1, c2, target)` is unchanged in behavior but is now conceptually just
the two-control special case of `mcx(&[c1, c2], target)`.

## Why these changes

- **Directly follows the plan.** `docs/design.md` listed these as steps 1 & 2,
  flagged "small, high-value" because the hard part (state-vector application)
  was already implemented and tested. This was the lowest-risk, highest-value
  next increment.
- **Closes obvious API gaps.** The gate *matrices* for `y`, `s`, `t` already
  existed in `src/gates.rs` but were not reachable from the `Circuit` builder,
  so a user could not actually place a Y/S/T gate in a circuit. These builders
  make the full documented gate set usable.
- **Generalizes instead of hard-coding.** `State::apply_multi_controlled_1q`
  already supported any number of controls, but `Circuit` only exposed the
  fixed 2-control `toffoli`. `mcx`/`mcu` expose that existing capability
  generally (n-controlled gates, e.g. a 3-controlled X), rather than adding a
  new one-off builder per control count.
- **`cz`/`cu` reuse existing controlled application.** `apply_controlled_1q`
  already handled an arbitrary 2×2 unitary, so controlled-Z and controlled-U
  are thin pass-throughs — no new state machinery, minimal surface area.

## How it was verified

`tests/builders.rs` (10 integration tests) asserts exact amplitudes/probabilities
against known truth tables, matching the repo's existing testing strategy:

- `y_maps_zero_to_i_one` — Y\|0> = i\|1>
- `s_phases_one_by_i` — S phases \|1> by i
- `t_squared_is_s` — T·T = S on \|1>
- `cz_phases_only_eleven` / `cz_is_symmetric` — CZ negates \|11> only; cz(a,b)==cz(b,a)
- `cu_with_x_is_cnot` — controlled-U with U=X reproduces CNOT
- `mcx_zero_controls_is_x`, `mcx_one_control_is_cnot`, `mcx_three_controls` —
  mcx degrades correctly to X / CNOT and gates only when *all* controls are set
- `mcu_with_z_phases_target` — mcu with Z phases \|111> by -1

All four project gates pass: `cargo fmt --all -- --check`,
`cargo clippy --all-targets -- -D warnings`, `cargo build`, `cargo test`
(38 tests). CI on PR #5 is green.

## What's left (for the next session)

Remaining items from `docs/design.md`, in priority order — none started:

3. **Circuit diagram printing** — a `Display`/ASCII rendering for `Circuit`
   (README roadmap item under v0.2).
4. **Richer CLI** (`src/main.rs`) — accept a gate list / simple program instead
   of the hard-coded Bell demo.
5. **GHZ + multi-qubit fixtures** — an explicit GHZ prep + sampling test.
6. **Performance (v0.3)** — benchmarks and a possible sparse fast path, only
   after the above.

## Conventions to keep (do not change silently)

Carried over from `docs/design.md` — reproducibility and tests depend on them:

- **Little-endian qubit order**: bit `q` of a basis index is qubit `q`.
- **Gate storage** is row-major `[[m00, m01], [m10, m11]]`, `Complex64`.
- **Rotations** are `exp(-i·θ/2·P)`; `R_axis(π)` equals its Pauli up to `-i`.
- **RNG is deterministic in its seed** (xorshift64 + SplitMix64 seeding).
- Controls must differ from the target (asserted in `State`).
- New functionality lands with its test in the same commit, and all four
  project gates must pass before committing.
