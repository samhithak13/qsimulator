# Fuzzing

[`cargo fuzz`](https://rust-fuzz.github.io/book/cargo-fuzz.html) targets for the
two parsers. The contract is simple: for *any* byte input, `qasm::parse` and
`program::parse` must return a `Result` and never panic. The targets only
parse — they never run a circuit — so no state vector is allocated regardless
of what qubit count the input declares.

## Running

Requires a nightly toolchain and `cargo-fuzz`:

```bash
cargo install cargo-fuzz
cargo +nightly fuzz run qasm_parse
cargo +nightly fuzz run program_parse
```

A time-boxed run (used in CI):

```bash
cargo +nightly fuzz run qasm_parse -- -max_total_time=60
```

If a target finds a crash it writes the reproducing input under
`fuzz/artifacts/`. The two parser panics fixed in the changelog (malformed
bracket order, and a huge register aborting on allocation) are exactly the
class of bug these guard against; the stable-toolchain smoke test in
`tests/robustness.rs` provides continuous coverage without nightly.
