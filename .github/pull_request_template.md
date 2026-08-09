## What this changes

<!-- One or two sentences. If it fixes an issue, link it. -->

## Why

<!-- What was wrong or missing. For a bug fix, what the wrong behaviour was. -->

## Checks

<!-- CI runs all of these; running them first is faster than a red build. -->

- [ ] `cargo fmt --all -- --check`
- [ ] `cargo clippy --all-targets -- -D warnings` (and with `--features parallel`)
- [ ] `cargo test`
- [ ] `RUSTDOCFLAGS="-D warnings" cargo doc --no-deps`
- [ ] Tests land in the same commit as the change; a bug fix has a regression test
- [ ] `CHANGELOG.md` updated under `[Unreleased]`

<!-- If this touches gates or OpenQASM, the Qiskit cross-validation is the
     check that matters -- see crossval/README.md. Say whether you ran it and
     with how many trials. -->
