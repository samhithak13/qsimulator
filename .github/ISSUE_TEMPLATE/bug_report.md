---
name: Bug report
about: A wrong result, a crash, or something that does not behave as documented
labels: bug
---

## What happened

<!-- If a circuit produced the wrong numbers, that is the most valuable kind
     of report this project can get. Say what you expected and why. -->

## Reproducing it

<!-- The smallest circuit that shows the problem. A program in the native text
     format or OpenQASM is ideal, since it can be pasted straight into the CLI:

     qsimulator --statevector - <<'END'
     qubits 2
     h 0
     cnot 0 1
     END
-->

```text

```

## Expected vs actual

<!-- Amplitudes, probabilities, or a sampled histogram. If you have a
     reference (Qiskit, pen and paper, a textbook identity), say which --
     it usually pins down whether the bug is in the engine or a convention
     mismatch. Note that qsimulator indexes qubits little-endian, so basis
     state |10> means qubit 1 is set. -->

## Environment

- qsimulator version or commit:
- `rustc --version`:
- OS:
- Features enabled (e.g. `parallel`):
