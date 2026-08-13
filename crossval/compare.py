#!/usr/bin/env python3
"""Cross-validate qsimulator against Qiskit over the OpenQASM 2.0 bridge.

Two phases, each run for `--trials` random circuits:

  * **gates** — generate a random OpenQASM program from the gate set both
    engines implement, and run it through qsimulator (via its `--statevector`
    CLI, which prints the final amplitudes as JSON) and Qiskit's reference
    `Statevector`. This checks that the two engines agree gate for gate.

  * **measure** — generate a random program containing a *mid-circuit*
    measurement, and compare qsimulator's sampled distribution against Qiskit
    Aer's over the exported program. A circuit whose state collapses partway
    has no single state vector, so this is the only phase that can check it;
    the comparison is statistical, over total variation distance.

  * **export** — generate a random program in qsimulator's native text format,
    including multi-controlled gates that OpenQASM 2 has no way to write
    directly. qsimulator runs it, while Qiskit runs the program as qsimulator
    exports it (`--emit-qasm`), i.e. after decomposition. This checks that the
    decomposition an export goes through means the same thing to another tool.

State vectors are compared up to global phase (fidelity), which is the
physically meaningful notion of state equality and is robust to per-gate
global-phase conventions. A single mismatch prints the offending program and
exits non-zero.

Usage:
    python crossval/compare.py [--trials N] [--seed S] [--tol T] [--binary PATH]
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import subprocess
import sys

import numpy as np

# Gates shared by qsimulator and Qiskit's qelib1, grouped by (qubits, angles).
ONE_QUBIT = ["id", "h", "x", "y", "z", "s", "t", "sdg", "tdg", "sx", "sxdg"]
ONE_QUBIT_1ANGLE = ["rx", "ry", "rz", "u1"]
TWO_QUBIT = ["cx", "cy", "cz", "ch", "swap", "csx"]
TWO_QUBIT_1ANGLE = ["crz", "cu1", "crx", "cry", "rxx", "rzz"]
THREE_QUBIT = ["ccx", "cswap", "rccx"]
FOUR_QUBIT = ["c3x", "c3sqrtx", "rc3x"]


def angle(rng: random.Random) -> str:
    # repr() round-trips exactly through both parsers.
    return repr(rng.uniform(-2.0 * np.pi, 2.0 * np.pi))


def random_qasm(n_qubits: int, n_gates: int, rng: random.Random) -> str:
    lines = ["OPENQASM 2.0;", 'include "qelib1.inc";', f"qreg q[{n_qubits}];"]
    for _ in range(n_gates):
        kind = rng.random()
        if kind < 0.35:
            g = rng.choice(ONE_QUBIT)
            q = rng.randrange(n_qubits)
            lines.append(f"{g} q[{q}];")
        elif kind < 0.55:
            g = rng.choice(ONE_QUBIT_1ANGLE)
            q = rng.randrange(n_qubits)
            lines.append(f"{g}({angle(rng)}) q[{q}];")
        elif kind < 0.70:
            q = rng.randrange(n_qubits)
            lines.append(f"u2({angle(rng)},{angle(rng)}) q[{q}];")
        elif kind < 0.82:
            q = rng.randrange(n_qubits)
            lines.append(f"u3({angle(rng)},{angle(rng)},{angle(rng)}) q[{q}];")
        elif n_qubits >= 2 and kind < 0.94:
            a, b = rng.sample(range(n_qubits), 2)
            pick = rng.random()
            if pick < 0.45:
                g = rng.choice(TWO_QUBIT)
                lines.append(f"{g} q[{a}],q[{b}];")
            elif pick < 0.8:
                g = rng.choice(TWO_QUBIT_1ANGLE)
                lines.append(f"{g}({angle(rng)}) q[{a}],q[{b}];")
            elif rng.random() < 0.5:
                lines.append(
                    f"cu3({angle(rng)},{angle(rng)},{angle(rng)}) q[{a}],q[{b}];"
                )
            else:
                lines.append(
                    f"cu({angle(rng)},{angle(rng)},{angle(rng)},{angle(rng)}) "
                    f"q[{a}],q[{b}];"
                )
        elif n_qubits >= 4 and kind < 0.97:
            picks = rng.sample(range(n_qubits), 4)
            g = rng.choice(FOUR_QUBIT)
            lines.append(f"{g} " + ",".join(f"q[{i}]" for i in picks) + ";")
        elif n_qubits >= 5 and kind < 0.98:
            picks = rng.sample(range(n_qubits), 5)
            lines.append("c4x " + ",".join(f"q[{i}]" for i in picks) + ";")
        elif n_qubits >= 3:
            a, b, c = rng.sample(range(n_qubits), 3)
            g = rng.choice(THREE_QUBIT)
            lines.append(f"{g} q[{a}],q[{b}],q[{c}];")
        else:
            q = rng.randrange(n_qubits)
            lines.append(f"x q[{q}];")
    return "\n".join(lines) + "\n"


def random_program(n_qubits: int, n_gates: int, rng: random.Random) -> str:
    """A random circuit in qsimulator's native text format.

    Weighted towards the multi-controlled instructions, since these are the
    ones whose OpenQASM form is a decomposition rather than a single gate; the
    rest are there to hand them a thoroughly entangled input state.
    """
    lines = [f"qubits {n_qubits}"]
    for _ in range(n_gates):
        kind = rng.random()
        if kind < 0.3:
            g = rng.choice(ONE_QUBIT)
            lines.append(f"{g} {rng.randrange(n_qubits)}")
        elif kind < 0.45:
            lines.append(
                f"u3 {angle(rng)} {angle(rng)} {angle(rng)} {rng.randrange(n_qubits)}"
            )
        elif n_qubits >= 2 and kind < 0.55:
            a, b = rng.sample(range(n_qubits), 2)
            lines.append(f"{rng.choice(['cnot', 'cy', 'cz', 'ch', 'swap'])} {a} {b}")
        elif n_qubits >= 2:
            # A multi-controlled gate over a random subset of the register.
            # One control shy of the full register exercises the borrowed-qubit
            # decomposition; the full register exercises the square-root one.
            width = rng.randint(2, n_qubits)
            qubits = rng.sample(range(n_qubits), width)
            operands = " ".join(str(q) for q in qubits)
            if rng.random() < 0.5:
                lines.append(f"mcx {operands}")
            else:
                # theta = 0 is a diagonal gate, which exports by a separate
                # (phase-only) path, so hit it deliberately every so often.
                theta = "0" if rng.random() < 0.3 else angle(rng)
                lines.append(f"mcu3 {theta} {angle(rng)} {angle(rng)} {operands}")
        else:
            lines.append(f"x {rng.randrange(n_qubits)}")
    return "\n".join(lines) + "\n"


def random_measured_program(n_qubits: int, n_gates: int, rng: random.Random) -> str:
    """A random circuit whose state collapses partway through.

    Guarantees at least one measurement with gates after it — the case that
    actually branches — and measures every qubit at the end so that Aer's
    classical register holds the same thing qsimulator samples.

    The measured qubit is deliberately put into superposition just before the
    collapse and taken out of the computational basis just after. Both are
    needed for the phase to have any teeth:

    - collapsing a qubit already in |0> or |1> does nothing, and
    - a collapse commutes with anything diagonal or permutation-like in that
      basis (X, Z, S, T, CNOT, SWAP, ...),

    so a circuit lacking either would give the same distribution whether or not
    the measurement happened — and would pass even with the collapse ignored
    entirely, which is the bug this phase exists to catch. Sandwiching the
    measurement between two basis-changing gates is what turns lost coherence
    into a difference in the histogram. The surrounding gates stay random, so
    the probe runs in varied entanglement contexts.
    """
    lines = [f"qubits {n_qubits}"]

    def gate() -> str:
        if n_qubits >= 2 and rng.random() < 0.4:
            a, b = rng.sample(range(n_qubits), 2)
            return f"{rng.choice(['cnot', 'cz', 'cy', 'ch', 'swap'])} {a} {b}"
        q = rng.randrange(n_qubits)
        if rng.random() < 0.5:
            return f"{rng.choice(ONE_QUBIT)} {q}"
        return f"u3 {angle(rng)} {angle(rng)} {angle(rng)} {q}"

    # Gates, a collapse, then more gates so the collapse actually matters.
    for _ in range(max(1, n_gates // 2)):
        lines.append(gate())
    measured = rng.randrange(n_qubits)
    # Superposition in, collapse, back out of the basis — see the note above.
    # `reset` collapses the same way and then forces |0>, so it belongs here
    # too: it is the other non-unitary operation this engine implements.
    lines.append(rng.choice([f"h {measured}", f"sx {measured}"]))
    lines.append(rng.choice([f"measure {measured}", f"reset {measured}"]))
    lines.append(rng.choice([f"h {measured}", f"sx {measured}"]))
    for _ in range(max(1, n_gates // 2)):
        lines.append(gate())
    # Read every qubit out, so both engines report the same register.
    lines += [f"measure {q}" for q in range(n_qubits)]
    return "\n".join(lines) + "\n"


def qsim_counts(program: str, binary: str, shots: int, seed: int) -> dict[str, float]:
    """qsimulator's sampled distribution, parsed from its histogram output."""
    out = qsim(["--shots", str(shots), "--seed", str(seed)], program, binary)
    counts = {m[1]: int(m[2]) for m in re.finditer(r"\|(\d+)>: (\d+) shots", out)}
    total = sum(counts.values())
    if total != shots:
        raise RuntimeError(f"parsed {total} of {shots} shots from:\n{out}")
    return {k: v / total for k, v in counts.items()}


def aer_counts(qasm: str, shots: int, seed: int) -> dict[str, float]:
    """Qiskit Aer's sampled distribution over the same program."""
    from qiskit import qasm2, transpile
    from qiskit_aer import AerSimulator

    circuit = qasm2.loads(qasm, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    sim = AerSimulator(seed_simulator=seed)
    result = sim.run(transpile(circuit, sim), shots=shots).result().get_counts()
    # Qiskit orders a count key c[n-1]...c[0], and our export writes qubit i to
    # c[i], so the two engines already agree on bit order.
    total = sum(result.values())
    return {k.replace(" ", ""): v / total for k, v in result.items()}


def total_variation(a: dict[str, float], b: dict[str, float]) -> float:
    """Half the L1 distance between two distributions: 0 identical, 1 disjoint."""
    return 0.5 * sum(abs(a.get(k, 0.0) - b.get(k, 0.0)) for k in set(a) | set(b))


def qsim(args: list[str], program: str, binary: str) -> str:
    result = subprocess.run(
        [binary, *args, "-"],
        input=program,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout


def qsim_statevector(qasm: str, binary: str) -> np.ndarray:
    data = json.loads(qsim(["--statevector"], qasm, binary))
    return np.array([complex(re, im) for re, im in data], dtype=complex)


def qiskit_statevector(qasm: str) -> np.ndarray:
    from qiskit import qasm2
    from qiskit.quantum_info import Statevector

    circuit = qasm2.loads(qasm, custom_instructions=qasm2.LEGACY_CUSTOM_INSTRUCTIONS)
    # Qiskit and qsimulator both index little-endian (qubit 0 is the least
    # significant bit), so no reordering is needed.
    return np.asarray(Statevector(circuit).data, dtype=complex)


def fidelity(a: np.ndarray, b: np.ndarray) -> float:
    """State fidelity |<a|b>|^2, which ignores an overall global phase."""
    return float(np.abs(np.vdot(a, b)) ** 2)


def build_binary() -> str:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    subprocess.run(["cargo", "build", "--quiet"], cwd=root, check=True)
    return os.path.join(root, "target", "debug", "qsimulator")


def compare(
    label: str,
    trial: int,
    source: str,
    a: np.ndarray,
    b: np.ndarray,
    tol: float,
) -> float | None:
    """Fidelity of the two state vectors, or None (after printing) on a
    mismatch."""
    if a.shape != b.shape:
        print(f"FAIL ({label} trial {trial}): shape {a.shape} vs {b.shape}\n{source}")
        return None
    f = fidelity(a, b)
    if 1.0 - f > tol:
        print(f"FAIL ({label} trial {trial}): fidelity {f:.3e} below 1 - {tol:g}")
        print(source)
        return None
    return f


def run_gate_phase(args, binary: str, rng: random.Random) -> float | None:
    """Both engines run the same OpenQASM program."""
    worst = 1.0
    for trial in range(args.trials):
        n_qubits = rng.randint(1, args.max_qubits)
        n_gates = rng.randint(2 * n_qubits, 6 * n_qubits)
        qasm = random_qasm(n_qubits, n_gates, rng)

        f = compare(
            "gates",
            trial,
            qasm,
            qsim_statevector(qasm, binary),
            qiskit_statevector(qasm),
            args.tol,
        )
        if f is None:
            return None
        worst = min(worst, f)
    return worst


def run_export_phase(args, binary: str, rng: random.Random) -> float | None:
    """qsimulator runs a native program; Qiskit runs what qsimulator exports,
    so any error in the export decomposition shows up as a mismatch."""
    worst = 1.0
    for trial in range(args.trials):
        n_qubits = rng.randint(2, args.max_qubits)
        n_gates = rng.randint(n_qubits, 3 * n_qubits)
        program = random_program(n_qubits, n_gates, rng)
        exported = qsim(["--emit-qasm"], program, binary)

        f = compare(
            "export",
            trial,
            f"{program}\n--- exported as ---\n{exported}",
            qsim_statevector(program, binary),
            qiskit_statevector(exported),
            args.tol,
        )
        if f is None:
            return None
        worst = min(worst, f)
    return worst


def run_measure_phase(args, binary: str, rng: random.Random) -> float | None:
    """Both engines sample a circuit that collapses partway through.

    Both sides are sampled, so the comparison carries shot noise on both. At the
    default 20000 shots the worst total variation measured over 150 passing
    trials was 0.017, while ignoring the collapse — the bug this phase exists to
    catch — produces around 0.5. The default tolerance sits about 3x above the
    observed noise and 30x below that signal.
    """
    worst = 0.0
    trials = max(1, args.trials // 10)
    for trial in range(trials):
        n_qubits = rng.randint(2, min(4, args.max_qubits))
        program = random_measured_program(n_qubits, rng.randint(2, 3 * n_qubits), rng)
        exported = qsim(["--emit-qasm"], program, binary)
        seed = rng.randrange(1 << 30)

        mine = qsim_counts(program, binary, args.shots, seed)
        theirs = aer_counts(exported, args.shots, seed)
        distance = total_variation(mine, theirs)
        worst = max(worst, distance)
        if distance > args.measure_tol:
            print(
                f"FAIL (measure trial {trial}): total variation {distance:.4f} "
                f"exceeds {args.measure_tol:g}"
            )
            print(f"{program}\n--- exported as ---\n{exported}")
            print(f"qsimulator {mine}\naer        {theirs}")
            return None
    return worst


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trials", type=int, default=500, help="trials per phase")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--tol", type=float, default=1e-9)
    parser.add_argument("--max-qubits", type=int, default=5)
    parser.add_argument("--binary", default=None, help="path to the qsimulator binary")
    parser.add_argument(
        "--shots", type=int, default=20000, help="shots per measure-phase trial"
    )
    parser.add_argument(
        "--measure-tol",
        type=float,
        default=0.05,
        help="total-variation tolerance for the sampled measure phase",
    )
    args = parser.parse_args()

    binary = args.binary or build_binary()
    rng = random.Random(args.seed)

    for label, phase in (("gates", run_gate_phase), ("export", run_export_phase)):
        worst = phase(args, binary, rng)
        if worst is None:
            return 1
        print(
            f"OK: {args.trials} {label} trials agree with Qiskit "
            f"(worst fidelity {worst:.15f}, tol {args.tol:g})"
        )

    worst = run_measure_phase(args, binary, rng)
    if worst is None:
        return 1
    print(
        f"OK: {max(1, args.trials // 10)} measure trials agree with Aer "
        f"(worst total variation {worst:.4f}, tol {args.measure_tol:g})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
