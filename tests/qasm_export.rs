//! Integration tests for OpenQASM 2.0 export (`Circuit::to_qasm`).
//!
//! The core property: exporting a circuit and re-importing it yields an
//! equivalent circuit (same final-state probabilities).

use approx::assert_relative_eq;
use qsimulator::gates;
use qsimulator::{qasm, Circuit, ExportError, State};

fn assert_same_probs(a: &State, b: &State) {
    assert_eq!(a.n_qubits(), b.n_qubits());
    for i in 0..a.amplitudes().len() {
        assert_relative_eq!(a.probability(i), b.probability(i), epsilon = 1e-12);
    }
}

/// A circuit exercising every exportable gate must survive a round trip
/// through QASM export + import unchanged.
#[test]
fn round_trips_a_mixed_circuit() {
    let mut c = Circuit::new(3);
    c.h(0)
        .x(1)
        .y(2)
        .z(0)
        .s(1)
        .t(2)
        .sdg(0)
        .tdg(1)
        .p(0.9, 2)
        .u2(0.3, 1.2, 0)
        .u3(0.5, -0.6, 0.7, 1)
        .rx(0.7, 0)
        .ry(1.1, 1)
        .rz(-0.4, 2)
        .cnot(0, 1)
        .cz(1, 2)
        .crz(0.4, 0, 1)
        .cp(1.3, 1, 2)
        .swap(0, 2)
        .toffoli(0, 1, 2);

    let qasm_src = c.to_qasm().expect("should export");
    let reimported = qasm::parse(&qasm_src).expect("should re-import");
    assert_same_probs(&reimported.run(), &c.run());
}

#[test]
fn export_has_header_and_qreg() {
    let mut c = Circuit::new(4);
    c.h(0);
    let q = c.to_qasm().unwrap();
    assert!(q.starts_with("OPENQASM 2.0;"), "{q}");
    assert!(q.contains("qreg q[4];"), "{q}");
    assert!(q.contains("h q[0];"), "{q}");
}

/// Rotation angles must survive at full precision.
#[test]
fn rotation_angle_round_trips_exactly() {
    let mut c = Circuit::new(1);
    c.rx(0.7, 0);
    let q = c.to_qasm().unwrap();
    assert!(q.contains("rx(0.7) q[0];"), "{q}");

    let reimported = qasm::parse(&q).unwrap();
    assert_same_probs(&reimported.run(), &c.run());
}

/// A QASM string exported and re-imported twice stays stable.
#[test]
fn export_is_idempotent_through_parse() {
    let src = "qreg q[2];\nh q[0];\ncx q[0],q[1];\n";
    let once = qasm::parse(src).unwrap().to_qasm().unwrap();
    let twice = qasm::parse(&once).unwrap().to_qasm().unwrap();
    assert_eq!(once, twice);
}

/// u3 exports with all three angles and round-trips exactly.
#[test]
fn u3_exports_with_three_angles() {
    let mut c = Circuit::new(1);
    c.u3(0.5, -0.6, 0.7, 0);
    let q = c.to_qasm().unwrap();
    assert!(q.contains("u3(0.5,-0.6,0.7) q[0];"), "{q}");

    let reimported = qasm::parse(&q).unwrap();
    assert_same_probs(&reimported.run(), &c.run());
}

/// Controlled rotations export as `crz` / `cu1` and round-trip.
#[test]
fn controlled_rotations_export() {
    let mut c = Circuit::new(2);
    c.crz(0.4, 0, 1).cp(1.3, 0, 1);
    let q = c.to_qasm().unwrap();
    assert!(q.contains("crz(0.4) q[0],q[1];"), "{q}");
    assert!(q.contains("cu1(1.3) q[0],q[1];"), "{q}");

    let reimported = qasm::parse(&q).unwrap();
    assert_same_probs(&reimported.run(), &c.run());
}

/// The phase gate exports as OpenQASM `u1` and round-trips.
#[test]
fn phase_gate_exports_as_u1() {
    let mut c = Circuit::new(1);
    c.p(0.9, 0);
    let q = c.to_qasm().unwrap();
    assert!(q.contains("u1(0.9) q[0];"), "{q}");

    let reimported = qasm::parse(&q).unwrap();
    assert_same_probs(&reimported.run(), &c.run());
}

/// An arbitrary controlled-U is decomposed into a control phase plus cu3, so
/// it exports and round-trips (rather than erroring as it used to).
#[test]
fn arbitrary_controlled_u_round_trips() {
    // Controlled-Hadamard, plus superpositions on both qubits so the control
    // and target phases are actually exercised.
    let mut c = Circuit::new(2);
    c.h(0).h(1).cu(gates::h(), 0, 1);

    let qasm = c.to_qasm().expect("cu should export via decomposition");
    assert!(qasm.contains("cu3("), "{qasm}");

    let reimported = qasm::parse(&qasm).expect("should re-import");
    assert_same_probs(&reimported.run(), &c.run());
}

/// cy and ch export as their named OpenQASM gates and round-trip; cswap
/// round-trips through its CNOT/Toffoli decomposition.
#[test]
fn cy_ch_cswap_export() {
    let mut c = Circuit::new(3);
    c.h(0).h(1).h(2).cy(0, 1).ch(1, 2).cswap(0, 1, 2);

    let qasm = c.to_qasm().unwrap();
    assert!(qasm.contains("cy q[0],q[1];"), "{qasm}");
    assert!(qasm.contains("ch q[1],q[2];"), "{qasm}");

    let reimported = qasm::parse(&qasm).expect("should re-import");
    assert_same_probs(&reimported.run(), &c.run());
}

/// controlled-U3 exports as `cu3` and round-trips.
#[test]
fn cu3_exports_and_round_trips() {
    let mut c = Circuit::new(2);
    c.h(0).h(1).cu3(0.5, -0.6, 0.7, 0, 1);

    let qasm = c.to_qasm().unwrap();
    assert!(qasm.contains("cu3(0.5,-0.6,0.7) q[0],q[1];"), "{qasm}");

    let reimported = qasm::parse(&qasm).unwrap();
    assert_same_probs(&reimported.run(), &c.run());
}

/// Spread the register over every basis state with distinct, unequal
/// amplitudes, so a decomposition that is only *nearly* right shows up.
fn scramble(c: &mut Circuit, n_qubits: usize) {
    for q in 0..n_qubits {
        let k = q as f64;
        c.u3(0.7 + 0.3 * k, 0.4 - 0.2 * k, 1.1 + 0.5 * k, q);
    }
}

/// Assert two states are equal amplitude by amplitude — a stricter check than
/// [`assert_same_probs`], since a multi-controlled decomposition has to
/// reproduce the relative phases it conditions on, not just the populations.
fn assert_same_amplitudes(a: &State, b: &State) {
    assert_eq!(a.n_qubits(), b.n_qubits());
    for (i, (x, y)) in a.amplitudes().iter().zip(b.amplitudes()).enumerate() {
        assert!(
            (x - y).norm() < 1e-12,
            "amplitude {i} differs: {x} vs {y}\n{a:?}\n{b:?}"
        );
    }
}

/// Export a circuit, re-import it, and assert the two agree exactly.
fn assert_export_round_trips(c: &Circuit) -> String {
    let qasm = c.to_qasm().expect("should export");
    let reimported = qasm::parse(&qasm).expect("should re-import");
    assert_same_amplitudes(&reimported.run(), &c.run());
    qasm
}

/// A multi-controlled X round-trips for every control count, both when the
/// register has a qubit to borrow (the Toffoli-ladder decomposition) and when
/// the controls plus the target are the whole register (the square-root
/// recursion).
#[test]
fn multi_controlled_x_round_trips_at_every_width() {
    for controls in 0..=5usize {
        for spare in 0..=1usize {
            let n = controls + 1 + spare;
            let mut c = Circuit::new(n);
            scramble(&mut c, n);
            c.mcx(&(0..controls).collect::<Vec<_>>(), controls);
            assert_export_round_trips(&c);
        }
    }
}

/// With a qubit to borrow, a wide multi-controlled X decomposes into Toffolis
/// alone — no square roots, so no `cu3` in the output.
#[test]
fn multi_controlled_x_with_a_spare_qubit_is_all_toffolis() {
    let mut c = Circuit::new(6);
    c.mcx(&[0, 1, 2, 3], 4); // q5 is free to borrow
    let qasm = assert_export_round_trips(&c);
    assert!(!qasm.contains("cu3"), "{qasm}");
    assert!(qasm.contains("ccx"), "{qasm}");
}

/// The borrowed qubit must be returned in the state it was found in, whatever
/// that state is — including entangled with the rest of the register.
#[test]
fn borrowed_qubit_is_restored() {
    let mut c = Circuit::new(6);
    scramble(&mut c, 6);
    c.cnot(0, 5).h(5).cnot(5, 1); // entangle the qubit that will be borrowed
    c.mcx(&[0, 1, 2, 3], 4);
    assert_export_round_trips(&c);
}

/// An arbitrary multi-controlled unitary round-trips, with and without a spare
/// qubit.
#[test]
fn multi_controlled_u_round_trips() {
    for n in [4usize, 5] {
        let mut c = Circuit::new(n);
        scramble(&mut c, n);
        c.mcu(gates::h(), &[0, 1, 2], 3);
        c.mcu(gates::u3(0.6, -1.2, 0.3), &[1, 2, 3], 0);
        assert_export_round_trips(&c);
    }
}

/// A diagonal multi-controlled gate (here a CCZ and a three-control T) is pure
/// phase, so as long as its inner Toffolis have a qubit to borrow it exports
/// without any Y rotation.
#[test]
fn multi_controlled_diagonal_exports_as_phases() {
    let mut c = Circuit::new(5);
    scramble(&mut c, 5);
    c.mcu(gates::z(), &[0, 1], 2).mcu(gates::t(), &[0, 1, 2], 3);
    let qasm = assert_export_round_trips(&c);
    assert!(!qasm.contains("ry("), "{qasm}");
    assert!(qasm.contains("cu1("), "{qasm}");
}

/// A single-control diagonal gate takes the phase path too, which is shorter
/// than the general `cu3` form: a phase on the control (only when the gate has
/// one, as `rz` does) plus a controlled phase.
#[test]
fn single_control_diagonal_exports_as_phases() {
    let mut c = Circuit::new(2);
    c.h(0).h(1).mcu(gates::rz(0.7), &[0], 1);
    let qasm = assert_export_round_trips(&c);
    assert!(qasm.contains("u1(-0.35) q[0];"), "{qasm}");
    assert!(qasm.contains("cu1(0.7) q[0],q[1];"), "{qasm}");
    assert!(!qasm.contains("cu3"), "{qasm}");
}

/// Degenerate control counts fall back to the plain gates.
#[test]
fn multi_controlled_with_zero_or_one_control() {
    let mut c = Circuit::new(2);
    c.h(0).mcx(&[], 1).mcx(&[0], 1).mcu(gates::h(), &[0], 1);
    let qasm = assert_export_round_trips(&c);
    assert!(qasm.contains("x q[1];"), "{qasm}");
    assert!(qasm.contains("cx q[0],q[1];"), "{qasm}");
}

/// An `mcu` with no controls is unconditional, so its global phase is dropped:
/// the re-imported circuit matches up to that phase, i.e. in probabilities.
#[test]
fn uncontrolled_mcu_round_trips_up_to_global_phase() {
    let mut c = Circuit::new(2);
    c.h(0).mcu(gates::t(), &[], 0);

    let qasm = c.to_qasm().expect("should export");
    let reimported = qasm::parse(&qasm).expect("should re-import");
    assert_same_probs(&reimported.run(), &c.run());
}

/// The only remaining export failure: an `Op` label with no OpenQASM name.
/// Unreachable through the builder API, so this checks the message instead.
#[test]
fn single_gate_export_error_reads_well() {
    let e = ExportError::SingleGate { label: "W" };
    assert_eq!(
        e.to_string(),
        "cannot export single-qubit gate `W` to OpenQASM 2"
    );
}
