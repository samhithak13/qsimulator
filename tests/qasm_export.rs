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

#[test]
fn error_exporting_arbitrary_controlled_u() {
    let mut c = Circuit::new(2);
    c.cu(gates::h(), 0, 1); // controlled-Hadamard: not in the QASM subset
    assert_eq!(c.to_qasm().unwrap_err(), ExportError::ControlledU);
}

#[test]
fn error_exporting_multi_controlled_u() {
    let mut c = Circuit::new(3);
    c.mcu(gates::z(), &[0, 1], 2);
    assert_eq!(c.to_qasm().unwrap_err(), ExportError::MultiControlledU);
}

#[test]
fn error_exporting_three_control_mcx() {
    let mut c = Circuit::new(4);
    c.mcx(&[0, 1, 2], 3); // C3X has no direct OpenQASM 2 gate
    assert_eq!(
        c.to_qasm().unwrap_err(),
        ExportError::MultiControlledX { controls: 3 }
    );
}
