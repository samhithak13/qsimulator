//! Integration tests for ASCII circuit-diagram rendering.

use qsimulator::Circuit;

#[test]
fn bell_diagram() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1);
    assert_eq!(c.diagram(), "q0: -H-*-\nq1: ---X-");
}

#[test]
fn toffoli_diagram() {
    let mut c = Circuit::new(3);
    c.toffoli(0, 1, 2);
    assert_eq!(c.diagram(), "q0: -*-\nq1: -*-\nq2: -X-");
}

#[test]
fn cnot_gap_draws_connector() {
    // A control on q0 and target on q2 must draw a `|` through q1.
    let mut c = Circuit::new(3);
    c.cnot(0, 2);
    assert_eq!(c.diagram(), "q0: -*-\nq1: -|-\nq2: -X-");
}

#[test]
fn swap_diagram() {
    let mut c = Circuit::new(3);
    c.swap(0, 2);
    assert_eq!(c.diagram(), "q0: -x-\nq1: -|-\nq2: -x-");
}

#[test]
fn multi_column_sequence() {
    // X on q0, then H on q1, then CNOT(0,1): three columns, left to right.
    let mut c = Circuit::new(2);
    c.x(0).h(1).cnot(0, 1);
    assert_eq!(c.diagram(), "q0: -X---*-\nq1: ---H-X-");
}

#[test]
fn display_matches_diagram() {
    let mut c = Circuit::new(2);
    c.h(0).cnot(0, 1);
    assert_eq!(format!("{c}"), c.diagram());
}

#[test]
fn wider_labels_stay_aligned() {
    // Mixing a 1-char and a 2-char label widens every cell to 2; rows must
    // remain equal length.
    let mut c = Circuit::new(2);
    c.h(0).rx(std::f64::consts::PI, 1);
    let diagram = c.diagram();
    let lines: Vec<&str> = diagram.lines().collect();
    assert_eq!(lines.len(), 2);
    assert_eq!(lines[0].len(), lines[1].len());
    assert!(lines[1].contains("RX"));
}

#[test]
fn empty_circuit_renders_bare_wires() {
    let c = Circuit::new(2);
    let d = c.diagram();
    assert!(d.contains("q0:"));
    assert!(d.contains("q1:"));
}
