//! Differential test of `lift_basis_normalized`. Same shape as
//! `lift_basis`, but the inputs are pre-normalised (B.P.z == 1 and
//! E.C == 1).
mod common;

use sqisign_ec::{lift_basis_normalized, JacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/lift_basis_normalized.json"
);

#[test]
fn lift_basis_normalized_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::lift_basis_normalized");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let b = common::ec_basis_from("b_in", &v.inputs);
        let e = common::ec_curve_from("e_in", &v.inputs);
        let exp_p = common::jac_point_from("p", &v.outputs);
        let exp_q = common::jac_point_from("q", &v.outputs);
        let exp_b = common::ec_basis_from("b_out", &v.outputs);
        let exp_e = common::ec_curve_from("e_out", &v.outputs);
        let exp_ret = common::u32_field("ret", &v.outputs);
        let mut p = JacPoint::zero();
        let mut q = JacPoint::zero();
        let got_ret = lift_basis_normalized(&mut p, &mut q, &b, &e);
        assert_eq!(p, exp_p, "vector {} p diverged", v.id);
        assert_eq!(q, exp_q, "vector {} q diverged", v.id);
        // The reference's C signature is non-const but the body never
        // writes back to the basis or curve; the recorded `b_out` and
        // `e_out` are therefore byte-identical to the inputs.
        assert_eq!(b, exp_b, "vector {} basis diverged", v.id);
        assert_eq!(e, exp_e, "vector {} curve diverged", v.id);
        assert_eq!(got_ret, exp_ret, "vector {} ret diverged", v.id);
    }
}
