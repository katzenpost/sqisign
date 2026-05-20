//! Differential test of `lift_basis`. Inputs: a basis and a curve (both
//! mutated). Outputs: P, Q (Jacobian), the mutated basis and curve, and
//! the u32 mask.
mod common;

use sqisign_ec::{lift_basis, JacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/lift_basis.json"
);

#[test]
fn lift_basis_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::lift_basis");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let mut b = common::ec_basis_from("b_in", &v.inputs);
        let mut e = common::ec_curve_from("e_in", &v.inputs);
        let exp_p = common::jac_point_from("p", &v.outputs);
        let exp_q = common::jac_point_from("q", &v.outputs);
        let exp_b = common::ec_basis_from("b_out", &v.outputs);
        let exp_e = common::ec_curve_from("e_out", &v.outputs);
        let exp_ret = common::u32_field("ret", &v.outputs);
        let mut p = JacPoint::zero();
        let mut q = JacPoint::zero();
        let got_ret = lift_basis(&mut p, &mut q, &mut b, &mut e);
        assert_eq!(p, exp_p, "vector {} p diverged", v.id);
        assert_eq!(q, exp_q, "vector {} q diverged", v.id);
        assert_eq!(b, exp_b, "vector {} basis diverged", v.id);
        assert_eq!(e, exp_e, "vector {} curve diverged", v.id);
        assert_eq!(got_ret, exp_ret, "vector {} ret diverged", v.id);
    }
}
