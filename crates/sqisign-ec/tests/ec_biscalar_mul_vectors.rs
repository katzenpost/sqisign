//! Differential test of `ec_biscalar_mul`.
mod common;

use sqisign_ec::{ec_biscalar_mul, EcPoint, NWORDS_ORDER};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_biscalar_mul.json"
);

#[test]
fn ec_biscalar_mul_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_biscalar_mul");
    assert!(file.vectors.len() >= 100);

    for v in &file.vectors {
        let sp = common::digits_field("scalar_p", NWORDS_ORDER, &v.inputs);
        let sq = common::digits_field("scalar_q", NWORDS_ORDER, &v.inputs);
        let kbits = common::i32_field("kbits", &v.inputs);
        let b = common::ec_basis_from("b", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp_c = common::ec_point_from("c", &v.outputs);
        let exp_r = common::i32_field("result", &v.outputs);
        let mut got = EcPoint::zero();
        let got_r = ec_biscalar_mul(&mut got, &sp, &sq, kbits, &b, &curve);
        assert_eq!(got_r, exp_r, "vector {} return diverged", v.id);
        if exp_r == 1 {
            assert_eq!(got, exp_c, "vector {} output diverged", v.id);
        }
    }
}
