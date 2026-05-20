//! Differential test of `xDBLMUL` (ported as `x_dbl_mul`).
mod common;

use sqisign_ec::{x_dbl_mul, EcPoint, NWORDS_ORDER};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xDBLMUL.json");

#[test]
fn x_dbl_mul_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xDBLMUL");
    assert!(file.vectors.len() >= 100);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let k = common::digits_field("k", NWORDS_ORDER, &v.inputs);
        let q = common::ec_point_from("q", &v.inputs);
        let l = common::digits_field("l", NWORDS_ORDER, &v.inputs);
        let pq = common::ec_point_from("pq", &v.inputs);
        let kbits = common::i32_field("kbits", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp_c = common::ec_point_from("c", &v.outputs);
        let exp_r = common::i32_field("result", &v.outputs);
        let mut got = EcPoint::zero();
        let got_r = x_dbl_mul(&mut got, &p, &k, &q, &l, &pq, kbits, &curve);
        assert_eq!(got_r, exp_r, "vector {} return diverged", v.id);
        if exp_r == 1 {
            assert_eq!(got, exp_c, "vector {} output diverged", v.id);
        }
    }
}
