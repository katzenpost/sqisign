//! Differential test of `xMUL` (ported as `x_mul`).
mod common;

use sqisign_ec::{x_mul, EcPoint, NWORDS_ORDER};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xMUL.json");

#[test]
fn x_mul_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xMUL");
    assert!(file.vectors.len() >= 100);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let k = common::digits_field("k", NWORDS_ORDER, &v.inputs);
        let kbits = common::i32_field("kbits", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        x_mul(&mut got, &p, &k, kbits, &curve);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
