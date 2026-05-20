//! Differential test of `ec_dbl_iter`.
mod common;

use sqisign_ec::{ec_dbl_iter, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_dbl_iter.json"
);

#[test]
fn ec_dbl_iter_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_dbl_iter");
    assert!(file.vectors.len() >= 100);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let n = common::i32_field("n", &v.inputs);
        let mut curve = common::ec_curve_from("curve_in", &v.inputs);
        let exp_c = common::ec_point_from("c", &v.outputs);
        let exp_curve = common::ec_curve_from("curve_out", &v.outputs);
        let mut got = EcPoint::zero();
        ec_dbl_iter(&mut got, n, &p, &mut curve);
        assert_eq!(got, exp_c, "vector {} output diverged", v.id);
        assert_eq!(curve, exp_curve, "vector {} curve diverged", v.id);
    }
}
