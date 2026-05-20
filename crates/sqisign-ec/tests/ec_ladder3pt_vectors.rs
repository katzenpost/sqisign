//! Differential test of `ec_ladder3pt`.
mod common;

use sqisign_ec::{ec_ladder3pt, EcPoint, NWORDS_ORDER};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_ladder3pt.json"
);

#[test]
fn ec_ladder3pt_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_ladder3pt");
    assert!(file.vectors.len() >= 100);

    for v in &file.vectors {
        let m = common::digits_field("m", NWORDS_ORDER, &v.inputs);
        let p = common::ec_point_from("p", &v.inputs);
        let q = common::ec_point_from("q", &v.inputs);
        let pq = common::ec_point_from("pq", &v.inputs);
        let e = common::ec_curve_from("e", &v.inputs);
        let exp_c = common::ec_point_from("c", &v.outputs);
        let exp_r = common::i32_field("result", &v.outputs);
        let mut got = EcPoint::zero();
        let got_r = ec_ladder3pt(&mut got, &m, &p, &q, &pq, &e);
        assert_eq!(got_r, exp_r, "vector {} return diverged", v.id);
        if exp_r == 1 {
            assert_eq!(got, exp_c, "vector {} output diverged", v.id);
        }
    }
}
