//! Differential test of `xDBL` (ported as `x_dbl`).
mod common;

use sqisign_ec::{x_dbl, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xDBL.json");

#[test]
fn x_dbl_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xDBL");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let ac = common::ec_curve_from("ac", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        x_dbl(&mut got, &p, &ac);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
