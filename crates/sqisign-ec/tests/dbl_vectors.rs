//! Differential test of `DBL` (Jacobian doubling).
mod common;

use sqisign_ec::{dbl, JacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/DBL.json");

#[test]
fn dbl_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::DBL");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("p", &v.inputs);
        let ac = common::ec_curve_from("ac", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        let mut got = JacPoint::zero();
        dbl(&mut got, &p, &ac);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
