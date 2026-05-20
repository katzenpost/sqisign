//! Differential test of `jac_neg`.
mod common;

use sqisign_ec::{jac_neg, JacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/jac_neg.json");

#[test]
fn jac_neg_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_neg");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("a", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        let mut got = JacPoint::zero();
        jac_neg(&mut got, &p);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
