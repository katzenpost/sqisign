//! Differential test of `copy_jac_point`.
mod common;

use sqisign_ec::{copy_jac_point, JacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/copy_jac_point.json"
);

#[test]
fn copy_jac_point_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::copy_jac_point");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let src = common::jac_point_from("a", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        let mut got = JacPoint::zero();
        copy_jac_point(&mut got, &src);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
