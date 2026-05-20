//! Differential test of `jac_to_xz`.
mod common;

use sqisign_ec::{jac_to_xz, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/jac_to_xz.json"
);

#[test]
fn jac_to_xz_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_to_xz");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let xy = common::jac_point_from("a", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        jac_to_xz(&mut got, &xy);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
