//! Differential test of `jac_is_equal`.
mod common;

use sqisign_ec::jac_is_equal;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/jac_is_equal.json"
);

#[test]
fn jac_is_equal_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_is_equal");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("p", &v.inputs);
        let q = common::jac_point_from("q", &v.inputs);
        let exp = common::u32_field("result", &v.outputs);
        let got = jac_is_equal(&p, &q);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
