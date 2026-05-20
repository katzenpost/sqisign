//! Differential test of `ec_normalize_point`. In-place mutator wrapped
//! via copy-then-call.
mod common;

use sqisign_ec::ec_normalize_point;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_normalize_point.json"
);

#[test]
fn ec_normalize_point_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_normalize_point");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let pre = common::ec_point_from("a", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = pre;
        ec_normalize_point(&mut got);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
