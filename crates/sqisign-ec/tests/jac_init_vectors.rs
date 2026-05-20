//! Differential test of `jac_init`.
mod common;

use sqisign_ec::jac_init;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/jac_init.json"
);

#[test]
fn jac_init_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_init");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let pre = common::jac_point_from("a", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        let mut got = pre;
        jac_init(&mut got);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
