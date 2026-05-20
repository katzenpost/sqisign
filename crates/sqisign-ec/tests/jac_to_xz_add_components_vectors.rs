//! Differential test of `jac_to_xz_add_components`.
mod common;

use sqisign_ec::{jac_to_xz_add_components, AddComponents};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/jac_to_xz_add_components.json"
);

#[test]
fn jac_to_xz_add_components_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_to_xz_add_components");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("p", &v.inputs);
        let q = common::jac_point_from("q", &v.inputs);
        let ac = common::ec_curve_from("ac", &v.inputs);
        let exp = common::add_components_from("c", &v.outputs);
        let mut got = AddComponents::zero();
        jac_to_xz_add_components(&mut got, &p, &q, &ac);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
