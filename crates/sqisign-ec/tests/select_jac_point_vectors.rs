//! Differential test of `select_jac_point`.
mod common;

use sqisign_ec::select_jac_point;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/select_jac_point.json"
);

#[test]
fn select_jac_point_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::select_jac_point");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let mut got = common::jac_point_from("q_in", &v.inputs);
        let p1 = common::jac_point_from("p1", &v.inputs);
        let p2 = common::jac_point_from("p2", &v.inputs);
        let option = common::u64_field("option", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        select_jac_point(&mut got, &p1, &p2, option);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
