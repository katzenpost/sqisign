//! Differential test of `jac_from_ws`.
//!
//! The reference writes `Q->y` and `Q->z` unconditionally, but only
//! writes `Q->x` when `curve->A` is non-zero; the recorded `q_in`
//! captures the prefill that survives the no-op case.
mod common;

use sqisign_ec::jac_from_ws;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/jac_from_ws.json"
);

#[test]
fn jac_from_ws_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::jac_from_ws");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let mut got = common::jac_point_from("q_in", &v.inputs);
        let p = common::jac_point_from("p", &v.inputs);
        let ao3 = common::fp2_from("ao3", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        jac_from_ws(&mut got, &p, &ao3, &curve);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
