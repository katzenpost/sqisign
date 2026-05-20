//! Differential test of `ADD` (Jacobian addition).
mod common;

use sqisign_ec::{add, JacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/ADD.json");

#[test]
fn add_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ADD");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::jac_point_from("p", &v.inputs);
        let q = common::jac_point_from("q", &v.inputs);
        let ac = common::ec_curve_from("ac", &v.inputs);
        let exp = common::jac_point_from("c", &v.outputs);
        let mut got = JacPoint::zero();
        add(&mut got, &p, &q, &ac);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
