//! Differential test of `xADD` (ported as `x_add`).
mod common;

use sqisign_ec::{x_add, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xADD.json");

#[test]
fn x_add_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xADD");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let q = common::ec_point_from("q", &v.inputs);
        let pq = common::ec_point_from("pq", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        x_add(&mut got, &p, &q, &pq);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
