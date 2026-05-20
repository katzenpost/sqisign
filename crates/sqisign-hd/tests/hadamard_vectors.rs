//! Differential test of `hadamard`.
mod common;

use sqisign_hd::{hadamard, ThetaPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/hadamard.json"
);

#[test]
fn hadamard_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::hadamard");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let a = common::theta_point_from("a", &v.inputs);
        let exp = common::theta_point_from("c", &v.outputs);
        let mut got = ThetaPoint::zero();
        hadamard(&mut got, &a);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
