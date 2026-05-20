//! Differential test of `pointwise_square`.
mod common;

use sqisign_hd::{pointwise_square, ThetaPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/pointwise_square.json"
);

#[test]
fn pointwise_square_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::pointwise_square");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let a = common::theta_point_from("a", &v.inputs);
        let exp = common::theta_point_from("c", &v.outputs);
        let mut got = ThetaPoint::zero();
        pointwise_square(&mut got, &a);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
