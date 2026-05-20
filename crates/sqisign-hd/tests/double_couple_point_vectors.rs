//! Differential test of `double_couple_point`.
mod common;

use sqisign_hd::{double_couple_point, ThetaCouplePoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/double_couple_point.json"
);

#[test]
fn double_couple_point_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::double_couple_point");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let inp = common::theta_couple_point_from("in", &v.inputs);
        let e12 = common::theta_couple_curve_from("e12", &v.inputs);
        let exp = common::theta_couple_point_from("out", &v.outputs);
        let mut got = ThetaCouplePoint::zero();
        double_couple_point(&mut got, &inp, &e12);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
