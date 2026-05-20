//! Differential test of `add_couple_jac_points`.
mod common;

use sqisign_hd::{add_couple_jac_points, ThetaCoupleJacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/add_couple_jac_points.json"
);

#[test]
fn add_couple_jac_points_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::add_couple_jac_points");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let t1 = common::theta_couple_jac_point_from("t1", &v.inputs);
        let t2 = common::theta_couple_jac_point_from("t2", &v.inputs);
        let e12 = common::theta_couple_curve_from("e12", &v.inputs);
        let exp = common::theta_couple_jac_point_from("out", &v.outputs);
        let mut got = ThetaCoupleJacPoint::zero();
        add_couple_jac_points(&mut got, &t1, &t2, &e12);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
