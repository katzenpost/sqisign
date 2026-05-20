//! Differential test of `double_couple_jac_point_iter`.
mod common;

use sqisign_hd::{double_couple_jac_point_iter, ThetaCoupleJacPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/double_couple_jac_point_iter.json"
);

#[test]
fn double_couple_jac_point_iter_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::double_couple_jac_point_iter");
    assert!(file.vectors.len() >= 200);

    for v in &file.vectors {
        let n = common::u32_field("n", &v.inputs);
        let inp = common::theta_couple_jac_point_from("in", &v.inputs);
        let e12 = common::theta_couple_curve_from("e12", &v.inputs);
        let exp = common::theta_couple_jac_point_from("out", &v.outputs);
        let mut got = ThetaCoupleJacPoint::zero();
        double_couple_jac_point_iter(&mut got, n, &inp, &e12);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
