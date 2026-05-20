//! Differential test of `couple_jac_to_xz`.
mod common;

use sqisign_hd::{couple_jac_to_xz, ThetaCouplePoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/couple_jac_to_xz.json"
);

#[test]
fn couple_jac_to_xz_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::couple_jac_to_xz");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let xy = common::theta_couple_jac_point_from("xyP", &v.inputs);
        let exp = common::theta_couple_point_from("P", &v.outputs);
        let mut got = ThetaCouplePoint::zero();
        couple_jac_to_xz(&mut got, &xy);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
