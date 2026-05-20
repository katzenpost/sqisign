//! Differential test of `double_point`.
mod common;

use sqisign_hd::{double_point, ThetaPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/double_point.json"
);

#[test]
fn double_point_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::double_point");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let mut a_in = common::theta_structure_from("a_in", &v.inputs);
        let inp = common::theta_point_from("in", &v.inputs);
        let exp_a = common::theta_structure_from("a_out", &v.outputs);
        let exp_out = common::theta_point_from("out", &v.outputs);
        // The harness forces precomputation=false to exercise the
        // precompute branch on every call; preserve that here.
        a_in.precomputation = false;
        let mut got = ThetaPoint::zero();
        double_point(&mut got, &mut a_in, &inp);
        assert_eq!(got, exp_out, "vector {} out diverged", v.id);
        assert_eq!(a_in, exp_a, "vector {} a_out diverged", v.id);
    }
}
