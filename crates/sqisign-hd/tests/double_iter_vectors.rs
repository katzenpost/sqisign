//! Differential test of `double_iter`.
mod common;

use sqisign_hd::{double_iter, ThetaPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/double_iter.json"
);

#[test]
fn double_iter_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_hd::double_iter");
    assert!(file.vectors.len() >= 500);

    for v in &file.vectors {
        let mut a_in = common::theta_structure_from("a_in", &v.inputs);
        let inp = common::theta_point_from("in", &v.inputs);
        let exp = common::theta_point_from("out", &v.outputs);
        let exp_a = common::theta_structure_from("a_out", &v.outputs);
        let exp_bits = common::u32_field("exp", &v.inputs) as i32;
        a_in.precomputation = false;
        let mut got = ThetaPoint::zero();
        double_iter(&mut got, &mut a_in, &inp, exp_bits);
        assert_eq!(got, exp, "vector {} out diverged", v.id);
        assert_eq!(a_in, exp_a, "vector {} a_out diverged", v.id);
    }
}
