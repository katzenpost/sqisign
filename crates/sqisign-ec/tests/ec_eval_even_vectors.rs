//! Differential test of `ec_eval_even`. Inputs: an `ec_isog_even_t`
//! and two points (mutated). Outputs: image curve, mutated points,
//! u32 mask.
mod common;

use sqisign_ec::{ec_eval_even, EcCurve};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_eval_even.json"
);

#[test]
fn ec_eval_even_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_eval_even");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let phi = common::ec_isog_even_from("phi_in", &v.inputs);
        let p0_in = common::ec_point_from("p0_in", &v.inputs);
        let p1_in = common::ec_point_from("p1_in", &v.inputs);
        let exp_image = common::ec_curve_from("image", &v.outputs);
        let exp_p0 = common::ec_point_from("p0_out", &v.outputs);
        let exp_p1 = common::ec_point_from("p1_out", &v.outputs);
        let exp_ret = common::u32_field("ret", &v.outputs);

        let mut image = EcCurve::zero();
        let mut points = [p0_in, p1_in];
        let got_ret = ec_eval_even(&mut image, &phi, &mut points);
        assert_eq!(image, exp_image, "vector {} image diverged", v.id);
        assert_eq!(points[0], exp_p0, "vector {} p0 diverged", v.id);
        assert_eq!(points[1], exp_p1, "vector {} p1 diverged", v.id);
        assert_eq!(got_ret, exp_ret, "vector {} ret diverged", v.id);
    }
}
