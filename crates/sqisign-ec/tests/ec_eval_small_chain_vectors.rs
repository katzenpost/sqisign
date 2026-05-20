//! Differential test of `ec_eval_small_chain`. Inputs: curve (mutated),
//! kernel, len, two points (mutated), special. Outputs: mutated curve,
//! mutated points, u32 mask.
mod common;

use sqisign_ec::ec_eval_small_chain;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_eval_small_chain.json"
);

#[test]
fn ec_eval_small_chain_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_eval_small_chain");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let mut curve = common::ec_curve_from("curve_in", &v.inputs);
        let kernel = common::ec_point_from("kernel", &v.inputs);
        let len = common::i32_field("len", &v.inputs);
        let p0_in = common::ec_point_from("p0_in", &v.inputs);
        let p1_in = common::ec_point_from("p1_in", &v.inputs);
        let special_bytes = decode("special", &v.inputs["special"]).expect("special hex");
        assert_eq!(special_bytes.len(), 1);
        let special = special_bytes[0] != 0;

        let exp_curve = common::ec_curve_from("curve_out", &v.outputs);
        let exp_p0 = common::ec_point_from("p0_out", &v.outputs);
        let exp_p1 = common::ec_point_from("p1_out", &v.outputs);
        let exp_ret = common::u32_field("ret", &v.outputs);

        let mut points = [p0_in, p1_in];
        let got_ret = ec_eval_small_chain(&mut curve, &kernel, len, &mut points, special);
        assert_eq!(curve, exp_curve, "vector {} curve diverged", v.id);
        assert_eq!(points[0], exp_p0, "vector {} p0 diverged", v.id);
        assert_eq!(points[1], exp_p1, "vector {} p1 diverged", v.id);
        assert_eq!(got_ret, exp_ret, "vector {} ret diverged", v.id);
    }
}
