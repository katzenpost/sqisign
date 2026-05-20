//! Differential test of `ec_curve_init_from_A` (ported as
//! `ec_curve_init_from_a`).
mod common;

use sqisign_ec::ec_curve_init_from_a;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_curve_init_from_A.json"
);

#[test]
fn ec_curve_init_from_a_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_curve_init_from_A");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let mut got = common::ec_curve_from("a_curve", &v.inputs);
        let a = common::fp2_from("a", &v.inputs);
        let exp_curve = common::ec_curve_from("c", &v.outputs);
        let exp_r = common::i32_field("result", &v.outputs);
        let got_r = ec_curve_init_from_a(&mut got, &a);
        assert_eq!(got, exp_curve, "vector {} curve diverged", v.id);
        assert_eq!(got_r, exp_r, "vector {} return diverged", v.id);
    }
}
