//! Differential test of `ec_curve_verify_A` (ported as
//! `ec_curve_verify_a`).
mod common;

use sqisign_ec::ec_curve_verify_a;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_curve_verify_A.json"
);

#[test]
fn ec_curve_verify_a_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_curve_verify_A");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let a = common::fp2_from("a", &v.inputs);
        let exp = common::i32_field("result", &v.outputs);
        let got = ec_curve_verify_a(&a);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
