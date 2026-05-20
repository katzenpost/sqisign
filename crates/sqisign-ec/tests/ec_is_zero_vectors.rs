//! Differential test of `ec_is_zero`.
mod common;

use sqisign_ec::ec_is_zero;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_is_zero.json"
);

#[test]
fn ec_is_zero_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_is_zero");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("a", &v.inputs);
        let exp = common::u32_field("result", &v.outputs);
        let got = ec_is_zero(&p);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
