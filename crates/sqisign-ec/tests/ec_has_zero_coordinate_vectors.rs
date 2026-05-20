//! Differential test of `ec_has_zero_coordinate`.
mod common;

use sqisign_ec::ec_has_zero_coordinate;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_has_zero_coordinate.json"
);

#[test]
fn ec_has_zero_coordinate_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_has_zero_coordinate");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("a", &v.inputs);
        let exp = common::u32_field("result", &v.outputs);
        let got = ec_has_zero_coordinate(&p);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
