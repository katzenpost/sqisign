//! Differential test of `ec_is_two_torsion`.
mod common;

use sqisign_ec::ec_is_two_torsion;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_is_two_torsion.json"
);

#[test]
fn ec_is_two_torsion_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_is_two_torsion");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let e = common::ec_curve_from("e", &v.inputs);
        let exp = common::u32_field("result", &v.outputs);
        let got = ec_is_two_torsion(&p, &e);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
