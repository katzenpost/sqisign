//! Differential test of `ec_curve_normalize_A24` (ported as
//! `ec_curve_normalize_a24`).
mod common;

use sqisign_ec::ec_curve_normalize_a24;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_curve_normalize_A24.json"
);

#[test]
fn ec_curve_normalize_a24_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_curve_normalize_A24");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let mut got = common::ec_curve_from("a", &v.inputs);
        let exp = common::ec_curve_from("c", &v.outputs);
        ec_curve_normalize_a24(&mut got);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
