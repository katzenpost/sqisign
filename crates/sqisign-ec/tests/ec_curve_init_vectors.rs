//! Differential test of `ec_curve_init`. Setter shape: prefill
//! ec_curve -> ec_curve.
mod common;

use sqisign_ec::{ec_curve_init, EcCurve};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_curve_init.json"
);

#[test]
fn ec_curve_init_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_curve_init");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let pre = common::ec_curve_from("a", &v.inputs);
        let exp = common::ec_curve_from("c", &v.outputs);
        let mut got: EcCurve = pre;
        ec_curve_init(&mut got);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
