//! Differential test of `ec_point_init` against the recorded C
//! battery. Setter shape: prefill ec_point -> ec_point.
mod common;

use sqisign_ec::{ec_point_init, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_point_init.json"
);

#[test]
fn ec_point_init_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_ec::ec_point_init");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let pre = common::ec_point_from("a", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got: EcPoint = pre;
        ec_point_init(&mut got);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
