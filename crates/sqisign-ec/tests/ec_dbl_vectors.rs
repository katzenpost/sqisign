//! Differential test of `ec_dbl`.
mod common;

use sqisign_ec::{ec_dbl, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/ec_dbl.json");

#[test]
fn ec_dbl_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_dbl");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        ec_dbl(&mut got, &p, &curve);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
