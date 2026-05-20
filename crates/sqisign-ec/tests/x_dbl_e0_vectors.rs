//! Differential test of `xDBL_E0` (ported as `x_dbl_e0`).
mod common;

use sqisign_ec::{x_dbl_e0, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xDBL_E0.json");

#[test]
fn x_dbl_e0_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xDBL_E0");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("a", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        x_dbl_e0(&mut got, &p);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
