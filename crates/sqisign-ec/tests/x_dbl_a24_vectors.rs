//! Differential test of `xDBL_A24` (ported as `x_dbl_a24`).
mod common;

use sqisign_ec::{x_dbl_a24, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/xDBL_A24.json"
);

#[test]
fn x_dbl_a24_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xDBL_A24");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let a24 = common::ec_point_from("a24", &v.inputs);
        let norm = common::u32_field("a24_normalized", &v.inputs) != 0;
        let exp = common::ec_point_from("c", &v.outputs);
        let mut got = EcPoint::zero();
        x_dbl_a24(&mut got, &p, &a24, norm);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
