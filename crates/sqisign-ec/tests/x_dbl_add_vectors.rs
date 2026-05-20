//! Differential test of `xDBLADD` (ported as `x_dbl_add`).
mod common;

use sqisign_ec::{x_dbl_add, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xDBLADD.json");

#[test]
fn x_dbl_add_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xDBLADD");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let q = common::ec_point_from("q", &v.inputs);
        let pq = common::ec_point_from("pq", &v.inputs);
        let a24 = common::ec_point_from("a24", &v.inputs);
        let norm = common::u32_field("a24_normalized", &v.inputs) != 0;
        let exp_r = common::ec_point_from("r", &v.outputs);
        let exp_s = common::ec_point_from("s", &v.outputs);
        let mut r = EcPoint::zero();
        let mut s = EcPoint::zero();
        x_dbl_add(&mut r, &mut s, &p, &q, &pq, &a24, norm);
        assert_eq!(r, exp_r, "vector {} R diverged", v.id);
        assert_eq!(s, exp_s, "vector {} S diverged", v.id);
    }
}
