//! Differential test of `xeval_4`.
mod common;

use sqisign_ec::{xeval_4, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xeval_4.json");

#[test]
fn xeval_4_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xeval_4");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let q = common::ec_point_from("q", &v.inputs);
        let kps = common::ec_kps4_from("kps", &v.inputs);
        let exp = common::ec_point_from("c", &v.outputs);
        let q_arr = [q];
        let mut got_arr = [EcPoint::zero()];
        xeval_4(&mut got_arr, &q_arr, &kps);
        assert_eq!(got_arr[0], exp, "vector {} diverged", v.id);
    }
}
