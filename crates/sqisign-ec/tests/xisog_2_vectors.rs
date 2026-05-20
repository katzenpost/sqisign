//! Differential test of `xisog_2`.
mod common;

use sqisign_ec::{xisog_2, EcKps2, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/ec/xisog_2.json");

#[test]
fn xisog_2_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xisog_2");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let p = common::ec_point_from("p", &v.inputs);
        let exp_kps = common::ec_kps2_from("kps", &v.outputs);
        let exp_b = common::ec_point_from("b", &v.outputs);
        let mut kps = EcKps2::zero();
        let mut b = EcPoint::zero();
        xisog_2(&mut kps, &mut b, p);
        assert_eq!(kps, exp_kps, "vector {} kps diverged", v.id);
        assert_eq!(b, exp_b, "vector {} B diverged", v.id);
    }
}
