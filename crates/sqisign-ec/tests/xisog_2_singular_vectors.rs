//! Differential test of `xisog_2_singular`.
mod common;

use sqisign_ec::{xisog_2_singular, EcKps2, EcPoint};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/xisog_2_singular.json"
);

#[test]
fn xisog_2_singular_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::xisog_2_singular");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let a24 = common::ec_point_from("a24", &v.inputs);
        let exp_kps = common::ec_kps2_from("kps", &v.outputs);
        let exp_b24 = common::ec_point_from("b24", &v.outputs);
        let mut kps = EcKps2::zero();
        let mut b24 = EcPoint::zero();
        xisog_2_singular(&mut kps, &mut b24, a24);
        assert_eq!(kps, exp_kps, "vector {} kps diverged", v.id);
        assert_eq!(b24, exp_b24, "vector {} B24 diverged", v.id);
    }
}
