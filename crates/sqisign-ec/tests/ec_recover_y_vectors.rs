//! Differential test of `ec_recover_y`. Inputs: an fp2 `px` and a
//! curve. Outputs: the recovered y (fp2) and the u32 mask.
mod common;

use sqisign_ec::ec_recover_y;
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_recover_y.json"
);

#[test]
fn ec_recover_y_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_recover_y");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let px = common::fp2_from("px", &v.inputs);
        let curve = common::ec_curve_from("curve", &v.inputs);
        let exp_y = common::fp2_from("y", &v.outputs);
        let exp_ret = common::u32_field("ret", &v.outputs);
        let mut y = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        let got_ret = ec_recover_y(&mut y, &px, &curve);
        assert_eq!(y, exp_y, "vector {} y diverged", v.id);
        assert_eq!(got_ret, exp_ret, "vector {} ret diverged", v.id);
    }
}
