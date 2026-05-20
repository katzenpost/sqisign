//! Differential test of `ec_isomorphism`. Inputs: two curves. Outputs:
//! the `ec_isom_t` (Nx, Nz, D) and the u32 error mask.
mod common;

use sqisign_ec::{ec_isomorphism, EcIsom};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_isomorphism.json"
);

#[test]
fn ec_isomorphism_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_isomorphism");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let from = common::ec_curve_from("from", &v.inputs);
        let to = common::ec_curve_from("to", &v.inputs);
        let exp_isom = common::ec_isom_from("isom", &v.outputs);
        let exp_ret = common::u32_field("ret", &v.outputs);
        let mut isom = EcIsom::zero();
        let got_ret = ec_isomorphism(&mut isom, &from, &to);
        assert_eq!(isom, exp_isom, "vector {} isom diverged", v.id);
        assert_eq!(got_ret, exp_ret, "vector {} ret diverged", v.id);
    }
}
