//! Differential test of `ec_j_inv`.
mod common;

use sqisign_ec::ec_j_inv;
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_j_inv.json"
);

#[test]
fn ec_j_inv_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_j_inv");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let e = common::ec_curve_from("e", &v.inputs);
        let exp = common::fp2_from("c", &v.outputs);
        let mut got = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        ec_j_inv(&mut got, &e);
        assert_eq!(got, exp, "vector {} diverged", v.id);
    }
}
