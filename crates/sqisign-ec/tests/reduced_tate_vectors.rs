//! Differential test of `reduced_tate`. Inputs: e (u32), P, Q, PQ
//! (ec_points), curve (mutated). Outputs: pairing value (fp2), mutated
//! curve.
mod common;

use sqisign_ec::reduced_tate;
use sqisign_gf::{Fp2, NWORDS_FIELD};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/reduced_tate.json"
);

#[test]
fn reduced_tate_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::reduced_tate");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let e_pow = common::u32_field("e", &v.inputs);
        let p = common::ec_point_from("p", &v.inputs);
        let q = common::ec_point_from("q", &v.inputs);
        let pq = common::ec_point_from("pq", &v.inputs);
        let mut e = common::ec_curve_from("e_in", &v.inputs);
        let exp_r = common::fp2_from("r", &v.outputs);
        let exp_e = common::ec_curve_from("e_out", &v.outputs);

        let mut r = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        reduced_tate(&mut r, e_pow, &p, &q, &pq, &mut e);
        assert_eq!(r, exp_r, "vector {} r diverged", v.id);
        assert_eq!(e, exp_e, "vector {} curve diverged", v.id);
    }
}
