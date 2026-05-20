//! Differential test of `ec_dlog_2_tate`. Same shape as
//! `ec_dlog_2_weil`, but PQ must be a full E[2^TORSION_EVEN_POWER]
//! basis and RS may have smaller order 2^e.
mod common;

use sqisign_ec::{ec_dlog_2_tate, NWORDS_ORDER};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_dlog_2_tate.json"
);

#[test]
fn ec_dlog_2_tate_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_dlog_2_tate");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let pq = common::ec_basis_from("pq_in", &v.inputs);
        let rs = common::ec_basis_from("rs", &v.inputs);
        let mut curve = common::ec_curve_from("e_in", &v.inputs);
        let e_pow = common::i32_field("e", &v.inputs);

        let exp_r1 = common::digits_field("r1", NWORDS_ORDER, &v.outputs);
        let exp_r2 = common::digits_field("r2", NWORDS_ORDER, &v.outputs);
        let exp_s1 = common::digits_field("s1", NWORDS_ORDER, &v.outputs);
        let exp_s2 = common::digits_field("s2", NWORDS_ORDER, &v.outputs);
        let exp_pq = common::ec_basis_from("pq_out", &v.outputs);
        let exp_e = common::ec_curve_from("e_out", &v.outputs);

        let mut r1 = vec![0u64; NWORDS_ORDER];
        let mut r2 = vec![0u64; NWORDS_ORDER];
        let mut s1 = vec![0u64; NWORDS_ORDER];
        let mut s2 = vec![0u64; NWORDS_ORDER];
        ec_dlog_2_tate(
            &mut r1, &mut r2, &mut s1, &mut s2, &pq, &rs, &mut curve, e_pow,
        );
        assert_eq!(r1, exp_r1, "vector {} r1 diverged", v.id);
        assert_eq!(r2, exp_r2, "vector {} r2 diverged", v.id);
        assert_eq!(s1, exp_s1, "vector {} s1 diverged", v.id);
        assert_eq!(s2, exp_s2, "vector {} s2 diverged", v.id);
        // ec_dlog_2_tate takes PQ as `const`; the recorded `pq_out` is
        // byte-identical to `pq_in` since the function never writes back.
        assert_eq!(pq, exp_pq, "vector {} basis diverged", v.id);
        assert_eq!(curve, exp_e, "vector {} curve diverged", v.id);
    }
}
