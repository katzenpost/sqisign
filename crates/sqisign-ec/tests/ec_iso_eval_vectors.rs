//! Differential test of `ec_iso_eval`. Inputs: a point and an
//! `ec_isom_t`. Output: the mutated point.
mod common;

use sqisign_ec::ec_iso_eval;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/ec_iso_eval.json"
);

#[test]
fn ec_iso_eval_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::ec_iso_eval");
    assert!(!file.vectors.is_empty());

    for v in &file.vectors {
        let p_in = common::ec_point_from("p_in", &v.inputs);
        let isom = common::ec_isom_from("isom", &v.inputs);
        let exp = common::ec_point_from("p_out", &v.outputs);
        let mut p = p_in;
        ec_iso_eval(&mut p, &isom);
        assert_eq!(p, exp, "vector {} diverged", v.id);
    }
}
