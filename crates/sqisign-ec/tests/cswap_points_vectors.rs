//! Differential test of `cswap_points`.
mod common;

use sqisign_ec::cswap_points;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/ec/cswap_points.json"
);

#[test]
fn cswap_points_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_ec::cswap_points");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let mut p = common::ec_point_from("p_in", &v.inputs);
        let mut q = common::ec_point_from("q_in", &v.inputs);
        let option = common::u64_field("option", &v.inputs);
        let exp_p = common::ec_point_from("p_out", &v.outputs);
        let exp_q = common::ec_point_from("q_out", &v.outputs);
        cswap_points(&mut p, &mut q, option);
        assert_eq!(p, exp_p, "vector {} P diverged", v.id);
        assert_eq!(q, exp_q, "vector {} Q diverged", v.id);
    }
}
