//! Differential test of `endomorphism_application_even_basis`.

mod common;

use sqisign_id2iso::endomorphism_application_even_basis;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/id2iso/endomorphism_application_even_basis.json"
);

#[test]
fn endomorphism_application_even_basis_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_id2iso::endomorphism_application_even_basis"
    );
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let mut pq = common::ec_basis_from("pq_in", &v.inputs);
        let e = common::ec_curve_from("e", &v.inputs);
        let idx = common::read_i32("idx", &v.inputs);
        let theta = common::read_quat_elem("theta", &v.inputs);
        let fparam = common::read_i32("f", &v.inputs);

        let exp_pq = common::ec_basis_from("pq_out", &v.outputs);

        endomorphism_application_even_basis(&mut pq, idx, &e, &theta, fparam);
        assert_eq!(pq, exp_pq, "vector {} basis diverged", v.id);
    }
}
