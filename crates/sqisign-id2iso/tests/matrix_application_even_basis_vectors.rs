//! Differential test of `matrix_application_even_basis`.

mod common;

use sqisign_id2iso::matrix_application_even_basis;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/id2iso/matrix_application_even_basis.json"
);

#[test]
fn matrix_application_even_basis_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_id2iso::matrix_application_even_basis");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let mut pq = common::ec_basis_from("pq_in", &v.inputs);
        let e = common::ec_curve_from("e", &v.inputs);
        let mut mat = common::read_mat2x2("mat_in", &v.inputs);
        let fparam = common::read_i32("f", &v.inputs);

        let exp_pq = common::ec_basis_from("pq_out", &v.outputs);
        let exp_mat = common::read_mat2x2("mat_out", &v.outputs);
        let exp_ok = common::read_i32("ok", &v.outputs);

        let ok = matrix_application_even_basis(&mut pq, &e, &mut mat, fparam);
        assert_eq!(ok, exp_ok, "vector {} ok bit", v.id);
        assert_eq!(pq, exp_pq, "vector {} basis diverged", v.id);
        assert!(
            common::mat2_eq(&mat, &exp_mat),
            "vector {} matrix diverged",
            v.id
        );
    }
}
