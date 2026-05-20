//! Differential test of `change_of_basis_matrix_tate_invert`.

mod common;

use sqisign_id2iso::change_of_basis_matrix_tate_invert;
use sqisign_quaternion::dim2::ibz_mat_2x2_new;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/id2iso/change_of_basis_matrix_tate_invert.json"
);

#[test]
fn change_of_basis_matrix_tate_invert_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_id2iso::change_of_basis_matrix_tate_invert"
    );
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let b1 = common::ec_basis_from("b1", &v.inputs);
        let b2 = common::ec_basis_from("b2", &v.inputs);
        let mut e = common::ec_curve_from("e_in", &v.inputs);
        let fparam = common::read_i32("f", &v.inputs);

        let exp_mat = common::read_mat2x2("mat", &v.outputs);
        let exp_e = common::ec_curve_from("e_out", &v.outputs);

        let mut mat = ibz_mat_2x2_new();
        change_of_basis_matrix_tate_invert(&mut mat, &b1, &b2, &mut e, fparam);
        assert!(
            common::mat2_eq(&mat, &exp_mat),
            "vector {} matrix diverged",
            v.id
        );
        assert_eq!(e, exp_e, "vector {} curve diverged", v.id);
    }
}
