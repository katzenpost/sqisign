//! Differential test of `quat_lattice_bound_parallelogram`.
//!
//! The dual-Gram reduction inside relies on `quat_lll_core`, so the box
//! and the change-of-basis matrix `U` are non-unique. We assert two
//! contracts:
//!
//!  * the `ok` (non-trivial) bit matches the reference;
//!  * the absolute determinant of `U` is 1 (the C reference asserts the
//!    same in debug builds);
//!
//! and report bit-exact agreement for visibility.

mod common;

use common::{
    ibz_mat_4x4_new_local, mat4x4_eq, read_i32, read_ibz, read_mat4x4, read_vec4, vec4_eq,
};
use sqisign_quaternion::{
    ibz_abs, ibz_mat_4x4_inv_with_det_as_denom, ibz_vec_4_new, quat_lattice_bound_parallelogram,
    Ibz,
};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_bound_parallelogram.json"
);

#[test]
fn quat_lattice_bound_parallelogram_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::quat_lattice_bound_parallelogram"
    );

    let mut bit_exact_box = 0usize;
    let mut bit_exact_u = 0usize;
    let mut total = 0usize;

    for v in &f.vectors {
        total += 1;
        let g = read_mat4x4("G", &v.inputs);
        let radius = read_ibz("radius", &v.inputs);

        let box_exp = read_vec4("box", &v.outputs);
        let u_exp = read_mat4x4("U", &v.outputs);
        let ok_expected = read_i32("ok", &v.outputs);

        let mut box_actual = ibz_vec_4_new();
        let mut u_actual = ibz_mat_4x4_new_local();
        let ok_rust = quat_lattice_bound_parallelogram(&mut box_actual, &mut u_actual, &g, &radius);
        assert_eq!(ok_rust, ok_expected, "vector {}: ok bit", v.id);
        assert_eq!(ok_rust, 1, "vector {}: ok bit was unexpectedly 0", v.id);

        // |det(U)| must equal 1 (U is unitary).
        let mut det = Ibz::zero();
        let _ = ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det), &u_actual);
        let mut det_abs = Ibz::zero();
        ibz_abs(&mut det_abs, &det);
        assert!(
            det_abs.0 == num_bigint::BigInt::from(1),
            "vector {}: |det(U)| = {} (must be 1)",
            v.id,
            det_abs.0
        );

        if vec4_eq(&box_actual, &box_exp) {
            bit_exact_box += 1;
        }
        if mat4x4_eq(&u_actual, &u_exp) {
            bit_exact_u += 1;
        }
    }
    eprintln!(
        "quat_lattice_bound_parallelogram: bit-exact box {}/{}, U {}/{}",
        bit_exact_box, total, bit_exact_u, total
    );
}
