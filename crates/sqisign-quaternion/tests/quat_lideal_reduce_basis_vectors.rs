//! Differential test of `quat_lideal_reduce_basis`.
//!
//! The reduced basis is non-unique; we assert structural validity and
//! report bit-exact agreement. The gram matrix the C reference returns is
//! scaled by `denom^2`, has its diagonal halved, and its strict upper
//! triangle zeroed; the Rust port mirrors this exactly.

mod common;

use common::{ibz_eq, ibz_mat_4x4_new_local, mat4x4_eq, read_i32, read_ibz, read_mat4x4};
use sqisign_quaternion::{
    ibz_set, quat_alg_elem_copy_ibz, quat_lattice_O0_set, quat_lideal_create_principal,
    quat_lideal_reduce_basis, quat_lll_set_ibq_parameters, quat_lll_verify, Ibq, Ibz, QuatAlg,
    QuatAlgElem, QuatLattice, QuatLeftIdeal,
};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lideal_reduce_basis.json"
);

fn read_quat_elem(
    prefix: &str,
    inputs: &std::collections::BTreeMap<String, String>,
) -> QuatAlgElem {
    let denom = read_ibz(&format!("{prefix}_denom"), inputs);
    let c0 = read_ibz(&format!("{prefix}_c0"), inputs);
    let c1 = read_ibz(&format!("{prefix}_c1"), inputs);
    let c2 = read_ibz(&format!("{prefix}_c2"), inputs);
    let c3 = read_ibz(&format!("{prefix}_c3"), inputs);
    let mut e = QuatAlgElem::new();
    quat_alg_elem_copy_ibz(&mut e, &denom, &c0, &c1, &c2, &c3);
    e
}

#[test]
fn quat_lideal_reduce_basis_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lideal_reduce_basis");

    let mut delta = Ibq::new();
    let mut eta = Ibq::new();
    quat_lll_set_ibq_parameters(&mut delta, &mut eta);

    let mut bit_exact_basis = 0usize;
    let mut bit_exact_gram = 0usize;
    let mut total = 0usize;

    for v in &f.vectors {
        total += 1;
        let x = read_quat_elem("x", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let alg = QuatAlg::init_set(&p);

        let red_exp = read_mat4x4("red", &v.outputs);
        let gram_exp = read_mat4x4("gram", &v.outputs);
        let valid_expected = read_i32("valid", &v.outputs);

        // Build O0 as the parent maximal order, build the principal ideal,
        // then reduce its basis.
        let mut o0 = QuatLattice::new();
        quat_lattice_O0_set(&mut o0);
        let mut lideal = QuatLeftIdeal::new();
        quat_lideal_create_principal(&mut lideal, &x, &o0, &alg);

        let mut red = ibz_mat_4x4_new_local();
        let mut gram = ibz_mat_4x4_new_local();
        quat_lideal_reduce_basis(&mut red, &mut gram, &lideal, &alg);

        let valid_rust = quat_lll_verify(&red, &delta, &eta, &alg);
        assert_eq!(valid_rust, valid_expected, "vector {}: validity bit", v.id);
        assert_eq!(
            valid_rust, 1,
            "vector {}: reduced basis was already invalid",
            v.id
        );

        if mat4x4_eq(&red, &red_exp) {
            bit_exact_basis += 1;
        }
        if mat4x4_eq(&gram, &gram_exp) {
            bit_exact_gram += 1;
        }
        // Sanity: gram is lower-triangular (strict upper is zero).
        #[allow(clippy::needless_range_loop)]
        for i in 0..4 {
            for j in (i + 1)..4 {
                let mut zero = Ibz::zero();
                ibz_set(&mut zero, 0);
                assert!(
                    ibz_eq(&gram[i][j], &zero),
                    "vector {}: gram[{}][{}] should be zero",
                    v.id,
                    i,
                    j
                );
            }
        }
    }
    eprintln!(
        "quat_lideal_reduce_basis: bit-exact basis {}/{}, bit-exact gram {}/{}",
        bit_exact_basis, total, bit_exact_gram, total
    );
}
