//! Differential test of `quat_lll_core`.
//!
//! L2 LLL reduction does not produce a unique output; many bases are
//! `(delta, eta)`-reduced for the same input lattice. The differential
//! contract is therefore structural validity, checked against the same
//! `quat_lll_verify` oracle the C harness uses to record the `valid` bit
//! on each vector. Bit-exact agreement with the C reference is reported
//! at the end as a bonus.

mod common;

use common::{mat4x4_eq, read_i32, read_ibz, read_mat4x4};
use sqisign_quaternion::{
    quat_lll_core, quat_lll_set_ibq_parameters, quat_lll_verify, Ibq, QuatAlg,
};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lll_core.json"
);

#[test]
fn quat_lll_core_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lll_core");

    let mut delta = Ibq::new();
    let mut eta = Ibq::new();
    quat_lll_set_ibq_parameters(&mut delta, &mut eta);

    let mut bit_exact_basis = 0usize;
    let mut bit_exact_gram = 0usize;
    let mut total = 0usize;

    for v in &f.vectors {
        total += 1;
        let mut g = read_mat4x4("G", &v.inputs);
        let mut basis = read_mat4x4("basis", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let alg = QuatAlg::init_set(&p);

        let g_exp = read_mat4x4("G_out", &v.outputs);
        let basis_exp = read_mat4x4("basis_out", &v.outputs);
        let valid_expected = read_i32("valid", &v.outputs);

        quat_lll_core(&mut g, &mut basis);

        let valid_rust = quat_lll_verify(&basis, &delta, &eta, &alg);
        assert_eq!(valid_rust, valid_expected, "vector {}: validity bit", v.id);
        assert_eq!(
            valid_rust, 1,
            "vector {}: reference produced an invalid LLL basis (impossible)",
            v.id
        );

        if mat4x4_eq(&basis, &basis_exp) {
            bit_exact_basis += 1;
        }
        if mat4x4_eq(&g, &g_exp) {
            bit_exact_gram += 1;
        }
    }
    // Bit-exact match rate is reported but not asserted: the test passes
    // when every Rust output is structurally valid.
    eprintln!(
        "quat_lll_core: bit-exact basis {}/{}, bit-exact gram {}/{}",
        bit_exact_basis, total, bit_exact_gram, total
    );
}
