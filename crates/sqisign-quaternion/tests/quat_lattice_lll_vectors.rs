//! Differential test of `quat_lattice_lll`.
//!
//! Like the `quat_lll_core` boundary, the output basis is not unique;
//! we assert structural validity via `quat_lll_verify` and report
//! bit-exact agreement with the C reference for visibility.

mod common;

use common::{ibz_mat_4x4_new_local, mat4x4_eq, read_i32, read_ibz, read_lattice, read_mat4x4};
use sqisign_quaternion::{
    quat_lattice_lll, quat_lll_set_ibq_parameters, quat_lll_verify, Ibq, QuatAlg,
};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_lll.json"
);

#[test]
fn quat_lattice_lll_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_lll");

    let mut delta = Ibq::new();
    let mut eta = Ibq::new();
    quat_lll_set_ibq_parameters(&mut delta, &mut eta);

    let mut bit_exact = 0usize;
    let mut total = 0usize;

    for v in &f.vectors {
        total += 1;
        let lat = read_lattice("lat", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let alg = QuatAlg::init_set(&p);

        let red_exp = read_mat4x4("red", &v.outputs);
        let valid_expected = read_i32("valid", &v.outputs);

        let mut red = ibz_mat_4x4_new_local();
        let _ = quat_lattice_lll(&mut red, &lat, &alg);

        let valid_rust = quat_lll_verify(&red, &delta, &eta, &alg);
        assert_eq!(valid_rust, valid_expected, "vector {}: validity bit", v.id);
        assert_eq!(
            valid_rust, 1,
            "vector {}: reference output was already invalid",
            v.id
        );

        if mat4x4_eq(&red, &red_exp) {
            bit_exact += 1;
        }
    }
    eprintln!(
        "quat_lattice_lll: bit-exact reduced basis {}/{}",
        bit_exact, total
    );
}
