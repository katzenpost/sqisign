//! Differential test of `quat_lattice_gram`.
mod common;
use common::{mat4x4_eq, read_ibz, read_lattice, read_mat4x4};
use sqisign_quaternion::{ibz_mat_4x4_new, quat_lattice_gram, QuatAlg};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_gram.json"
);

#[test]
fn quat_lattice_gram_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_gram");
    for v in &f.vectors {
        let l = read_lattice("l", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let exp = read_mat4x4("G", &v.outputs);
        let alg = QuatAlg::init_set(&p);
        let mut g = ibz_mat_4x4_new();
        quat_lattice_gram(&mut g, &l, &alg);
        assert!(mat4x4_eq(&g, &exp), "vector {}", v.id);
    }
}
