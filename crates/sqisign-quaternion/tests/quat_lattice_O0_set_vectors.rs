//! Differential test of `quat_lattice_O0_set`.
#![allow(non_snake_case)]
mod common;
use common::{ibz_eq, mat4x4_eq, read_lattice};
use sqisign_quaternion::{quat_lattice_O0_set, QuatLattice};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_O0_set.json"
);

#[test]
fn quat_lattice_O0_set_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_O0_set");
    assert_eq!(f.vectors.len(), 1);
    let v = &f.vectors[0];
    let exp = read_lattice("O0", &v.outputs);
    let mut o0 = QuatLattice::new();
    quat_lattice_O0_set(&mut o0);
    assert!(ibz_eq(&o0.denom, &exp.denom), "denom");
    assert!(mat4x4_eq(&o0.basis, &exp.basis), "basis");
}
