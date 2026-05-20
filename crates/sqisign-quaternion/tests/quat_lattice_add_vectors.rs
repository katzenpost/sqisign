//! Differential test of `quat_lattice_add`.
mod common;
use common::{ibz_eq, mat4x4_eq, read_lattice};
use sqisign_quaternion::{quat_lattice_add, QuatLattice};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_add.json"
);

#[test]
fn quat_lattice_add_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_add");
    for v in &f.vectors {
        let a = read_lattice("a", &v.inputs);
        let b = read_lattice("b", &v.inputs);
        let exp = read_lattice("r", &v.outputs);
        let mut r = QuatLattice::new();
        quat_lattice_add(&mut r, &a, &b);
        assert!(ibz_eq(&r.denom, &exp.denom), "vector {}: denom", v.id);
        assert!(mat4x4_eq(&r.basis, &exp.basis), "vector {}: basis", v.id);
    }
}
