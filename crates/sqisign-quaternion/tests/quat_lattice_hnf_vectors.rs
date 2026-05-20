//! Differential test of `quat_lattice_hnf`.
mod common;
use common::{ibz_eq, mat4x4_eq, read_lattice};
use sqisign_quaternion::quat_lattice_hnf;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_hnf.json"
);

#[test]
fn quat_lattice_hnf_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_hnf");
    for v in &f.vectors {
        let mut l = read_lattice("l", &v.inputs);
        let exp = read_lattice("r", &v.outputs);
        quat_lattice_hnf(&mut l);
        assert!(ibz_eq(&l.denom, &exp.denom), "vector {}: denom", v.id);
        assert!(mat4x4_eq(&l.basis, &exp.basis), "vector {}: basis", v.id);
    }
}
