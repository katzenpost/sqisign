//! Differential test of `quat_lattice_equal`.
mod common;
use common::{read_i32, read_lattice};
use sqisign_quaternion::quat_lattice_equal;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_equal.json"
);

#[test]
fn quat_lattice_equal_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_equal");
    for v in &f.vectors {
        let a = read_lattice("a", &v.inputs);
        let b = read_lattice("b", &v.inputs);
        let exp = read_i32("res", &v.outputs);
        let r = quat_lattice_equal(&a, &b);
        assert_eq!(r & 0xff, exp & 0xff, "vector {}", v.id);
    }
}
