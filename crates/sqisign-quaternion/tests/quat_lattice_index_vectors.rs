//! Differential test of `quat_lattice_index`.
mod common;
use common::{ibz_eq, read_ibz, read_lattice};
use sqisign_quaternion::{quat_lattice_index, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_index.json"
);

#[test]
fn quat_lattice_index_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_index");
    for v in &f.vectors {
        let sub_ = read_lattice("sub", &v.inputs);
        let over = read_lattice("over", &v.inputs);
        let exp = read_ibz("idx", &v.outputs);
        let mut idx = Ibz::zero();
        quat_lattice_index(&mut idx, &sub_, &over);
        assert!(ibz_eq(&idx, &exp), "vector {}", v.id);
    }
}
