//! Differential test of `ibz_mat_4x4_transpose`.
mod common;
use common::{mat4x4_eq, read_mat4x4};
use sqisign_quaternion::{ibz_mat_4x4_new, ibz_mat_4x4_transpose};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_4x4_transpose.json"
);

#[test]
fn ibz_mat_4x4_transpose_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_4x4_transpose");
    for v in &f.vectors {
        let a = read_mat4x4("a", &v.inputs);
        let exp = read_mat4x4("r", &v.outputs);
        let mut r = ibz_mat_4x4_new();
        ibz_mat_4x4_transpose(&mut r, &a);
        assert!(mat4x4_eq(&r, &exp), "vector {}", v.id);
    }
}
