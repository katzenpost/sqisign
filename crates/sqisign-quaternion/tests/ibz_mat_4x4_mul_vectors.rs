//! Differential test of `ibz_mat_4x4_mul`.
mod common;
use common::{mat4x4_eq, read_mat4x4};
use sqisign_quaternion::{ibz_mat_4x4_mul, ibz_mat_4x4_new};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_4x4_mul.json"
);

#[test]
fn ibz_mat_4x4_mul_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_4x4_mul");
    for v in &f.vectors {
        let a = read_mat4x4("a", &v.inputs);
        let b = read_mat4x4("b", &v.inputs);
        let exp = read_mat4x4("r", &v.outputs);
        let mut r = ibz_mat_4x4_new();
        ibz_mat_4x4_mul(&mut r, &a, &b);
        assert!(mat4x4_eq(&r, &exp), "vector {}", v.id);
    }
}
