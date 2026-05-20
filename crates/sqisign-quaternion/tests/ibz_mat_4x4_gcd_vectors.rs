//! Differential test of `ibz_mat_4x4_gcd`.
mod common;
use common::{ibz_eq, read_ibz, read_mat4x4};
use sqisign_quaternion::{ibz_mat_4x4_gcd, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_4x4_gcd.json"
);

#[test]
fn ibz_mat_4x4_gcd_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_4x4_gcd");
    for v in &f.vectors {
        let a = read_mat4x4("a", &v.inputs);
        let exp = read_ibz("g", &v.outputs);
        let mut g = Ibz::zero();
        ibz_mat_4x4_gcd(&mut g, &a);
        assert!(ibz_eq(&g, &exp), "vector {}", v.id);
    }
}
