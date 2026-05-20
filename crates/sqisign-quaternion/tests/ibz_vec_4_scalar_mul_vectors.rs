//! Differential test of `ibz_vec_4_scalar_mul`.
mod common;
use common::{read_ibz, read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_vec_4_new, ibz_vec_4_scalar_mul};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_vec_4_scalar_mul.json"
);

#[test]
fn ibz_vec_4_scalar_mul_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_vec_4_scalar_mul");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let s = read_ibz("s", &v.inputs);
        let vv = read_vec4("v", &v.inputs);
        let exp = read_vec4("r", &v.outputs);
        let mut r = ibz_vec_4_new();
        ibz_vec_4_scalar_mul(&mut r, &s, &vv);
        assert!(vec4_eq(&r, &exp), "vector {}", v.id);
    }
}
