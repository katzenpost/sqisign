//! Differential test of `ibz_vec_4_add`.
mod common;
use common::{read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_vec_4_add, ibz_vec_4_new};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_vec_4_add.json"
);

#[test]
fn ibz_vec_4_add_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_vec_4_add");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let a = read_vec4("a", &v.inputs);
        let b = read_vec4("b", &v.inputs);
        let exp = read_vec4("r", &v.outputs);
        let mut r = ibz_vec_4_new();
        ibz_vec_4_add(&mut r, &a, &b);
        assert!(vec4_eq(&r, &exp), "vector {}", v.id);
    }
}
