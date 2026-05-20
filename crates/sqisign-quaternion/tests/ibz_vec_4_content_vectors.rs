//! Differential test of `ibz_vec_4_content`.
mod common;
use common::{ibz_eq, read_ibz, read_vec4};
use sqisign_quaternion::{ibz_vec_4_content, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_vec_4_content.json"
);

#[test]
fn ibz_vec_4_content_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_vec_4_content");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let vv = read_vec4("v", &v.inputs);
        let exp = read_ibz("c", &v.outputs);
        let mut c = Ibz::zero();
        ibz_vec_4_content(&mut c, &vv);
        assert!(ibz_eq(&c, &exp), "vector {}", v.id);
    }
}
