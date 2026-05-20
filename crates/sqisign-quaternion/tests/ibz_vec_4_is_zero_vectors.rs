//! Differential test of `ibz_vec_4_is_zero`.
mod common;
use common::{read_i32, read_vec4};
use sqisign_quaternion::ibz_vec_4_is_zero;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_vec_4_is_zero.json"
);

#[test]
fn ibz_vec_4_is_zero_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_vec_4_is_zero");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let x = read_vec4("x", &v.inputs);
        let exp = read_i32("res", &v.outputs);
        let r = ibz_vec_4_is_zero(&x);
        assert_eq!(r & 0xff, exp & 0xff, "vector {}", v.id);
    }
}
