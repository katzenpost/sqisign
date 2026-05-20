//! Differential test of `ibz_mat_4x4_is_hnf`.
mod common;
use common::{read_i32, read_mat4x4};
use sqisign_quaternion::ibz_mat_4x4_is_hnf;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_4x4_is_hnf.json"
);

#[test]
fn ibz_mat_4x4_is_hnf_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_4x4_is_hnf");
    for v in &f.vectors {
        let m = read_mat4x4("m", &v.inputs);
        let exp = read_i32("res", &v.outputs);
        let r = ibz_mat_4x4_is_hnf(&m);
        assert_eq!(r & 0xff, exp & 0xff, "vector {}", v.id);
    }
}
