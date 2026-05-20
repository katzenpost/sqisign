//! Differential test of `ibz_mat_4x4_eval`.
mod common;
use common::{read_mat4x4, read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_mat_4x4_eval, ibz_vec_4_new};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_4x4_eval.json"
);

#[test]
fn ibz_mat_4x4_eval_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_4x4_eval");
    for v in &f.vectors {
        let m = read_mat4x4("m", &v.inputs);
        let vv = read_vec4("v", &v.inputs);
        let exp = read_vec4("r", &v.outputs);
        let mut r = ibz_vec_4_new();
        ibz_mat_4x4_eval(&mut r, &m, &vv);
        assert!(vec4_eq(&r, &exp), "vector {}", v.id);
    }
}
