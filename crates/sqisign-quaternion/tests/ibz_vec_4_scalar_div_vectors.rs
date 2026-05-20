//! Differential test of `ibz_vec_4_scalar_div`.
mod common;
use common::{read_i32, read_ibz, read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_vec_4_new, ibz_vec_4_scalar_div};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_vec_4_scalar_div.json"
);

#[test]
fn ibz_vec_4_scalar_div_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_vec_4_scalar_div");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let s = read_ibz("s", &v.inputs);
        let vv = read_vec4("v", &v.inputs);
        let exp_r = read_vec4("r", &v.outputs);
        let exp_ok = read_i32("ok", &v.outputs);
        let mut r = ibz_vec_4_new();
        let ok = ibz_vec_4_scalar_div(&mut r, &s, &vv);
        assert_eq!(ok & 0xff, exp_ok & 0xff, "vector {}: ok", v.id);
        assert!(vec4_eq(&r, &exp_r), "vector {}: r", v.id);
    }
}
