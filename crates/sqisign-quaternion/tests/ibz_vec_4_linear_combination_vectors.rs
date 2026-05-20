//! Differential test of `ibz_vec_4_linear_combination`.
mod common;
use common::{read_ibz, read_vec4, vec4_eq};
use sqisign_quaternion::{ibz_vec_4_linear_combination, ibz_vec_4_new};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_vec_4_linear_combination.json"
);

#[test]
fn ibz_vec_4_linear_combination_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::ibz_vec_4_linear_combination"
    );
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let ca = read_ibz("ca", &v.inputs);
        let va = read_vec4("va", &v.inputs);
        let cb = read_ibz("cb", &v.inputs);
        let vb = read_vec4("vb", &v.inputs);
        let exp = read_vec4("r", &v.outputs);
        let mut r = ibz_vec_4_new();
        ibz_vec_4_linear_combination(&mut r, &ca, &va, &cb, &vb);
        assert!(vec4_eq(&r, &exp), "vector {}", v.id);
    }
}
