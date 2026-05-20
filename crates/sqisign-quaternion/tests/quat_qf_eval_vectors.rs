//! Differential test of `quat_qf_eval`.
mod common;
use common::{ibz_eq, read_ibz, read_mat4x4, read_vec4};
use sqisign_quaternion::{quat_qf_eval, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_qf_eval.json"
);

#[test]
fn quat_qf_eval_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_qf_eval");
    for v in &f.vectors {
        let qf = read_mat4x4("qf", &v.inputs);
        let c = read_vec4("c", &v.inputs);
        let exp = read_ibz("r", &v.outputs);
        let mut r = Ibz::zero();
        quat_qf_eval(&mut r, &qf, &c);
        assert!(ibz_eq(&r, &exp), "vector {}", v.id);
    }
}
