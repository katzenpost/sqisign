//! Differential test of `ibz_mat_4x4_inv_with_det_as_denom`.
mod common;
use common::{ibz_eq, mat4x4_eq, read_i32, read_ibz, read_mat4x4};
use sqisign_quaternion::{ibz_mat_4x4_inv_with_det_as_denom, ibz_mat_4x4_new, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_4x4_inv_with_det_as_denom.json"
);

#[test]
fn ibz_mat_4x4_inv_with_det_as_denom_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(
        f.boundary,
        "sqisign_quaternion::ibz_mat_4x4_inv_with_det_as_denom"
    );
    for v in &f.vectors {
        let m = read_mat4x4("m", &v.inputs);
        let exp_inv = read_mat4x4("inv", &v.outputs);
        let exp_det = read_ibz("det", &v.outputs);
        let exp_ok = read_i32("ok", &v.outputs);
        let mut inv = ibz_mat_4x4_new();
        let mut det = Ibz::zero();
        let ok = ibz_mat_4x4_inv_with_det_as_denom(Some(&mut inv), Some(&mut det), &m);
        assert_eq!(ok & 0xff, exp_ok & 0xff, "vector {}: ok", v.id);
        assert!(ibz_eq(&det, &exp_det), "vector {}: det", v.id);
        assert!(mat4x4_eq(&inv, &exp_inv), "vector {}: inv", v.id);
    }
}
