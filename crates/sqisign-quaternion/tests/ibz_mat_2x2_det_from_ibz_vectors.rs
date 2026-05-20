//! Differential test of `ibz_mat_2x2_det_from_ibz`.
mod common;
use common::{ibz_eq, read_ibz};
use sqisign_quaternion::{ibz_mat_2x2_det_from_ibz, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mat_2x2_det_from_ibz.json"
);

#[test]
fn ibz_mat_2x2_det_from_ibz_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mat_2x2_det_from_ibz");
    for v in &f.vectors {
        let a11 = read_ibz("a11", &v.inputs);
        let a12 = read_ibz("a12", &v.inputs);
        let a21 = read_ibz("a21", &v.inputs);
        let a22 = read_ibz("a22", &v.inputs);
        let exp = read_ibz("det", &v.outputs);
        let mut d = Ibz::zero();
        ibz_mat_2x2_det_from_ibz(&mut d, &a11, &a12, &a21, &a22);
        assert!(ibz_eq(&d, &exp), "vector {}", v.id);
    }
}
