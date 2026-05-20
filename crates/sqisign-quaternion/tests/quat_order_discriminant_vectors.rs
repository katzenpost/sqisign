//! Differential test of `quat_order_discriminant`.
mod common;
use common::{ibz_eq, read_i32, read_ibz, read_lattice};
use sqisign_quaternion::{quat_order_discriminant, Ibz, QuatAlg};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_order_discriminant.json"
);

#[test]
fn quat_order_discriminant_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_order_discriminant");
    for v in &f.vectors {
        let o = read_lattice("o", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let exp_d = read_ibz("disc", &v.outputs);
        let exp_ok = read_i32("ok", &v.outputs);
        let alg = QuatAlg::init_set(&p);
        let mut disc = Ibz::zero();
        let ok = quat_order_discriminant(&mut disc, &o, &alg);
        assert_eq!(ok & 0xff, exp_ok & 0xff, "vector {}: ok", v.id);
        assert!(ibz_eq(&disc, &exp_d), "vector {}: disc", v.id);
    }
}
