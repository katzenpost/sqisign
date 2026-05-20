//! Differential test of `ec_biscalar_mul_ibz_vec`.

mod common;

use sqisign_ec::{ec_point_init, EcPoint};
use sqisign_id2iso::ec_biscalar_mul_ibz_vec;
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/id2iso/ec_biscalar_mul_ibz_vec.json"
);

#[test]
fn ec_biscalar_mul_ibz_vec_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_id2iso::ec_biscalar_mul_ibz_vec");
    assert!(!f.vectors.is_empty());
    for v in &f.vectors {
        let sc = common::read_vec2("sc", &v.inputs);
        let pq = common::ec_basis_from("pq", &v.inputs);
        let e = common::ec_curve_from("e", &v.inputs);
        let fparam = common::read_i32("f", &v.inputs);
        let exp = common::ec_point_from("r", &v.outputs);

        let mut r = EcPoint::zero();
        ec_point_init(&mut r);
        ec_biscalar_mul_ibz_vec(&mut r, &sc, fparam, &pq, &e);
        assert_eq!(r, exp, "vector {} diverged", v.id);
    }
}
