//! Differential test of `quat_alg_coord_mul`.
use sqisign_quaternion::{quat_alg_coord_mul, Ibz, QuatAlg};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_alg_coord_mul.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}

#[test]
fn quat_alg_coord_mul_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::quat_alg_coord_mul");
    assert!(f.vectors.len() >= 100);
    for v in &f.vectors {
        let a = [
            read_ibz("a_c0", &v.inputs["a_c0"]),
            read_ibz("a_c1", &v.inputs["a_c1"]),
            read_ibz("a_c2", &v.inputs["a_c2"]),
            read_ibz("a_c3", &v.inputs["a_c3"]),
        ];
        let b = [
            read_ibz("b_c0", &v.inputs["b_c0"]),
            read_ibz("b_c1", &v.inputs["b_c1"]),
            read_ibz("b_c2", &v.inputs["b_c2"]),
            read_ibz("b_c3", &v.inputs["b_c3"]),
        ];
        let p = read_ibz("p", &v.inputs["p"]);
        let exp = [
            read_ibz("r_c0", &v.outputs["r_c0"]),
            read_ibz("r_c1", &v.outputs["r_c1"]),
            read_ibz("r_c2", &v.outputs["r_c2"]),
            read_ibz("r_c3", &v.outputs["r_c3"]),
        ];
        let alg = QuatAlg::init_set(&p);
        let mut r = [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()];
        quat_alg_coord_mul(&mut r, &a, &b, &alg);
        for k in 0..4 {
            assert_eq!(r[k].0, exp[k].0, "vector {}: coord {} mismatch", v.id, k);
        }
    }
}
