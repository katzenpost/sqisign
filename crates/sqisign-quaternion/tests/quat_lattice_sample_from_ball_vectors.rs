//! Differential test of `quat_lattice_sample_from_ball` against
//! C-derived vectors.
//!
//! Each vector seeds a fresh KAT-only [`CtrDrbg`] with a recorded 48-byte
//! entropy block (no personalization), reads the input lattice (denom +
//! 4x4 basis), prime `p`, and radius, then calls
//! `quat_lattice_sample_from_ball(&mut drbg, &mut res, &lat, &alg, &radius)`.
//! The recorded `ok` flag and the produced `res` (denominator + 4
//! coordinates) must agree with the C reference byte-for-byte.

mod common;

use common::{ibz_eq, read_ibz, read_lattice};
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{quat_lattice_sample_from_ball, Ibz, QuatAlg, QuatAlgElem};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/quat_lattice_sample_from_ball.json"
);

fn fixed48(v: &[u8]) -> [u8; 48] {
    let mut a = [0u8; 48];
    assert_eq!(v.len(), 48, "entropy must be 48 bytes, got {}", v.len());
    a.copy_from_slice(v);
    a
}

fn read_ok(name: &str, m: &std::collections::BTreeMap<String, String>) -> i32 {
    let b = decode(name, &m[name]).expect("hex decode");
    assert_eq!(b.len(), 1, "{name}: expected 1-byte ok flag");
    b[0] as i8 as i32
}

#[test]
fn quat_lattice_sample_from_ball_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(f.boundary, "sqisign_quaternion::quat_lattice_sample_from_ball");
    assert!(!f.vectors.is_empty(), "vector battery empty");

    for v in &f.vectors {
        let entropy = fixed48(&decode("entropy", &v.inputs["entropy"]).expect("entropy hex"));
        let lat = read_lattice("lat", &v.inputs);
        let p = read_ibz("p", &v.inputs);
        let radius = read_ibz("radius", &v.inputs);

        let ok_exp = read_ok("ok", &v.outputs);
        let res_denom_exp = read_ibz("res_denom", &v.outputs);
        let res_c_exp: [Ibz; 4] = [
            read_ibz("res_c0", &v.outputs),
            read_ibz("res_c1", &v.outputs),
            read_ibz("res_c2", &v.outputs),
            read_ibz("res_c3", &v.outputs),
        ];

        let alg = QuatAlg::init_set(&p);
        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut res = QuatAlgElem::new();
        let ok_got = quat_lattice_sample_from_ball(&mut drbg, &mut res, &lat, &alg, &radius);

        assert_eq!(ok_got, ok_exp, "vector {}: ok flag", v.id);
        assert!(
            ibz_eq(&res.denom, &res_denom_exp),
            "vector {}: res.denom diverged (got {}, want {})",
            v.id,
            res.denom.0,
            res_denom_exp.0
        );
        for k in 0..4 {
            assert!(
                ibz_eq(&res.coord[k], &res_c_exp[k]),
                "vector {}: res.coord[{}] diverged (got {}, want {})",
                v.id,
                k,
                res.coord[k].0,
                res_c_exp[k].0
            );
        }
    }
}
