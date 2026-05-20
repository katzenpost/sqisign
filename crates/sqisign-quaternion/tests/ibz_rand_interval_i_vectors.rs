//! Differential test of `ibz_rand_interval_i` against C-derived vectors.

mod common;

use common::{read_i32, read_ibz};
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{ibz_rand_interval_i, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_rand_interval_i.json"
);

fn fixed48(label: &str, v: &[u8]) -> [u8; 48] {
    let mut a = [0u8; 48];
    assert_eq!(v.len(), 48, "{label} must be 48 bytes, got {}", v.len());
    a.copy_from_slice(v);
    a
}

fn read_ok(name: &str, m: &std::collections::BTreeMap<String, String>) -> i32 {
    let b = decode(name, &m[name]).expect("hex decode");
    assert_eq!(b.len(), 1);
    b[0] as i8 as i32
}

#[test]
fn ibz_rand_interval_i_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_rand_interval_i");
    assert!(!f.vectors.is_empty());

    for v in &f.vectors {
        let entropy = fixed48(
            "entropy",
            &decode("entropy", &v.inputs["entropy"]).expect("entropy hex"),
        );
        let a = read_i32("a", &v.inputs);
        let b = read_i32("b", &v.inputs);
        let ok_exp = read_ok("ok", &v.outputs);
        let r_exp = read_ibz("r", &v.outputs);

        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut r = Ibz::zero();
        let ok_got = ibz_rand_interval_i(&mut drbg, &mut r, a, b);

        assert_eq!(ok_got, ok_exp, "vector {}: ok flag", v.id);
        assert_eq!(
            r.0, r_exp.0,
            "vector {}: ibz_rand_interval_i({a}, {b}) diverged from C reference",
            v.id
        );
    }
}
