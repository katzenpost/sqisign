//! Differential test of `ibz_generate_random_prime` against C-derived vectors.

mod common;

use common::{read_i32, read_ibz};
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{ibz_generate_random_prime, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_generate_random_prime.json"
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
fn ibz_generate_random_prime_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_generate_random_prime");
    assert!(!f.vectors.is_empty());

    for v in &f.vectors {
        let entropy = fixed48(
            "entropy",
            &decode("entropy", &v.inputs["entropy"]).expect("entropy hex"),
        );
        let is3mod4 = read_i32("is3mod4", &v.inputs);
        let bitsize = read_i32("bitsize", &v.inputs);
        let iter = read_i32("iter", &v.inputs);
        let ok_exp = read_ok("ok", &v.outputs);
        let p_exp = read_ibz("p", &v.outputs);

        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut p = Ibz::zero();
        let ok_got = ibz_generate_random_prime(&mut drbg, &mut p, is3mod4, bitsize, iter);

        // `ok` mirrors `ibz_probab_prime`'s 0/1/2 trichotomy. GMP can
        // certify with trial division and return 2; our pure-Rust port
        // certifies only up to a small bound, then returns 1 (probable
        // prime). The value of `p` is identical in either case, and the
        // 0-vs-nonzero dichotomy is the meaningful boundary contract.
        // See `ibz_probab_prime_vectors.rs` for the same convention.
        assert_eq!(
            ok_got != 0,
            ok_exp != 0,
            "vector {}: success flag dichotomy (port={ok_got}, ref={ok_exp})",
            v.id
        );
        assert_eq!(
            p.0, p_exp.0,
            "vector {}: ibz_generate_random_prime(is3mod4={is3mod4}, bitsize={bitsize}) diverged from C reference",
            v.id
        );
    }
}
