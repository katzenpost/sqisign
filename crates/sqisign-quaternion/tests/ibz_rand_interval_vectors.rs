//! Differential test of `ibz_rand_interval` against C-derived vectors.
//!
//! Each vector seeds a fresh KAT-only [`CtrDrbg`] with a recorded
//! 48-byte entropy block (no personalization), calls
//! `ibz_rand_interval(&mut drbg, rand, a, b)`, and asserts the
//! canonical-bytes representation of `rand` matches the C reference's
//! output. The recorded `ok` flag must also agree (1 on success; the
//! reference's only failure path is a `randombytes` error, which the
//! [`RngSource`] contract treats as a panic).

mod common;

use common::{read_i32, read_ibz};
use sqisign_common::CtrDrbg;
use sqisign_quaternion::{ibz_rand_interval, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_rand_interval.json"
);

fn fixed48(label: &str, v: &[u8]) -> [u8; 48] {
    let mut a = [0u8; 48];
    assert_eq!(v.len(), 48, "{label} must be 48 bytes, got {}", v.len());
    a.copy_from_slice(v);
    a
}

#[test]
fn ibz_rand_interval_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_rand_interval");
    assert!(!f.vectors.is_empty(), "vector battery empty");

    for v in &f.vectors {
        let entropy = fixed48(
            "entropy",
            &decode("entropy", &v.inputs["entropy"]).expect("entropy hex"),
        );
        let a = read_ibz("a", &v.inputs);
        let b = read_ibz("b", &v.inputs);

        let ok_exp = read_i32_or_byte("ok", &v.outputs);
        let r_exp = read_ibz("r", &v.outputs);

        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut r = Ibz::zero();
        let ok_got = ibz_rand_interval(&mut drbg, &mut r, &a, &b);

        assert_eq!(ok_got, ok_exp, "vector {}: ok flag", v.id);
        assert_eq!(
            r.0, r_exp.0,
            "vector {}: ibz_rand_interval(a, b) diverged from C reference",
            v.id
        );
    }
}

/// The `ok` field in the recorded vectors is a single byte. Decode it
/// as a signed int8 widened to i32 so this test mirrors the convention
/// `ibz_rand_*` returns.
fn read_i32_or_byte(name: &str, m: &std::collections::BTreeMap<String, String>) -> i32 {
    let b = decode(name, &m[name]).expect("hex decode");
    match b.len() {
        1 => b[0] as i8 as i32,
        4 => read_i32(name, m),
        n => panic!("{name}: unexpected byte length {n}"),
    }
}
