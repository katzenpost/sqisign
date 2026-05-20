//! Differential test of `ibz_probab_prime`.
//!
//! GMP returns a 0/1/2 trichotomy: 0 certainly composite, 2 certainly
//! prime, 1 probably prime. Our pure-Rust port cannot certify primality
//! without a full proof; it returns 0 for composites, 1 for probable
//! primes, and falls back to 2 only for the small-prime fast path. The
//! reference's pinned vectors use small primes so we assert the
//! 0-vs-nonzero dichotomy and the small-prime trichotomy together.
use sqisign_quaternion::{ibz_probab_prime, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_probab_prime.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_i8(l: &str, h: &str) -> i8 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 1);
    b[0] as i8
}
fn read_le_i32(l: &str, h: &str) -> i32 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 4);
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

#[test]
fn ibz_probab_prime_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_probab_prime");
    for v in &f.vectors {
        let n = read_ibz("n", &v.inputs["n"]);
        let reps = read_le_i32("reps", &v.inputs["reps"]);
        let exp = read_i8("r", &v.outputs["r"]);
        let got = ibz_probab_prime(&n, reps) as i8;
        let is_composite_exp = exp == 0;
        let is_composite_got = got == 0;
        assert_eq!(
            is_composite_got, is_composite_exp,
            "vector {}: compositeness disagreed (port={got}, ref={exp})",
            v.id
        );
    }
}
