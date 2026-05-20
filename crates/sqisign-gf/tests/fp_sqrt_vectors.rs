//! Differential test of the ported `fp_sqrt` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! `fp_sqrt` is the *in-place* one-liner `modsqrt(*a, NULL, *a)` (see
//! `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:640..644`). The
//! cdump harness snapshots the input into a separate `a` and the
//! in-place result into `c`, so the recorded record carries both sides
//! of the in-place mutation; the port wrapper snapshots `*a` into a
//! local to resolve the borrow conflict and then writes back, exactly
//! the same way `fp_inv` does.
//!
//! For `p == 3 mod 4` (the case for `p5248`), `sqrt(x) == x^((p+1)/4) ==
//! progenitor(x) * x mod p`. On a non-residue input the reference makes
//! no defensive check; whatever `progenitor(x) * x` evaluates to is
//! returned as "garbage but deterministic" output. The differential pin
//! is bit-equality to *that* recorded output regardless of QR status:
//! the boundary is the byte-for-byte behaviour of the reference, not a
//! mathematically-justified-only-on-residues claim.

use sqisign_gf::{fp_sqrt, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/gf/fp_sqrt.json");

fn fp_from(label: &str, bytes: &[u8]) -> Fp {
    assert_eq!(
        bytes.len(),
        NWORDS_FIELD * 8,
        "{label} must be exactly {NWORDS_FIELD} u64 limbs"
    );
    let mut limbs = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks_exact(8).enumerate() {
        limbs[i] = u64::from_le_bytes(chunk.try_into().unwrap());
    }
    limbs
}

#[test]
fn fp_sqrt_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_sqrt");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = fp_from("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        let mut x = a;
        fp_sqrt(&mut x);
        assert_eq!(
            x, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );
    }
}

#[test]
fn fp_sqrt_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical limbs. modsqrt is a
    // fixed addition chain over modpro/modmul with no randomness, no
    // global state, and no hidden inputs; even on a non-residue input
    // the "garbage" output is deterministic for a given input.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let mut x1 = a;
        let mut x2 = a;
        fp_sqrt(&mut x1);
        fp_sqrt(&mut x2);
        assert_eq!(
            x1, x2,
            "fp_sqrt not deterministic on input from vector {}",
            v.id
        );
    }
}
