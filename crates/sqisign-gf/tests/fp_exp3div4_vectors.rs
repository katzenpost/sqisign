//! Differential test of the ported `fp_exp3div4` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! The representation is redundant and non-canonical, so the assertion is
//! deliberately *bit-equality to the recorded reference output*, not a
//! congruence modulo `p`: the port must reproduce exactly what the
//! reference's `modpro` leaves in memory, including its full 64-bit
//! truncated (unmasked) limb 4. `fp_exp3div4` is the one-liner
//! `modpro(*a, *out)` (see
//! `the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:652..656`): the
//! progenitor, `a^((p-3)/4) mod p`, the building block the rest of the
//! chain (`fp_inv`, `fp_sqrt`, `fp_is_square`) folds further to obtain
//! inverses, square roots, and the Euler criterion. The unary record
//! shape (one input `a`, one output `c`) is the same as `fp_neg`/`fp_sqr`.

use sqisign_gf::{fp_exp3div4, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_exp3div4.json"
);

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
fn fp_exp3div4_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_exp3div4");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = fp_from("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        let mut c = [0u64; NWORDS_FIELD];
        fp_exp3div4(&mut c, &a);
        assert_eq!(
            c, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );
    }
}

#[test]
fn fp_exp3div4_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second call
    // on identical input yields identical limbs. modpro is a fixed
    // addition chain over modsqr/modmul/modnsqr, with no randomness, no
    // global state, and no hidden inputs; this assertion is cheap and
    // catches any future state leak.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let mut c1 = [0u64; NWORDS_FIELD];
        let mut c2 = [0u64; NWORDS_FIELD];
        fp_exp3div4(&mut c1, &a);
        fp_exp3div4(&mut c2, &a);
        assert_eq!(
            c1, c2,
            "fp_exp3div4 not deterministic on input from vector {}",
            v.id
        );
    }
}
