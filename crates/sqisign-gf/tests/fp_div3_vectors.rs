//! Differential test of the ported `fp_div3` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! The representation is redundant and non-canonical, so the assertion is
//! deliberately *bit-equality to the recorded reference output*, not a
//! congruence modulo `p`: the port must reproduce exactly what the
//! reference's `modmul(THREE_INV, a, out)` leaves in memory, including
//! its full 64-bit truncated (unmasked) limb 4.
//!
//! `fp_div3` is the one-liner `modmul(THREE_INV, *a, *out)` (see
//! `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:658..662`), so the
//! record shape is one input (`a`) and one output (`c`), the unary
//! cadence shared with `fp_neg`, `fp_sqr`, and `fp_half`.

use sqisign_gf::{fp_div3, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/gf/fp_div3.json");

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
fn fp_div3_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_div3");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = fp_from("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        let mut c = [0u64; NWORDS_FIELD];
        fp_div3(&mut c, &a);
        assert_eq!(
            c, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );
    }
}

#[test]
fn fp_div3_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second call
    // on identical inputs yields identical limbs. fp_div3 is the
    // one-liner modmul(THREE_INV, a, out); modmul has no randomness, no
    // global state, and no hidden inputs, and THREE_INV is a fixed const,
    // so the assertion is cheap and catches any future state leak.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let mut c1 = [0u64; NWORDS_FIELD];
        let mut c2 = [0u64; NWORDS_FIELD];
        fp_div3(&mut c1, &a);
        fp_div3(&mut c2, &a);
        assert_eq!(
            c1, c2,
            "fp_div3 not deterministic on input from vector {}",
            v.id
        );
    }
}
