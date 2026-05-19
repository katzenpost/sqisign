//! Differential test of the ported `fp_neg` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` representation: five little-endian 8-byte limbs (the
//! reference's `digit_t = uint64_t` memory layout, `NWORDS_FIELD == 5`).
//!
//! The representation is redundant and non-canonical, so the assertion is
//! deliberately *bit-equality to the recorded reference output*, not a
//! congruence modulo `p`: the port must reproduce exactly what the
//! reference's `modneg` leaves in memory, including its unmasked limb 4.
//!
//! `modneg` is the unary analogue of `modsub` (`0 - b[i]` limbwise), so
//! the record shape is one input (`a`) and one output (`c`), unlike the
//! two-input `fp_sub`/`fp_add` vectors.

use sqisign_gf::{fp_neg, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/gf/fp_neg.json");

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
fn fp_neg_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_neg");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    // The all-zero edge pattern is emitted as the first record; pin that
    // modneg of the canonical zero is the bit-exact canonical zero, both
    // as the reference recorded it and as a count assertion (precedent
    // mp_neg/mp_mul).
    let mut zero_inputs = 0usize;
    let mut zero_from_nonzero = 0usize;

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = fp_from("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        let mut c = [0u64; NWORDS_FIELD];
        fp_neg(&mut c, &a);
        assert_eq!(
            c, expected,
            "vector {} diverged from the C reference at the fp_t boundary",
            v.id
        );

        let a_zero = a == [0u64; NWORDS_FIELD];
        let c_zero = c == [0u64; NWORDS_FIELD];
        if a_zero {
            zero_inputs += 1;
            assert_eq!(
                c, [0u64; NWORDS_FIELD],
                "fp_neg(0) must be the bit-exact canonical zero (vector {})",
                v.id
            );
        } else if c_zero {
            zero_from_nonzero += 1;
        }
    }

    // Exactly one all-zero input in the battery: the first edge pattern.
    // It is the only vector whose output is the all-zero representative
    // (no nonzero full-width input negates to bit-exact zero across the
    // 1012-vector battery). Pinning the counts guards against the edge
    // pattern silently dropping out of the emitter.
    assert_eq!(zero_inputs, 1, "expected exactly one all-zero input vector");
    assert_eq!(
        zero_from_nonzero, 0,
        "no nonzero input should negate to the bit-exact all-zero representative"
    );
}
