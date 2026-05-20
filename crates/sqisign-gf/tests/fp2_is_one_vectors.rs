//! Differential test of the ported `fp2_is_one` against the
//! committed C-derived vectors. Unary predicate shape:
//! (a_re, a_im) -> result (4-byte LE u32 mask).
//!
//! `fp2_is_one(a)` returns the all-ones mask iff
//! `a == (Montgomery_ONE, 0)`. The 12 x 12 edge battery contains the
//! `(Montgomery_ONE, zero)` pair (when, separately, the real half
//! happens to equal the precomputed ONE), but more commonly all 144
//! pattern pairs are *not* one; the seeded sweep almost never hits the
//! positive case. The bulk of records therefore test the negative
//! outcome; bit-equality on every recorded vector is the only check.

use sqisign_gf::{fp2_is_one, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_is_one.json"
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
fn fp2_is_one_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_is_one");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let a_re = fp_from(
            "a_re",
            &decode("a_re", &v.inputs["a_re"]).expect("a_re hex"),
        );
        let a_im = fp_from(
            "a_im",
            &decode("a_im", &v.inputs["a_im"]).expect("a_im hex"),
        );
        let r_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(r_bytes.len(), 4);
        let expected = u32::from_le_bytes(r_bytes.try_into().unwrap());

        let a = Fp2 { re: a_re, im: a_im };
        let got = fp2_is_one(&a);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference at the fp2_is_one result mask",
            v.id
        );
    }
}
