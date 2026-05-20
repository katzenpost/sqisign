//! Differential test of the ported `fp2_is_equal` against the
//! committed C-derived vectors. Binary predicate shape:
//! (a_re, a_im, b_re, b_im) -> result (4-byte LE u32 mask).
//!
//! `fp2_is_equal(a, b)` ANDs the two per-component `fp_is_equal` masks;
//! it is `0xFFFFFFFF` iff both halves are equal modulo `p` (the per-
//! half `redc` canonicalisation happens inside `fp_is_equal`). The
//! battery includes explicit non-canonical equality witnesses (zero
//! vs the radix-2^51 encoding of `p`) so the redc-first structure is
//! pinned.

use sqisign_gf::{fp2_is_equal, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_is_equal.json"
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
fn fp2_is_equal_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_is_equal");
    let total = file.vectors.len();
    assert!(total >= 1000);

    let mut pos = 0usize;
    let mut neg = 0usize;
    for v in &file.vectors {
        let a_re = fp_from(
            "a_re",
            &decode("a_re", &v.inputs["a_re"]).expect("a_re hex"),
        );
        let a_im = fp_from(
            "a_im",
            &decode("a_im", &v.inputs["a_im"]).expect("a_im hex"),
        );
        let b_re = fp_from(
            "b_re",
            &decode("b_re", &v.inputs["b_re"]).expect("b_re hex"),
        );
        let b_im = fp_from(
            "b_im",
            &decode("b_im", &v.inputs["b_im"]).expect("b_im hex"),
        );
        let r_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(r_bytes.len(), 4);
        let expected = u32::from_le_bytes(r_bytes.try_into().unwrap());

        let a = Fp2 { re: a_re, im: a_im };
        let b = Fp2 { re: b_re, im: b_im };
        let got = fp2_is_equal(&a, &b);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference at the fp2_is_equal result mask",
            v.id
        );
        match expected {
            0xFFFF_FFFF => pos += 1,
            0 => neg += 1,
            other => panic!("vector {} recorded non-mask result {other:#x}", v.id),
        }
    }
    assert_eq!(pos + neg, total);
    assert!(
        pos > 0 && neg > 0,
        "both outcomes must be exercised: pos={pos}, neg={neg}"
    );
}
