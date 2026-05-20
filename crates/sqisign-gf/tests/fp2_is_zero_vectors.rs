//! Differential test of the ported `fp2_is_zero` against the
//! committed C-derived vectors. Unary predicate shape:
//! (a_re, a_im) -> result (4-byte LE u32 mask, `0xFFFFFFFF` on
//! positive outcome, `0` on negative).

use sqisign_gf::{fp2_is_zero, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_is_zero.json"
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
fn fp2_is_zero_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_is_zero");
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
        let r_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(r_bytes.len(), 4);
        let expected = u32::from_le_bytes(r_bytes.try_into().unwrap());

        let a = Fp2 { re: a_re, im: a_im };
        let got = fp2_is_zero(&a);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference at the fp2_is_zero result mask",
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
