//! Differential test of the ported `fp2_sqrt_verify` against the
//! committed C-derived vectors. Combined unary + result shape:
//! (a_in_re, a_in_im) -> (a_out_re, a_out_im, result).
//!
//! `fp2_sqrt_verify(a)` snapshots the input, calls `fp2_sqrt(a)` in
//! place, squares the result, and returns `fp2_is_equal(snap, square)`.
//! The post-call `a` is the sqrt routine's output; `result` is
//! `0xFFFFFFFF` iff the original input was a square in `Fp^2`.

use sqisign_gf::{fp2_sqrt_verify, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_sqrt_verify.json"
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
fn fp2_sqrt_verify_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_sqrt_verify");
    let total = file.vectors.len();
    assert!(total >= 1000);

    let mut pos = 0usize;
    let mut neg = 0usize;
    for v in &file.vectors {
        let a_re = fp_from(
            "a_in_re",
            &decode("a_in_re", &v.inputs["a_in_re"]).expect("a_in_re hex"),
        );
        let a_im = fp_from(
            "a_in_im",
            &decode("a_in_im", &v.inputs["a_in_im"]).expect("a_in_im hex"),
        );
        let exp_re = fp_from(
            "a_out_re",
            &decode("a_out_re", &v.outputs["a_out_re"]).expect("a_out_re hex"),
        );
        let exp_im = fp_from(
            "a_out_im",
            &decode("a_out_im", &v.outputs["a_out_im"]).expect("a_out_im hex"),
        );
        let r_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(r_bytes.len(), 4);
        let expected_result = u32::from_le_bytes(r_bytes.try_into().unwrap());

        let mut x = Fp2 { re: a_re, im: a_im };
        let got = fp2_sqrt_verify(&mut x);
        assert_eq!(
            got, expected_result,
            "vector {} diverged from the C reference on result",
            v.id
        );
        assert_eq!(x.re, exp_re, "vector {} diverged on a_out_re", v.id);
        assert_eq!(x.im, exp_im, "vector {} diverged on a_out_im", v.id);

        match expected_result {
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
