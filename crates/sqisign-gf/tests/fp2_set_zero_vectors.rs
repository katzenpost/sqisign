//! Differential test of the ported `fp2_set_zero` against the
//! committed C-derived vectors. The differential boundary is the raw
//! internal `fp2_t` representation: two `fp_t` halves, each five
//! little-endian 8-byte limbs (the reference's `digit_t = uint64_t`
//! memory layout, `NWORDS_FIELD == 5` per half).
//!
//! `fp2_set_zero(x)` zeroes both halves of `x`. The record carries
//! the destination's pre-call limbs as `prefill_re` / `prefill_im`
//! and the post-call limbs as `out_re` / `out_im`; a port that
//! quietly skipped either half would leave the corresponding
//! pre-fill bytes visible and diverge on at least one of the two
//! output limb vectors.

use sqisign_gf::{fp2_set_zero, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_set_zero.json"
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
fn fp2_set_zero_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_set_zero");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let pre_re = fp_from(
            "prefill_re",
            &decode("prefill_re", &v.inputs["prefill_re"]).expect("prefill_re hex"),
        );
        let pre_im = fp_from(
            "prefill_im",
            &decode("prefill_im", &v.inputs["prefill_im"]).expect("prefill_im hex"),
        );
        let exp_re = fp_from(
            "out_re",
            &decode("out_re", &v.outputs["out_re"]).expect("out_re hex"),
        );
        let exp_im = fp_from(
            "out_im",
            &decode("out_im", &v.outputs["out_im"]).expect("out_im hex"),
        );

        let mut x = Fp2 {
            re: pre_re,
            im: pre_im,
        };
        fp2_set_zero(&mut x);
        assert_eq!(
            x.re, exp_re,
            "vector {} diverged from the C reference on out_re",
            v.id
        );
        assert_eq!(
            x.im, exp_im,
            "vector {} diverged from the C reference on out_im",
            v.id
        );
    }
}
