//! Differential test of the ported `fp2_set_one` against the
//! committed C-derived vectors. Setter shape; same record layout as
//! `fp2_set_zero` (prefill_re, prefill_im -> out_re, out_im).
//!
//! `fp2_set_one(x)` writes the Montgomery representative of `1` into
//! the real half and zero into the imaginary half. The same
//! `MONTGOMERY_ONE` bit pattern `fp_set_one` produces appears in
//! `x.re`; this differential gate pins both halves bit-for-bit.

use sqisign_gf::{fp2_set_one, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_set_one.json"
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
fn fp2_set_one_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_set_one");
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
        fp2_set_one(&mut x);
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
