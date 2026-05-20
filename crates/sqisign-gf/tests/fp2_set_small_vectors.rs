//! Differential test of the ported `fp2_set_small` against the
//! committed C-derived vectors. Setter-with-a-value shape:
//! (prefill_re, prefill_im, val) -> (out_re, out_im).
//!
//! `fp2_set_small(x, val)` writes the int32-narrowed `val` into the
//! real half via `fp_set_small` and zeroes the imaginary half. Only
//! `out_re` is affected by `val`; `out_im` is always the canonical
//! zero regardless of `prefill_im`. The pre-fills exercise the
//! complete-overwrite property on both halves.

use sqisign_gf::{fp2_set_small, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_set_small.json"
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
fn fp2_set_small_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_set_small");
    assert!(file.vectors.len() >= 1000);

    for v in &file.vectors {
        let pre_re = fp_from(
            "prefill_re",
            &decode("prefill_re", &v.inputs["prefill_re"]).expect("prefill_re hex"),
        );
        let pre_im = fp_from(
            "prefill_im",
            &decode("prefill_im", &v.inputs["prefill_im"]).expect("prefill_im hex"),
        );
        let val_bytes = decode("val", &v.inputs["val"]).expect("val hex");
        assert_eq!(val_bytes.len(), 8);
        let val = u64::from_le_bytes(val_bytes.try_into().unwrap());
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
        fp2_set_small(&mut x, val);
        assert_eq!(x.re, exp_re, "vector {} diverged on out_re", v.id);
        assert_eq!(x.im, exp_im, "vector {} diverged on out_im", v.id);
    }
}
