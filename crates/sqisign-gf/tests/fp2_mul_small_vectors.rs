//! Differential test of the ported `fp2_mul_small` against the
//! committed C-derived vectors. Field-and-value shape:
//! (a_re, a_im, val) -> (out_re, out_im).
//!
//! `fp2_mul_small(x, y, n)` applies the scalar `n` to BOTH halves via
//! `fp_mul_small` (unlike `fp2_set_small`, which only writes the real
//! half). The recorded `val` is 8 bytes for shape uniformity with
//! `fp2_set_small`; only the low 32 bits are consumed at the boundary.

use sqisign_gf::{fp2_mul_small, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_mul_small.json"
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
fn fp2_mul_small_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_mul_small");
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
        let val_bytes = decode("val", &v.inputs["val"]).expect("val hex");
        assert_eq!(val_bytes.len(), 8);
        let val_u64 = u64::from_le_bytes(val_bytes.try_into().unwrap());
        let exp_re = fp_from(
            "out_re",
            &decode("out_re", &v.outputs["out_re"]).expect("out_re hex"),
        );
        let exp_im = fp_from(
            "out_im",
            &decode("out_im", &v.outputs["out_im"]).expect("out_im hex"),
        );

        let a = Fp2 { re: a_re, im: a_im };
        let mut x = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        fp2_mul_small(&mut x, &a, val_u64 as u32);
        assert_eq!(x.re, exp_re, "vector {} diverged on out_re", v.id);
        assert_eq!(x.im, exp_im, "vector {} diverged on out_im", v.id);
    }
}
