//! Differential test of the ported `fp2_sqrt` against the committed
//! C-derived vectors. Unary shape: (a_re, a_im) -> (c_re, c_im).
//!
//! `fp2_sqrt(x)` is in-place; the cdump harness uses the
//! `fp2_sqrt_adapter` to snapshot the input as `a` and the post-call
//! buffer as `c`. The replay copies the recorded inputs into a fresh
//! `Fp2` and reads the post-call limbs. The sqrt algorithm is the
//! Aardal et al construction (eprint 2024/1563) the reference
//! transcribes in `fp2.c:148..202`; bit-equality on every recorded
//! vector (squares and non-squares alike) is the pin.

use sqisign_gf::{fp2_sqrt, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_sqrt.json"
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
fn fp2_sqrt_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_sqrt");
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
        let exp_re = fp_from(
            "c_re",
            &decode("c_re", &v.outputs["c_re"]).expect("c_re hex"),
        );
        let exp_im = fp_from(
            "c_im",
            &decode("c_im", &v.outputs["c_im"]).expect("c_im hex"),
        );

        let mut x = Fp2 { re: a_re, im: a_im };
        fp2_sqrt(&mut x);
        assert_eq!(x.re, exp_re, "vector {} diverged on c_re", v.id);
        assert_eq!(x.im, exp_im, "vector {} diverged on c_im", v.id);
    }
}
