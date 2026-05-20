//! Differential test of the ported `fp2_sqr` against the committed
//! C-derived vectors. Unary shape: (a_re, a_im) -> (c_re, c_im).
//! `fp2_sqr(c, a)` uses the difference-of-squares identity (sum, diff,
//! cross product, doubling) the reference transcribes in `fp2.c`.

use sqisign_gf::{fp2_sqr, Fp, Fp2, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/gf/fp2_sqr.json");

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
fn fp2_sqr_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_sqr");
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

        let a = Fp2 { re: a_re, im: a_im };
        let mut c = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        fp2_sqr(&mut c, &a);
        assert_eq!(c.re, exp_re, "vector {} diverged on c_re", v.id);
        assert_eq!(c.im, exp_im, "vector {} diverged on c_im", v.id);
    }
}
