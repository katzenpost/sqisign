//! Differential test of the ported `fp2_decode` against the committed
//! C-derived vectors. Record shape: src (64 bytes) ->
//! (d_re, d_im, result).
//!
//! `fp2_decode(d, src)` reads 64 bytes, calls `fp_decode` on each
//! 32-byte half, and returns `re_mask & im_mask`. On any non-canonical
//! half the corresponding per-half d limbs are zeroed by `fp_decode`'s
//! `& res` mask; the combined `result` is `0` whenever either half is
//! out of range.

use sqisign_gf::{fp2_decode, Fp, Fp2, FP2_ENCODED_BYTES, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_decode.json"
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
fn fp2_decode_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_decode");
    let total = file.vectors.len();
    assert!(total >= 1000);

    let mut canonical = 0usize;
    let mut non_canonical = 0usize;

    for v in &file.vectors {
        let src_bytes = decode("src", &v.inputs["src"]).expect("src hex");
        assert_eq!(src_bytes.len(), FP2_ENCODED_BYTES);
        let mut src = [0u8; FP2_ENCODED_BYTES];
        src.copy_from_slice(&src_bytes);

        let exp_re = fp_from(
            "d_re",
            &decode("d_re", &v.outputs["d_re"]).expect("d_re hex"),
        );
        let exp_im = fp_from(
            "d_im",
            &decode("d_im", &v.outputs["d_im"]).expect("d_im hex"),
        );
        let r_bytes = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(r_bytes.len(), 4);
        let expected_result = u32::from_le_bytes(r_bytes.try_into().unwrap());

        let mut d = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        let got = fp2_decode(&mut d, &src);
        assert_eq!(
            got, expected_result,
            "vector {} diverged from the C reference on the result mask",
            v.id
        );
        assert_eq!(d.re, exp_re, "vector {} diverged on d_re", v.id);
        assert_eq!(d.im, exp_im, "vector {} diverged on d_im", v.id);

        match expected_result {
            0xFFFF_FFFF => canonical += 1,
            0 => non_canonical += 1,
            other => panic!(
                "vector {}: recorded mask {other:#x} is neither 0 nor 0xFFFFFFFF",
                v.id
            ),
        }
    }
    assert_eq!(canonical + non_canonical, total);
    assert!(
        canonical > 0 && non_canonical > 0,
        "both classes must be exercised: canonical={canonical}, non_canonical={non_canonical}"
    );
}
