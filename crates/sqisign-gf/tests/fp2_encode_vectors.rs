//! Differential test of the ported `fp2_encode` against the committed
//! C-derived vectors. Record shape: (a_re, a_im) -> dst (64 bytes).
//!
//! `fp2_encode(dst, a)` writes `fp_encode(a.re)` into the low 32
//! bytes and `fp_encode(a.im)` into the high 32 bytes. Each per-half
//! encode runs `redc + modshr` (see `fp_encode`'s pin), so
//! non-canonical inputs collapse before serialisation.

use sqisign_gf::{fp2_encode, Fp, Fp2, FP2_ENCODED_BYTES, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp2_encode.json"
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
fn fp2_encode_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp2_encode");
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
        let expected_bytes = decode("dst", &v.outputs["dst"]).expect("dst hex");
        assert_eq!(
            expected_bytes.len(),
            FP2_ENCODED_BYTES,
            "dst is exactly {} bytes (vector {})",
            FP2_ENCODED_BYTES,
            v.id
        );
        let mut expected = [0u8; FP2_ENCODED_BYTES];
        expected.copy_from_slice(&expected_bytes);

        let a = Fp2 { re: a_re, im: a_im };
        let mut dst = [0u8; FP2_ENCODED_BYTES];
        fp2_encode(&mut dst, &a);
        assert_eq!(
            dst, expected,
            "vector {} diverged from the C reference at the fp2_encode boundary",
            v.id
        );
    }
}
