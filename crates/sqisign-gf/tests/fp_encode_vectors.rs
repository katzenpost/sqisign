//! Differential test of the ported `fp_encode` against the committed
//! C-derived vectors. The differential boundary is the raw internal
//! `fp_t` input (five little-endian 8-byte limbs) plus the 32-byte
//! little-endian output buffer the reference writes to `void *dst`.
//!
//! `fp_encode` is the canonical serialization the reference defines as
//! `redc` followed by a 32-iteration `c[0] & 0xff; modshr(8, c)` loop
//! (see `the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:664..675`).
//! Two non-canonical encodings of the same field element yield the same
//! 32 bytes; the edge battery in
//! `tools/cdump/src/dump_main.c::emit_fp_encode_edges` includes the
//! radix-2^51 encoding of `p` (which `redc`s to zero) alongside the
//! canonical zero to exercise the redc-first structure.

use sqisign_gf::{fp_encode, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_encode.json"
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
fn fp_encode_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_encode");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected_bytes = decode("dst", &v.outputs["dst"]).expect("dst hex");
        assert_eq!(
            expected_bytes.len(),
            32,
            "dst is exactly 32 bytes (vector {})",
            v.id
        );
        let mut expected = [0u8; 32];
        expected.copy_from_slice(&expected_bytes);

        let mut dst = [0u8; 32];
        fp_encode(&mut dst, &a);
        assert_eq!(
            dst, expected,
            "vector {} diverged from the C reference at the fp_encode boundary",
            v.id
        );
    }
}

#[test]
fn fp_encode_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical bytes. redc + modshr is
    // a fixed chain with no randomness, no global state, and no hidden
    // inputs.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let a = fp_from("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let mut d1 = [0u8; 32];
        let mut d2 = [0u8; 32];
        fp_encode(&mut d1, &a);
        fp_encode(&mut d2, &a);
        assert_eq!(
            d1, d2,
            "fp_encode not deterministic on input from vector {}",
            v.id
        );
    }
}
