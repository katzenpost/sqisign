//! Differential test of the ported `fp_decode_reduce` against the
//! committed C-derived vectors. The differential boundary is the
//! variable-length input byte slice plus the explicitly recorded
//! 8-byte little-endian length, mapping to the raw internal `fp_t`
//! (five little-endian 8-byte limbs) the reference writes.
//!
//! `fp_decode_reduce` is the arbitrary-length reducer the reference
//! defines as a two-phase fold (see
//! `the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:752..791`): the
//! trailing partial block (`len % 32`) is decoded via `fp_decode` after
//! zero-padding, then each preceding 32-byte block is partial-reduced
//! via the `5 * 2^248 == 1 mod p` identity, re-encoded, decoded via
//! `fp_decode` again, and added to `d` after `d` has been multiplied by
//! `R2 == 2^256 mod p`.
//!
//! The differential battery sweeps representative lengths
//! (`0, 1, 16, 31, 32, 33, 63, 64, 100, 200`) crossed with twelve byte
//! patterns plus a 1000-record random sweep over lengths in `0..=512`.
//! The test counts how many records hit each length to ensure both the
//! empty-input early-return and the multi-block descending fold are
//! actually exercised, not just one of the two.

use std::collections::BTreeMap;

use sqisign_gf::{fp_decode_reduce, Fp, NWORDS_FIELD};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/gf/fp_decode_reduce.json"
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
fn fp_decode_reduce_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_gf::fp_decode_reduce");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    // Length-class histogram: groups records by the path inside
    // fp_decode_reduce that the recorded `len` exercises.
    let mut empty_count = 0usize;
    let mut partial_only_count = 0usize; // 0 < len < 32
    let mut exact_block_count = 0usize; // len % 32 == 0 && len > 0
    let mut multi_mixed_count = 0usize; // len >= 32 && len % 32 != 0

    for v in &file.vectors {
        let src = decode("src", &v.inputs["src"]).expect("src hex");
        let len_bytes = decode("len", &v.inputs["len"]).expect("len hex");
        assert_eq!(
            len_bytes.len(),
            8,
            "len is an 8-byte little-endian u64 (vector {})",
            v.id
        );
        let recorded_len = u64::from_le_bytes(len_bytes.try_into().unwrap()) as usize;
        assert_eq!(
            recorded_len,
            src.len(),
            "vector {}: recorded len {recorded_len} disagrees with src length {}",
            v.id,
            src.len()
        );

        let expected = fp_from("d", &decode("d", &v.outputs["d"]).expect("d hex"));

        let mut d = [0u64; NWORDS_FIELD];
        fp_decode_reduce(&mut d, &src);
        assert_eq!(
            d, expected,
            "vector {} (len {}) diverged from the C reference at the fp_decode_reduce boundary",
            v.id, recorded_len
        );

        match (recorded_len, recorded_len % 32) {
            (0, _) => empty_count += 1,
            (l, _) if l < 32 => partial_only_count += 1,
            (_, 0) => exact_block_count += 1,
            _ => multi_mixed_count += 1,
        }
    }

    // The edge battery has 12 records at len=0; the random sweep hits
    // len=0 with low probability. Pinning a non-zero count for every
    // path the reducer's structure distinguishes guards against the
    // edge patterns silently dropping out of the emitter, and ensures
    // both the partial-only and multi-block branches see real
    // coverage.
    assert!(empty_count >= 12, "expected the 12 edge records at len=0");
    assert!(
        partial_only_count > 0,
        "expected at least one partial-only record (len in 1..32)"
    );
    assert!(
        exact_block_count > 0,
        "expected at least one exact-block record (len > 0, len % 32 == 0)"
    );
    assert!(
        multi_mixed_count > 0,
        "expected at least one multi-block-with-remainder record (len > 32, len % 32 != 0)"
    );
}

#[test]
fn fp_decode_reduce_deterministic() {
    // Pins the determinism the rest of the suite relies on: a second
    // call on identical input yields identical output. partial_reduce
    // + fp_decode + fp_mul + fp_add is a fixed chain with no randomness,
    // no global state, and no hidden inputs.
    let file = load(VECTORS).expect("vector file must load");
    for v in &file.vectors {
        let src = decode("src", &v.inputs["src"]).expect("src hex");
        let mut d1 = [0u64; NWORDS_FIELD];
        let mut d2 = [0u64; NWORDS_FIELD];
        fp_decode_reduce(&mut d1, &src);
        fp_decode_reduce(&mut d2, &src);
        assert_eq!(
            d1, d2,
            "fp_decode_reduce not deterministic on input from vector {}",
            v.id
        );
    }
}

#[test]
fn fp_decode_reduce_length_distribution() {
    // Pin the length distribution of the edge battery (10 lengths from
    // {0, 1, 16, 31, 32, 33, 63, 64, 100, 200} x 12 patterns = 120 edge
    // records, each length appearing exactly 12 times). The sweep adds
    // 1000 random-length records over 0..=512, so each edge length
    // shows up at least 12 times. Pinning the minimum guards against
    // the edge battery silently dropping a length.
    let file = load(VECTORS).expect("vector file must load");
    let mut histogram: BTreeMap<usize, usize> = BTreeMap::new();
    for v in &file.vectors {
        let len_bytes = decode("len", &v.inputs["len"]).expect("len hex");
        let recorded_len = u64::from_le_bytes(len_bytes.try_into().unwrap()) as usize;
        *histogram.entry(recorded_len).or_insert(0) += 1;
    }
    for edge_len in [0usize, 1, 16, 31, 32, 33, 63, 64, 100, 200] {
        let count = histogram.get(&edge_len).copied().unwrap_or(0);
        assert!(
            count >= 12,
            "edge length {edge_len} should appear at least 12 times (the 12-pattern cross-product); got {count}"
        );
    }
}
