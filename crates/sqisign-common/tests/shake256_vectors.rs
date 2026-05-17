//! Phase 0 end-to-end proof of the differential harness.
//!
//! This test closes the loop the plan describes:
//!
//! ```text
//! instrumented C reference -> raw dump -> vector-gen -> canonical JSON
//!   -> sqisign-vectors::load -> bit-compare here
//! ```
//!
//! There is deliberately no `sqisign-common` code under test yet. The oracle
//! is the known-good `sha3` crate. When `common`'s own SHAKE lands in
//! Phase 1, only the `oracle_shake256` body changes; the harness, the
//! vectors, and this assertion stay exactly as they are. That is the point:
//! the scaffolding is proven correct before any primitive depends on it.

use sha3::digest::{ExtendableOutput, Update, XofReader};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/common/shake256.json"
);

fn oracle_shake256(input: &[u8], out_len: usize) -> Vec<u8> {
    let mut hasher = sha3::Shake256::default();
    hasher.update(input);
    let mut out = vec![0u8; out_len];
    hasher.finalize_xof().read(&mut out);
    out
}

fn le_usize(bytes: &[u8]) -> usize {
    let mut acc = 0u64;
    for (i, b) in bytes.iter().enumerate() {
        acc |= (*b as u64) << (8 * i);
    }
    acc as usize
}

#[test]
fn shake256_matches_reference_vectors() {
    // `load` enforces canonical form and the UPSTREAM.md pin before we see
    // any data; a stale or hand-edited vector file fails here, loudly.
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");

    assert_eq!(file.boundary, "sqisign_common::shake256");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    let mut checked = 0usize;
    for v in &file.vectors {
        let input = decode("input", &v.inputs["input"]).expect("input hex");
        let outlen = le_usize(&decode("outlen", &v.inputs["outlen"]).expect("outlen hex"));
        let expected = decode("output", &v.outputs["output"]).expect("output hex");

        assert_eq!(
            expected.len(),
            outlen,
            "vector {}: recorded output length disagrees with recorded outlen",
            v.id
        );

        let got = oracle_shake256(&input, outlen);
        assert_eq!(
            got,
            expected,
            "vector {} diverged from the C reference (inlen={}, outlen={})",
            v.id,
            input.len(),
            outlen
        );
        checked += 1;
    }
    assert_eq!(checked, file.vectors.len());
}
