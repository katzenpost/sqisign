//! Differential test of the ported NIST CTR-DRBG against the committed
//! C-derived vectors.
//!
//! Each vector records the 48-byte entropy, an optional personalization
//! string (a zero-length `pers` blob means none), and the exact sequence of
//! request lengths the reference was driven with (`req_splits`, packed
//! little-endian u32). The replay seeds an identical instance and services
//! the identical request sequence, bit-comparing the concatenated output:
//! it is the stateful DRBG path, evolving after every draw, that is tested,
//! not a single call.

use sqisign_common::CtrDrbg;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/common/ctr_drbg.json"
);

fn unpack_splits(bytes: &[u8]) -> Vec<usize> {
    assert_eq!(bytes.len() % 4, 0, "split list not a whole number of u32");
    bytes
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize)
        .collect()
}

fn fixed48(label: &str, v: &[u8]) -> [u8; 48] {
    let mut a = [0u8; 48];
    assert_eq!(v.len(), 48, "{label} must be 48 bytes, got {}", v.len());
    a.copy_from_slice(v);
    a
}

#[test]
fn ctr_drbg_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_common::ctr_drbg");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let entropy = fixed48(
            "entropy",
            &decode("entropy", &v.inputs["entropy"]).expect("entropy hex"),
        );
        let pers_bytes = decode("pers", &v.inputs["pers"]).expect("pers hex");
        let splits =
            unpack_splits(&decode("req_splits", &v.inputs["req_splits"]).expect("req_splits hex"));
        let expected = decode("output", &v.outputs["output"]).expect("output hex");

        let pers = if pers_bytes.is_empty() {
            None
        } else {
            Some(fixed48("pers", &pers_bytes))
        };
        let mut drbg = CtrDrbg::new(&entropy, pers.as_ref());

        let mut got = vec![0u8; splits.iter().sum()];
        let mut off = 0;
        for &n in &splits {
            drbg.fill(&mut got[off..off + n]);
            off += n;
        }
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference CTR-DRBG",
            v.id
        );
    }
}
