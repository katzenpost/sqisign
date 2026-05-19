//! Differential test of the ported incremental SHAKE API
//! (`Shake256Absorb`/`Shake128Absorb`) against the committed C-derived
//! vectors.
//!
//! Each vector records the *exact* absorb and squeeze fragmentation the C
//! reference's incremental API was driven with (`absorb_splits`,
//! `squeeze_splits`, packed little-endian u32). The replay here feeds the
//! identical chunking through our incremental API and bit-compares: this
//! tests the incremental code path itself, not merely the one-shot result.
//!
//! Both inc-ctx structs are layout-identical in C and the two Rust types
//! differ only in Keccak rate, so a single helper drives both boundaries.

use sqisign_common::{Shake128Absorb, Shake256Absorb};
use sqisign_vectors::{decode, load, VectorFile};

fn unpack_splits(bytes: &[u8]) -> Vec<usize> {
    assert_eq!(bytes.len() % 4, 0, "split list not a whole number of u32");
    bytes
        .chunks_exact(4)
        .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]]) as usize)
        .collect()
}

fn load_boundary(path: &str, name: &str) -> VectorFile {
    let file = load(path).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, name);
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );
    file
}

/// Replays one vector's recorded absorb/squeeze chunking through `absorb`
/// (a closure picking the right inc type) and returns the squeezed bytes.
fn replay<F>(file: &VectorFile, drive: F)
where
    F: Fn(&[u8], &[usize], &[usize]) -> Vec<u8>,
{
    for v in &file.vectors {
        let input = decode("input", &v.inputs["input"]).expect("input hex");
        let asplits = unpack_splits(
            &decode("absorb_splits", &v.inputs["absorb_splits"]).expect("absorb_splits hex"),
        );
        let ssplits = unpack_splits(
            &decode("squeeze_splits", &v.inputs["squeeze_splits"]).expect("squeeze_splits hex"),
        );
        let expected = decode("output", &v.outputs["output"]).expect("output hex");

        assert_eq!(
            asplits.iter().sum::<usize>(),
            input.len(),
            "vector {}: absorb splits do not sum to input length",
            v.id
        );
        assert_eq!(
            ssplits.iter().sum::<usize>(),
            expected.len(),
            "vector {}: squeeze splits do not sum to output length",
            v.id
        );

        let got = drive(&input, &asplits, &ssplits);
        assert_eq!(
            got, expected,
            "vector {} diverged from the C reference incremental path",
            v.id
        );
    }
}

const S256: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/common/shake256_inc.json"
);
const S128: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/common/shake128_inc.json"
);

#[test]
fn shake256_inc_matches_reference_vectors() {
    let file = load_boundary(S256, "sqisign_common::shake256_inc");
    replay(&file, |input, asplits, ssplits| {
        let mut a = Shake256Absorb::new();
        let mut off = 0;
        for &n in asplits {
            a.absorb(&input[off..off + n]);
            off += n;
        }
        let mut sq = a.finalize();
        let mut out = vec![0u8; ssplits.iter().sum()];
        let mut off = 0;
        for &n in ssplits {
            sq.squeeze(&mut out[off..off + n]);
            off += n;
        }
        out
    });
}

#[test]
fn shake128_inc_matches_reference_vectors() {
    let file = load_boundary(S128, "sqisign_common::shake128_inc");
    replay(&file, |input, asplits, ssplits| {
        let mut a = Shake128Absorb::new();
        let mut off = 0;
        for &n in asplits {
            a.absorb(&input[off..off + n]);
            off += n;
        }
        let mut sq = a.finalize();
        let mut out = vec![0u8; ssplits.iter().sum()];
        let mut off = 0;
        for &n in ssplits {
            sq.squeeze(&mut out[off..off + n]);
            off += n;
        }
        out
    });
}
