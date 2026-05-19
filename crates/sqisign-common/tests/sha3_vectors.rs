//! Differential test of the ported fixed-output SHA3 digests against the
//! committed C-derived vectors. Same harness as the SHAKE boundaries; the
//! fixed-output records carry no `outlen` field because the digest size is
//! intrinsic to each function.

use sqisign_common::{sha3_256, sha3_384, sha3_512};
use sqisign_vectors::{decode, load};

fn check(path: &str, name: &str, digest: impl Fn(&[u8]) -> Vec<u8>) {
    let file = load(path).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, name);
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );
    for v in &file.vectors {
        let input = decode("input", &v.inputs["input"]).expect("input hex");
        let expected = decode("output", &v.outputs["output"]).expect("output hex");
        assert_eq!(
            digest(&input),
            expected,
            "vector {} diverged from the C reference (inlen={})",
            v.id,
            input.len()
        );
    }
}

const DIR: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/common/");

#[test]
fn sha3_256_matches_reference_vectors() {
    check(
        &format!("{DIR}sha3_256.json"),
        "sqisign_common::sha3_256",
        |i| sha3_256(i).to_vec(),
    );
}

#[test]
fn sha3_384_matches_reference_vectors() {
    check(
        &format!("{DIR}sha3_384.json"),
        "sqisign_common::sha3_384",
        |i| sha3_384(i).to_vec(),
    );
}

#[test]
fn sha3_512_matches_reference_vectors() {
    check(
        &format!("{DIR}sha3_512.json"),
        "sqisign_common::sha3_512",
        |i| sha3_512(i).to_vec(),
    );
}
