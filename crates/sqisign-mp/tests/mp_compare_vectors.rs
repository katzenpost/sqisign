//! Differential test of the ported `mp_compare` against the committed
//! C-derived vectors. Result is the int -1/0/1 (4-byte little-endian
//! two's-complement). Also independently checks it equals `sign(a - b)`
//! over the little-endian values, so a future upstream re-pin changing
//! the convention would be noticed.

use sqisign_mp::mp_compare;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_compare.json"
);

fn limbs(label: &str, bytes: &[u8]) -> Vec<u64> {
    assert_eq!(
        bytes.len() % 8,
        0,
        "{label} not a whole number of u64 limbs"
    );
    bytes
        .chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

fn le_i32(b: &[u8]) -> i32 {
    assert_eq!(b.len(), 4, "result must be an i32");
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

/// sign(a - b) over equal-length little-endian values, scanning high down.
fn sign(a: &[u64], b: &[u64]) -> i32 {
    for i in (0..a.len()).rev() {
        if a[i] > b[i] {
            return 1;
        }
        if a[i] < b[i] {
            return -1;
        }
    }
    0
}

#[test]
fn mp_compare_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_compare");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = limbs("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let expected = le_i32(&decode("result", &v.outputs["result"]).expect("result hex"));

        let got = mp_compare(&a, &b);
        assert_eq!(got, expected, "vector {} diverged from the reference", v.id);
        assert_eq!(
            got,
            sign(&a, &b),
            "vector {}: not sign(a-b) (upstream convention changed?)",
            v.id
        );
    }
}
