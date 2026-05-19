//! Differential test of the ported `select_ct` against the committed
//! C-derived vectors. The battery includes arbitrary masks (not just
//! 0 / all-ones), pinning the exact bitwise blend; independently checked
//! to equal `((a^b)&mask)^a`.

use sqisign_mp::select_ct;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/select_ct.json"
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

fn le_u64(b: &[u8]) -> u64 {
    assert_eq!(b.len(), 8, "mask is a u64");
    u64::from_le_bytes(b.try_into().unwrap())
}

#[test]
fn select_ct_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::select_ct");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = limbs("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let mask = le_u64(&decode("mask", &v.inputs["mask"]).expect("mask hex"));
        let expected = limbs("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        let mut c = vec![0xdead_u64; a.len()];
        select_ct(&mut c, &a, &b, mask);
        assert_eq!(c, expected, "vector {} diverged from the reference", v.id);

        let blend: Vec<u64> = (0..a.len())
            .map(|i| ((a[i] ^ b[i]) & mask) ^ a[i])
            .collect();
        assert_eq!(
            c, blend,
            "vector {}: not ((a^b)&mask)^a (upstream changed?)",
            v.id
        );
    }
}
