//! Differential test of the ported `mp_sub` against the committed
//! C-derived vectors. Same little-endian limb layout as `mp_add`; the
//! battery includes `a < b` cases that pin the wrapping (borrow-discard)
//! behaviour.

use sqisign_mp::mp_sub;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/mp/mp_sub.json");

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

fn le_u32(b: &[u8]) -> usize {
    assert_eq!(b.len(), 4, "nwords must be a u32");
    u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize
}

#[test]
fn mp_sub_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_sub");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let b = limbs("b", &decode("b", &v.inputs["b"]).expect("b hex"));
        let nwords = le_u32(&decode("nwords", &v.inputs["nwords"]).expect("nwords hex"));
        let expected = limbs("c", &decode("c", &v.outputs["c"]).expect("c hex"));

        assert_eq!(a.len(), nwords, "vector {}: a vs nwords", v.id);

        let mut c = vec![0u64; nwords];
        mp_sub(&mut c, &a, &b);
        assert_eq!(
            c, expected,
            "vector {} diverged from the C reference (nwords={})",
            v.id, nwords
        );
    }
}
