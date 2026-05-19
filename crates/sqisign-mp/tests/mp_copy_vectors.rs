//! Differential test of the ported `mp_copy` against the committed
//! C-derived vectors. A plain copy; every vector has `b == a`. Thin, but
//! it keeps the boundary under the same oracle discipline as the rest.

use sqisign_mp::mp_copy;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../vectors/mp/mp_copy.json");

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

#[test]
fn mp_copy_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_copy");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = limbs("b", &decode("b", &v.outputs["b"]).expect("b hex"));

        let mut b = vec![0xdead_u64; a.len()];
        mp_copy(&mut b, &a);
        assert_eq!(b, expected, "vector {} diverged from the reference", v.id);
        assert_eq!(b, a, "vector {}: mp_copy must be the identity", v.id);
    }
}
