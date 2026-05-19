//! Differential test of the ported `mp_is_zero` against the committed
//! C-derived vectors. Result is a single byte 0/1; also independently
//! checked to equal "all limbs zero".

use sqisign_mp::mp_is_zero;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_is_zero.json"
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

#[test]
fn mp_is_zero_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_is_zero");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let a = limbs("a", &decode("a", &v.inputs["a"]).expect("a hex"));
        let expected = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(expected.len(), 1, "result is a single byte");
        let want = expected[0] != 0;

        let got = mp_is_zero(&a);
        assert_eq!(got, want, "vector {} diverged from the reference", v.id);
        assert_eq!(
            got,
            a.iter().all(|&x| x == 0),
            "vector {}: not the all-zero predicate (upstream changed?)",
            v.id
        );
    }
}
