//! Differential test of the ported `mp_is_one` against the committed
//! C-derived vectors. Result is a single byte 0/1; also independently
//! checked to equal `x[0] == 1 && x[1..] all zero`.

use sqisign_mp::mp_is_one;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/mp/mp_is_one.json"
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
fn mp_is_one_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_mp::mp_is_one");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let x = limbs("x", &decode("x", &v.inputs["x"]).expect("x hex"));
        let expected = decode("result", &v.outputs["result"]).expect("result hex");
        assert_eq!(expected.len(), 1, "result is a single byte");
        let want = expected[0] != 0;

        let got = mp_is_one(&x);
        assert_eq!(got, want, "vector {} diverged from the reference", v.id);
        assert_eq!(
            got,
            x[0] == 1 && x[1..].iter().all(|&t| t == 0),
            "vector {}: not the canonical-one predicate (upstream changed?)",
            v.id
        );
    }
}
