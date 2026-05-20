//! Differential test of the ported `ibz_add` against the committed
//! C-derived vectors. The boundary is value-level: the C side serialises
//! each `ibz_t` in canonical signed-bytes form
//! (1-byte sign, u32 LE length, big-endian magnitude). The port reads the
//! inputs, computes via `num-bigint`, and compares the canonical
//! re-serialisation of its result to the recorded output bytes.

use sqisign_quaternion::{ibz_add, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_add.json"
);

fn read_ibz(label: &str, hex: &str) -> Ibz {
    let bytes = decode(label, hex).expect("hex decode");
    Ibz::from_canonical_bytes(&bytes).expect("canonical bytes parse")
}

#[test]
fn ibz_add_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load");
    assert_eq!(file.boundary, "sqisign_quaternion::ibz_add");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );
    for v in &file.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let b = read_ibz("b", &v.inputs["b"]);
        let expected = read_ibz("r", &v.outputs["r"]);
        let mut got = Ibz::zero();
        ibz_add(&mut got, &a, &b);
        assert_eq!(got.0, expected.0, "vector {} diverged", v.id);
    }
}
