//! Differential test of `ibz_two_adic`.
use sqisign_quaternion::{ibz_two_adic, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_two_adic.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_le_i32(l: &str, h: &str) -> i32 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 4);
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

#[test]
fn ibz_two_adic_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_two_adic");
    assert!(f.vectors.len() >= 200);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let exp = read_le_i32("r", &v.outputs["r"]);
        let got = ibz_two_adic(&a);
        assert_eq!(got, exp, "vector {}", v.id);
    }
}
