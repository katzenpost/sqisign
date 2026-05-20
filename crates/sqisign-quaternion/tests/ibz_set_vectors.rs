//! Differential test of `ibz_set`.
use sqisign_quaternion::{ibz_set, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_set.json"
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
fn ibz_set_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_set");
    assert!(f.vectors.len() >= 50);
    for v in &f.vectors {
        let x = read_le_i32("x", &v.inputs["x"]);
        let expected = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        ibz_set(&mut r, x);
        assert_eq!(r.0, expected.0, "vector {}", v.id);
    }
}
