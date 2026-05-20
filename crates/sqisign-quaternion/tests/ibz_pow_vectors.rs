//! Differential test of `ibz_pow`.
use sqisign_quaternion::{ibz_pow, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_pow.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_le_u32(l: &str, h: &str) -> u32 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 4);
    u32::from_le_bytes([b[0], b[1], b[2], b[3]])
}

#[test]
fn ibz_pow_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_pow");
    assert!(f.vectors.len() >= 200);
    for v in &f.vectors {
        let x = read_ibz("x", &v.inputs["x"]);
        let e = read_le_u32("e", &v.inputs["e"]);
        let expected = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        ibz_pow(&mut r, &x, e);
        assert_eq!(r.0, expected.0, "vector {}", v.id);
    }
}
