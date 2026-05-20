//! Differential test of `ibz_mul` against committed C-derived vectors.
use sqisign_quaternion::{ibz_mul, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mul.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}

#[test]
fn ibz_mul_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mul");
    assert!(f.vectors.len() >= 1000);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let b = read_ibz("b", &v.inputs["b"]);
        let expected = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        ibz_mul(&mut r, &a, &b);
        assert_eq!(r.0, expected.0, "vector {}", v.id);
    }
}
