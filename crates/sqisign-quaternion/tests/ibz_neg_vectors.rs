//! Differential test of `ibz_neg` against committed C-derived vectors.
use sqisign_quaternion::{ibz_neg, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_neg.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}

#[test]
fn ibz_neg_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_neg");
    assert!(f.vectors.len() >= 200);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let expected = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        ibz_neg(&mut r, &a);
        assert_eq!(r.0, expected.0, "vector {}", v.id);
    }
}
