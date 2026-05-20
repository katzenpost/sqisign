//! Differential test of `ibz_is_zero` against committed C-derived vectors.
use sqisign_quaternion::{ibz_is_zero, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_is_zero.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_i8(l: &str, h: &str) -> i8 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 1);
    b[0] as i8
}

#[test]
fn ibz_is_zero_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_is_zero");
    assert!(f.vectors.len() >= 200);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let exp = read_i8("r", &v.outputs["r"]);
        let got = ibz_is_zero(&a);
        let got_norm: i8 = if got > 0 {
            1
        } else if got < 0 {
            -1
        } else {
            0
        };
        assert_eq!(got_norm, exp, "vector {}", v.id);
    }
}
