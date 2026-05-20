//! Differential test of `ibz_cmp` against committed C-derived vectors.
use sqisign_quaternion::{ibz_cmp, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_cmp.json"
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
fn ibz_cmp_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_cmp");
    assert!(f.vectors.len() >= 1000);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let b = read_ibz("b", &v.inputs["b"]);
        let exp = read_i8("r", &v.outputs["r"]);
        let got = ibz_cmp(&a, &b);
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
