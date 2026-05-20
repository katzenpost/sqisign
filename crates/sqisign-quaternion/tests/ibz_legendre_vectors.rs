//! Differential test of `ibz_legendre`.
use sqisign_quaternion::{ibz_legendre, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_legendre.json"
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
fn ibz_legendre_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_legendre");
    assert!(f.vectors.len() >= 100);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let p = read_ibz("p", &v.inputs["p"]);
        let exp = read_i8("r", &v.outputs["r"]);
        let got = ibz_legendre(&a, &p) as i8;
        assert_eq!(got, exp, "vector {}", v.id);
    }
}
