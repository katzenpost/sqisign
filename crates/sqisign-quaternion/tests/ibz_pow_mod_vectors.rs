//! Differential test of `ibz_pow_mod`.
use sqisign_quaternion::{ibz_pow_mod, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_pow_mod.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}

#[test]
fn ibz_pow_mod_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_pow_mod");
    assert!(f.vectors.len() >= 500);
    for v in &f.vectors {
        let x = read_ibz("x", &v.inputs["x"]);
        let e = read_ibz("e", &v.inputs["e"]);
        let m = read_ibz("m", &v.inputs["m"]);
        let expected = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        ibz_pow_mod(&mut r, &x, &e, &m);
        assert_eq!(r.0, expected.0, "vector {}", v.id);
    }
}
