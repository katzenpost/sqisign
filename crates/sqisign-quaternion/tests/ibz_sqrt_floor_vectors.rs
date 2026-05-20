//! Differential test of `ibz_sqrt_floor`.
use sqisign_quaternion::{ibz_sqrt_floor, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_sqrt_floor.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}

#[test]
fn ibz_sqrt_floor_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_sqrt_floor");
    assert!(f.vectors.len() >= 50);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let exp = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        ibz_sqrt_floor(&mut r, &a);
        assert_eq!(r.0, exp.0, "vector {}", v.id);
    }
}
