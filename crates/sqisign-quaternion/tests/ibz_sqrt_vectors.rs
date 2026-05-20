//! Differential test of `ibz_sqrt`.
use sqisign_quaternion::{ibz_sqrt, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_sqrt.json"
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
fn ibz_sqrt_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_sqrt");
    assert!(f.vectors.len() >= 50);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let exp_ok = read_i8("ok", &v.outputs["ok"]);
        let exp_r = read_ibz("r", &v.outputs["r"]);
        let mut r = Ibz::zero();
        let ok = ibz_sqrt(&mut r, &a) as i8;
        assert_eq!(ok, exp_ok, "vector {}: ok", v.id);
        if exp_ok == 1 {
            assert_eq!(r.0, exp_r.0, "vector {}: r value", v.id);
        }
    }
}
