//! Differential test of `ibz_div` against committed C-derived vectors.
use sqisign_quaternion::{ibz_div, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_div.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}

#[test]
fn ibz_div_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_div");
    assert!(f.vectors.len() >= 1000);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let b = read_ibz("b", &v.inputs["b"]);
        let exp_q = read_ibz("q", &v.outputs["q"]);
        let exp_r = read_ibz("r", &v.outputs["r"]);
        let mut q = Ibz::zero();
        let mut r = Ibz::zero();
        ibz_div(&mut q, &mut r, &a, &b);
        assert_eq!(q.0, exp_q.0, "vector {} q", v.id);
        assert_eq!(r.0, exp_r.0, "vector {} r", v.id);
    }
}
