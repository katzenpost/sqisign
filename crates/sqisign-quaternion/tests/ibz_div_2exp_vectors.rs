//! Differential test of `ibz_div_2exp` against committed C-derived vectors.
use sqisign_quaternion::{ibz_div_2exp, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_div_2exp.json"
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
fn ibz_div_2exp_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_div_2exp");
    assert!(f.vectors.len() >= 1000);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let exp = read_le_u32("exp", &v.inputs["exp"]);
        let expected = read_ibz("q", &v.outputs["q"]);
        let mut q = Ibz::zero();
        ibz_div_2exp(&mut q, &a, exp);
        assert_eq!(q.0, expected.0, "vector {}", v.id);
    }
}
