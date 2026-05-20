//! Differential test of `ibz_cmp_int32`.
use sqisign_quaternion::{ibz_cmp_int32, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_cmp_int32.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_le_i32(l: &str, h: &str) -> i32 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 4);
    i32::from_le_bytes([b[0], b[1], b[2], b[3]])
}
fn read_i8(l: &str, h: &str) -> i8 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 1);
    b[0] as i8
}

#[test]
fn ibz_cmp_int32_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_cmp_int32");
    assert!(f.vectors.len() >= 500);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let y = read_le_i32("y", &v.inputs["y"]);
        let exp = read_i8("r", &v.outputs["r"]);
        let got = ibz_cmp_int32(&a, y);
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
