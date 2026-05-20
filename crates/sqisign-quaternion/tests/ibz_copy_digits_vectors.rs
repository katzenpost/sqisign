//! Differential test of `ibz_copy_digits`.
use sqisign_quaternion::{ibz_copy_digits, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_copy_digits.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_le_u32(l: &str, h: &str) -> u32 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 4);
    u32::from_le_bytes([b[0], b[1], b[2], b[3]])
}
fn read_words(l: &str, h: &str) -> Vec<u64> {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len() % 8, 0);
    b.chunks_exact(8)
        .map(|c| u64::from_le_bytes(c.try_into().unwrap()))
        .collect()
}

#[test]
fn ibz_copy_digits_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_copy_digits");
    assert!(f.vectors.len() >= 100);
    for v in &f.vectors {
        let dig = read_words("dig", &v.inputs["dig"]);
        let nwords = read_le_u32("nwords", &v.inputs["nwords"]) as usize;
        let expected = read_ibz("r", &v.outputs["r"]);
        assert_eq!(dig.len(), nwords, "vector {}: dig len != nwords", v.id);
        let mut r = Ibz::zero();
        ibz_copy_digits(&mut r, &dig);
        assert_eq!(r.0, expected.0, "vector {}", v.id);
    }
}
