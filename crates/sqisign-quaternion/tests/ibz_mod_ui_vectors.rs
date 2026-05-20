//! Differential test of `ibz_mod_ui`.
use sqisign_quaternion::{ibz_mod_ui, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_mod_ui.json"
);

fn read_ibz(l: &str, h: &str) -> Ibz {
    Ibz::from_canonical_bytes(&decode(l, h).unwrap()).unwrap()
}
fn read_le_u64(l: &str, h: &str) -> u64 {
    let b = decode(l, h).unwrap();
    assert_eq!(b.len(), 8);
    u64::from_le_bytes(b.try_into().unwrap())
}

#[test]
fn ibz_mod_ui_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_mod_ui");
    assert!(f.vectors.len() >= 1000);
    for v in &f.vectors {
        let n = read_ibz("n", &v.inputs["n"]);
        let d = read_le_u64("d", &v.inputs["d"]);
        let expected = read_le_u64("r", &v.outputs["r"]);
        let got = ibz_mod_ui(&n, d);
        assert_eq!(got, expected, "vector {}", v.id);
    }
}
