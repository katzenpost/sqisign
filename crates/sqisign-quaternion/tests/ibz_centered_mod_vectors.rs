//! Differential test of `ibz_centered_mod`.
mod common;
use common::{ibz_eq, read_ibz};
use sqisign_quaternion::{ibz_centered_mod, Ibz};
use sqisign_vectors::load;

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_centered_mod.json"
);

#[test]
fn ibz_centered_mod_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_centered_mod");
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs);
        let m = read_ibz("m", &v.inputs);
        let exp = read_ibz("r", &v.outputs);
        let mut r = Ibz::zero();
        ibz_centered_mod(&mut r, &a, &m);
        assert!(ibz_eq(&r, &exp), "vector {}", v.id);
    }
}
