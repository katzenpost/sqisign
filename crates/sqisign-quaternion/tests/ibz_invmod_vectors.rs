//! Differential test of `ibz_invmod`.
use sqisign_quaternion::{ibz_invmod, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_invmod.json"
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
fn ibz_invmod_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_invmod");
    assert!(f.vectors.len() >= 100);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let m = read_ibz("m", &v.inputs["m"]);
        let exp_ok = read_i8("ok", &v.outputs["ok"]);
        let exp_inv = read_ibz("inv", &v.outputs["inv"]);
        let mut inv = Ibz::zero();
        let ok = ibz_invmod(&mut inv, &a, &m) as i8;
        assert_eq!(ok, exp_ok, "vector {}: ok flag", v.id);
        if exp_ok == 1 {
            assert_eq!(inv.0, exp_inv.0, "vector {}: inv value", v.id);
        }
        // When ok == 0, the C reference leaves `inv` in scratch state.
        // We do not pin its value.
    }
}
