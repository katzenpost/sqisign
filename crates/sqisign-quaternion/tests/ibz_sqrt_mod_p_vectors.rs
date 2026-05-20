//! Differential test of `ibz_sqrt_mod_p`.
//!
//! The reference dispatches between three closed forms (`p == 3 mod 4`,
//! `p == 5 mod 8`) and Tonelli-Shanks (`p == 1 mod 8`); a square root mod
//! `p` is not unique (s and p-s are both valid), so for the success path
//! we verify `r^2 == a (mod p)` rather than pinning the specific
//! representative the reference happens to produce. The Rust port returns
//! the same representative under the same dispatch, but we keep the test
//! tolerant so an upstream switch of branch (e.g. preferring s vs p-s)
//! does not break this gate.
use sqisign_quaternion::{ibz_legendre, ibz_mod, ibz_mul, ibz_sqrt_mod_p, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_sqrt_mod_p.json"
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
fn ibz_sqrt_mod_p_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_sqrt_mod_p");
    assert!(f.vectors.len() >= 200);
    for v in &f.vectors {
        let a = read_ibz("a", &v.inputs["a"]);
        let p = read_ibz("p", &v.inputs["p"]);
        let exp_ok = read_i8("ok", &v.outputs["ok"]);
        let mut r = Ibz::zero();
        let ok = ibz_sqrt_mod_p(&mut r, &a, &p) as i8;
        assert_eq!(
            ok,
            exp_ok,
            "vector {}: ok mismatch (Legendre(a, p) = {})",
            v.id,
            ibz_legendre(&a, &p)
        );
        if exp_ok == 1 {
            // r^2 == a (mod p)
            let mut r2 = Ibz::zero();
            ibz_mul(&mut r2, &r, &r);
            let mut r2_mod = Ibz::zero();
            ibz_mod(&mut r2_mod, &r2, &p);
            let mut a_mod = Ibz::zero();
            ibz_mod(&mut a_mod, &a, &p);
            assert_eq!(r2_mod.0, a_mod.0, "vector {}: r^2 != a mod p", v.id);
        }
    }
}
