//! Differential test of `ibz_cornacchia_prime`.
//!
//! Solving `x^2 + n y^2 == p` has a positive-integer solution (x, y) only
//! when the standard Cornacchia conditions hold. The C reference and the
//! Rust port both return 1/0 to signal solvability; on success the
//! representative (x, y) is unique up to signs and may flip between
//! implementations of the Euclidean step. We assert the success flag
//! exactly and verify `x^2 + n*y^2 == p` on the success path rather than
//! pinning the specific (x, y) the reference produced.
use sqisign_quaternion::{ibz_add, ibz_cornacchia_prime, ibz_mul, Ibz};
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/quaternion/ibz_cornacchia_prime.json"
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
fn ibz_cornacchia_prime_matches_reference_vectors() {
    let f = load(VECTORS).unwrap();
    assert_eq!(f.boundary, "sqisign_quaternion::ibz_cornacchia_prime");
    for v in &f.vectors {
        let n = read_ibz("n", &v.inputs["n"]);
        let p = read_ibz("p", &v.inputs["p"]);
        let exp_ok = read_i8("ok", &v.outputs["ok"]);
        let mut x = Ibz::zero();
        let mut y = Ibz::zero();
        let ok = ibz_cornacchia_prime(&mut x, &mut y, &n, &p) as i8;
        assert_eq!(
            ok, exp_ok,
            "vector {}: solvability flag mismatch (port={ok}, ref={exp_ok})",
            v.id
        );
        if exp_ok == 1 {
            let mut xx = Ibz::zero();
            ibz_mul(&mut xx, &x, &x);
            let mut yy = Ibz::zero();
            ibz_mul(&mut yy, &y, &y);
            let mut nyy = Ibz::zero();
            ibz_mul(&mut nyy, &n, &yy);
            let mut sum = Ibz::zero();
            ibz_add(&mut sum, &xx, &nyy);
            assert_eq!(sum.0, p.0, "vector {}: x^2 + n*y^2 != p", v.id);
        }
    }
}
