//! Property tests for `mp_mul`.
//!
//! `mp_mul` is the low-half product. For `nwords >= 2` it is exactly
//! `(a*b) mod 2^(64n)` and obeys the ring laws. For `nwords == 1` it
//! faithfully reproduces the reference's column-0 double-count, so the
//! single-limb invariant is the *defect's* law, `2*(a*b) mod 2^64`, not
//! ordinary multiplication. Both are pinned here so a regression in
//! either regime is caught.

use proptest::prelude::*;
use sqisign_mp::{mp_add, mp_mul};

fn mul(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; a.len()];
    mp_mul(&mut c, a, b);
    c
}

proptest! {
    // Single limb: the reproduced upstream defect is exactly 2*(a*b).
    #[test]
    fn single_limb_is_doubled_product(x in any::<u64>(), y in any::<u64>()) {
        let low = (x as u128 * y as u128) as u64;
        prop_assert_eq!(mul(&[x], &[y]), vec![low.wrapping_mul(2)]);
    }

    // nwords >= 2: commutative.
    #[test]
    fn multilimb_commutes(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 2..24)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        prop_assert_eq!(mul(&a, &b), mul(&b, &a));
    }

    // nwords >= 2: multiplying by one is the identity (low half).
    #[test]
    fn multilimb_times_one(a in proptest::collection::vec(any::<u64>(), 2..24)) {
        let mut one = vec![0u64; a.len()];
        one[0] = 1;
        prop_assert_eq!(mul(&a, &one), a.clone());
    }

    // nwords >= 2: left-distributes over addition in Z/2^(64n).
    #[test]
    fn multilimb_distributes(v in proptest::collection::vec(
            (any::<u64>(), any::<u64>(), any::<u64>()), 2..16)) {
        let a: Vec<u64> = v.iter().map(|&(x, _, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y, _)| y).collect();
        let c: Vec<u64> = v.iter().map(|&(_, _, z)| z).collect();
        let mut bc = vec![0u64; a.len()];
        mp_add(&mut bc, &b, &c);
        let lhs = mul(&a, &bc);
        let mut rhs = vec![0u64; a.len()];
        mp_add(&mut rhs, &mul(&a, &b), &mul(&a, &c));
        prop_assert_eq!(lhs, rhs);
    }

    // nwords >= 2: agrees with a u128 reference for two limbs.
    #[test]
    fn two_limb_is_u128_low_half(a0 in any::<u64>(), a1 in any::<u64>(),
                                 b0 in any::<u64>(), b1 in any::<u64>()) {
        let a = (a0 as u128) | ((a1 as u128) << 64);
        let b = (b0 as u128) | ((b1 as u128) << 64);
        let prod = a.wrapping_mul(b); // low 128 bits
        let c = mul(&[a0, a1], &[b0, b1]);
        prop_assert_eq!(c[0], prod as u64);
        prop_assert_eq!(c[1], (prod >> 64) as u64);
    }
}
