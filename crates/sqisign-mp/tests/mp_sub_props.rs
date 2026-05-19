//! Property tests for `mp_sub`, the inverse of `mp_add` in Z/2^(64*n).

use proptest::prelude::*;
use sqisign_mp::{mp_add, mp_sub};

fn add(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; a.len()];
    mp_add(&mut c, a, b);
    c
}
fn sub(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; a.len()];
    mp_sub(&mut c, a, b);
    c
}

proptest! {
    // (a + b) - b == a in the ring (add and sub are mutual inverses).
    #[test]
    fn add_sub_roundtrip(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        prop_assert_eq!(sub(&add(&a, &b), &b), a);
    }

    // a - a == 0.
    #[test]
    fn self_difference_is_zero(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        prop_assert_eq!(sub(&a, &a), vec![0u64; a.len()]);
    }

    // a - 0 == a.
    #[test]
    fn minus_zero_is_identity(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        let zero = vec![0u64; a.len()];
        prop_assert_eq!(sub(&a, &zero), a.clone());
    }

    // Single limb agrees with native wrapping subtraction.
    #[test]
    fn single_limb_is_wrapping_sub(x in any::<u64>(), y in any::<u64>()) {
        prop_assert_eq!(sub(&[x], &[y]), vec![x.wrapping_sub(y)]);
    }

    // 0 - a == two's-complement negation: (0 - a) + a == 0.
    #[test]
    fn negation_closes(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        let zero = vec![0u64; a.len()];
        let neg = sub(&zero, &a);
        prop_assert_eq!(add(&neg, &a), zero);
    }
}
