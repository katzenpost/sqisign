//! Property tests for `mp_is_zero`: the all-zero predicate.

use proptest::prelude::*;
use sqisign_mp::mp_is_zero;

proptest! {
    // An all-zero buffer of any length is zero.
    #[test]
    fn zeros_are_zero(n in 0usize..256) {
        prop_assert!(mp_is_zero(&vec![0u64; n]));
    }

    // Any buffer with a nonzero limb is not zero.
    #[test]
    fn nonzero_limb_is_not_zero(mut a in proptest::collection::vec(any::<u64>(), 1..128),
                                idx in 0usize..128, val in 1u64..) {
        let i = idx % a.len();
        a[i] = val;
        prop_assert!(!mp_is_zero(&a));
    }

    // Matches the direct predicate over arbitrary input.
    #[test]
    fn matches_all_zero(a in proptest::collection::vec(any::<u64>(), 0..256)) {
        prop_assert_eq!(mp_is_zero(&a), a.iter().all(|&x| x == 0));
    }
}
