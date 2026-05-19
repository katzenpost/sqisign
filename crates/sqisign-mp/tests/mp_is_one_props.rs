//! Property tests for `mp_is_one`: the canonical-one predicate.

use proptest::prelude::*;
use sqisign_mp::mp_is_one;

proptest! {
    // The canonical encoding of 1 (low limb 1, rest zero) is one.
    #[test]
    fn canonical_one_is_one(n in 1usize..256) {
        let mut x = vec![0u64; n];
        x[0] = 1;
        prop_assert!(mp_is_one(&x));
    }

    // x[0] != 1 is never one, whatever the rest.
    #[test]
    fn wrong_low_limb_is_not_one(mut x in proptest::collection::vec(any::<u64>(), 1..64),
                                 lo in any::<u64>()) {
        prop_assume!(lo != 1);
        x[0] = lo;
        prop_assert!(!mp_is_one(&x));
    }

    // x[0] == 1 but some higher limb nonzero is never one.
    #[test]
    fn high_limb_set_is_not_one(mut x in proptest::collection::vec(any::<u64>(), 2..64),
                                idx in 1usize..64, val in 1u64..) {
        x[0] = 1;
        let i = 1 + (idx % (x.len() - 1));
        x[i] = val;
        prop_assert!(!mp_is_one(&x));
    }

    // Matches the direct predicate over arbitrary input.
    #[test]
    fn matches_direct(x in proptest::collection::vec(any::<u64>(), 1..256)) {
        prop_assert_eq!(
            mp_is_one(&x),
            x[0] == 1 && x[1..].iter().all(|&t| t == 0)
        );
    }
}
