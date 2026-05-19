//! Property tests for `mp_neg`.
//!
//! The reference omits carry propagation past limb 0, so it is true
//! two's-complement negation exactly when `a[0] != 0`. These pin the
//! faithful model and that characterisation.

use proptest::prelude::*;
use sqisign_mp::{mp_add, mp_neg};

proptest! {
    // The exact faithful model: complement all, +1 on limb 0 only.
    #[test]
    fn equals_faithful_model(v in proptest::collection::vec(any::<u64>(), 1..40)) {
        let mut got = v.clone();
        mp_neg(&mut got);
        let mut model: Vec<u64> = v.iter().map(|x| !x).collect();
        model[0] = model[0].wrapping_add(1);
        prop_assert_eq!(got, model);
    }

    // a[0] != 0  =>  a + mp_neg(a) == 0 (true negation in the ring).
    #[test]
    fn true_negation_when_low_nonzero(v in proptest::collection::vec(any::<u64>(), 1..40)) {
        prop_assume!(v[0] != 0);
        let mut neg = v.clone();
        mp_neg(&mut neg);
        let mut sum = vec![0u64; v.len()];
        mp_add(&mut sum, &v, &neg);
        prop_assert!(sum.iter().all(|&w| w == 0), "a + (-a) != 0");
    }

    // a[0] == 0 and some higher limb non-zero => NOT true negation
    // (the dropped carry), and specifically limb 0 of the result is 0.
    #[test]
    fn quirk_when_low_zero(mut v in proptest::collection::vec(any::<u64>(), 2..40)) {
        v[0] = 0;
        prop_assume!(v.iter().any(|&w| w != 0));
        let mut neg = v.clone();
        mp_neg(&mut neg);
        prop_assert_eq!(neg[0], 0, "~0 + 1 wraps to 0");
        let mut sum = vec![0u64; v.len()];
        mp_add(&mut sum, &v, &neg);
        prop_assert!(sum.iter().any(|&w| w != 0), "should NOT be true -a");
    }

    // Negating twice returns the original iff a[0] != 0 (true involution
    // only in the non-quirk regime).
    #[test]
    fn double_negation_in_nonquirk(v in proptest::collection::vec(any::<u64>(), 1..40)) {
        prop_assume!(v[0] != 0);
        let mut x = v.clone();
        mp_neg(&mut x);
        // After one negation a[0] = (~v[0]).wrapping_add(1); since
        // v[0] != 0, ~v[0] != MAX so x[0] != 0, the second neg is also
        // true negation, hence an involution here.
        mp_neg(&mut x);
        prop_assert_eq!(x, v);
    }
}
