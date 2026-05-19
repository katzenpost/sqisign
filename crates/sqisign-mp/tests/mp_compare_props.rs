//! Property tests for `mp_compare`: it is a total order consistent with
//! the integer value, antisymmetric, and reflexive.

use proptest::prelude::*;
use sqisign_mp::{mp_add, mp_compare};

proptest! {
    // Reflexive: a vs a is 0.
    #[test]
    fn reflexive(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        prop_assert_eq!(mp_compare(&a, &a), 0);
    }

    // Antisymmetric: compare(a,b) == -compare(b,a).
    #[test]
    fn antisymmetric(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        prop_assert_eq!(mp_compare(&a, &b), -mp_compare(&b, &a));
    }

    // Consistent with the value: result sign equals the sign of the
    // big-integer difference (checked via the top differing limb).
    #[test]
    fn matches_value_order(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        let mut expect = 0i32;
        for i in (0..a.len()).rev() {
            if a[i] != b[i] {
                expect = if a[i] > b[i] { 1 } else { -1 };
                break;
            }
        }
        prop_assert_eq!(mp_compare(&a, &b), expect);
    }

    // Adding 1 (without overflow of the whole value) makes a strictly
    // greater: compare(a+1, a) == 1 when a is not all-ones.
    #[test]
    fn successor_is_greater(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        prop_assume!(!a.iter().all(|&w| w == u64::MAX));
        let mut one = vec![0u64; a.len()];
        one[0] = 1;
        let mut ap1 = vec![0u64; a.len()];
        mp_add(&mut ap1, &a, &one);
        prop_assert_eq!(mp_compare(&ap1, &a), 1);
    }

    // Transitivity on a sorted triple.
    #[test]
    fn transitive(x in any::<u64>(), y in any::<u64>(), z in any::<u64>()) {
        let mut v = [x, y, z];
        v.sort_unstable();
        let lo = [v[0]];
        let mid = [v[1]];
        let hi = [v[2]];
        if mp_compare(&lo, &mid) <= 0 && mp_compare(&mid, &hi) <= 0 {
            prop_assert!(mp_compare(&lo, &hi) <= 0);
        }
    }
}
