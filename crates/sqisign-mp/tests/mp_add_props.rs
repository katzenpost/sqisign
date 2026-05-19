//! Property tests for `mp_add`.
//!
//! `mp_add` is addition in the ring Z/2^(64*n): the algebraic laws of that
//! ring are the invariants. The oracle for the value is `u128`/`BigUint`-
//! free arithmetic done limb-wise here only as cross-checks of structure,
//! not as a second implementation of the reference (the C-derived vectors
//! are the canonical correctness check).

use proptest::prelude::*;
use sqisign_mp::mp_add;

fn add(a: &[u64], b: &[u64]) -> Vec<u64> {
    let mut c = vec![0u64; a.len()];
    mp_add(&mut c, a, b);
    c
}

proptest! {
    // Commutativity: a + b == b + a.
    #[test]
    fn commutative(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        prop_assert_eq!(add(&a, &b), add(&b, &a));
    }

    // Zero is the additive identity.
    #[test]
    fn identity(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        let zero = vec![0u64; a.len()];
        prop_assert_eq!(add(&a, &zero), a.clone());
    }

    // Associativity in the ring: (a + b) + c == a + (b + c).
    #[test]
    fn associative(v in proptest::collection::vec((any::<u64>(), any::<u64>(), any::<u64>()), 1..32)) {
        let a: Vec<u64> = v.iter().map(|&(x, _, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y, _)| y).collect();
        let c: Vec<u64> = v.iter().map(|&(_, _, z)| z).collect();
        prop_assert_eq!(add(&add(&a, &b), &c), add(&a, &add(&b, &c)));
    }

    // Single-limb agrees with native wrapping u64 addition.
    #[test]
    fn single_limb_is_wrapping_add(x in any::<u64>(), y in any::<u64>()) {
        prop_assert_eq!(add(&[x], &[y]), vec![x.wrapping_add(y)]);
    }

    // The result is exactly nwords limbs: the final carry is dropped, never
    // appended (the reference's truncating contract).
    #[test]
    fn output_width_is_nwords(a in proptest::collection::vec(any::<u64>(), 1..40)) {
        let b = vec![u64::MAX; a.len()];
        prop_assert_eq!(add(&a, &b).len(), a.len());
    }
}
