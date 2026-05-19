//! Property tests for `swap_ct`: branchless conditional swap.

use proptest::prelude::*;
use sqisign_mp::swap_ct;

proptest! {
    // option 0 is a no-op.
    #[test]
    fn option_zero_is_noop(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a0: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b0: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        let (mut a, mut b) = (a0.clone(), b0.clone());
        swap_ct(&mut a, &mut b, 0);
        prop_assert_eq!(a, a0);
        prop_assert_eq!(b, b0);
    }

    // option all-ones swaps a and b exactly.
    #[test]
    fn option_ones_swaps(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a0: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b0: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        let (mut a, mut b) = (a0.clone(), b0.clone());
        swap_ct(&mut a, &mut b, u64::MAX);
        prop_assert_eq!(a, b0);
        prop_assert_eq!(b, a0);
    }

    // Arbitrary option is the per-bit conditional swap.
    #[test]
    fn arbitrary_option_is_bit_swap(
        v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40),
        o in any::<u64>(),
    ) {
        let a0: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b0: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        let (mut a, mut b) = (a0.clone(), b0.clone());
        swap_ct(&mut a, &mut b, o);
        for i in 0..a0.len() {
            prop_assert_eq!(a[i], (a0[i] & !o) | (b0[i] & o));
            prop_assert_eq!(b[i], (b0[i] & !o) | (a0[i] & o));
        }
    }

    // Swapping twice with the same option returns the originals
    // (involution for any option: temp is idempotent under double xor).
    #[test]
    fn double_swap_is_identity(
        v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40),
        o in any::<u64>(),
    ) {
        let a0: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b0: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        let (mut a, mut b) = (a0.clone(), b0.clone());
        swap_ct(&mut a, &mut b, o);
        swap_ct(&mut a, &mut b, o);
        prop_assert_eq!(a, a0);
        prop_assert_eq!(b, b0);
    }
}
