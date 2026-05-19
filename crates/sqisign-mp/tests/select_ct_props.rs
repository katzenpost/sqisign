//! Property tests for `select_ct`: branchless bitwise conditional select.

use proptest::prelude::*;
use sqisign_mp::select_ct;

fn sel(a: &[u64], b: &[u64], m: u64) -> Vec<u64> {
    let mut c = vec![0u64; a.len()];
    select_ct(&mut c, a, b, m);
    c
}

proptest! {
    // mask 0 selects a exactly.
    #[test]
    fn mask_zero_selects_a(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        prop_assert_eq!(sel(&a, &b, 0), a.clone());
    }

    // mask all-ones selects b exactly.
    #[test]
    fn mask_ones_selects_b(v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40)) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        prop_assert_eq!(sel(&a, &b, u64::MAX), b.clone());
    }

    // Arbitrary mask is the per-bit blend (a where 0, b where 1).
    #[test]
    fn arbitrary_mask_is_bit_blend(
        v in proptest::collection::vec((any::<u64>(), any::<u64>()), 1..40),
        m in any::<u64>(),
    ) {
        let a: Vec<u64> = v.iter().map(|&(x, _)| x).collect();
        let b: Vec<u64> = v.iter().map(|&(_, y)| y).collect();
        let c = sel(&a, &b, m);
        for i in 0..a.len() {
            prop_assert_eq!(c[i], (a[i] & !m) | (b[i] & m));
        }
    }

    // Selecting between equal operands is a no-op for any mask.
    #[test]
    fn equal_operands_invariant(a in proptest::collection::vec(any::<u64>(), 1..40), m in any::<u64>()) {
        prop_assert_eq!(sel(&a, &a, m), a.clone());
    }
}
