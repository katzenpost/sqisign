//! Property tests for `fp_cswap`: branchless constant-time conditional
//! swap that consults only the LSB of `ctl`.
//!
//! Four properties:
//!
//! 1. **`ctl & 1 == 0` is a no-op.** For arbitrary five-limb `a` and
//!    `b` and any `ctl` with LSB clear, `fp_cswap(a, b, ctl)` leaves
//!    both operands unchanged, limb for limb.
//! 2. **`ctl & 1 == 1` swaps `a` and `b`.** Symmetrically, for
//!    arbitrary five-limb inputs and any `ctl` with LSB set,
//!    `fp_cswap(a, b, ctl)` leaves `a` equal to the pre-call `b` and
//!    `b` equal to the pre-call `a`, limb for limb.
//! 3. **Double swap with the same `ctl` is the identity.** Two
//!    successive `fp_cswap(a, b, ctl)` calls return both operands to
//!    their pre-call values regardless of the `ctl` value: at `ctl &
//!    1 == 0` each call is a no-op, at `ctl & 1 == 1` the second swap
//!    undoes the first.
//! 4. **Only the LSB matters.** For arbitrary five-limb inputs,
//!    `fp_cswap(a, b, 0x42)` (LSB clear, high bits set) leaves
//!    everything unchanged exactly as `fp_cswap(a, b, 0)` would, and
//!    `fp_cswap(a, b, 0x43)` (LSB set, high bits set) swaps exactly as
//!    `fp_cswap(a, b, 1)` would. This pins the reference's
//!    `(int)(ctl & 0x1)` narrowing at the call boundary.

use proptest::prelude::*;
use sqisign_gf::{fp_cswap, Fp, NWORDS_FIELD};

fn uniform5() -> impl Strategy<Value = Fp> {
    (
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
    )
        .prop_map(|(a, b, c, d, e)| [a, b, c, d, e])
}

proptest! {
    // (1) ctl with LSB clear leaves a and b unchanged.
    #[test]
    fn lsb_zero_is_noop(
        a in uniform5(),
        b in uniform5(),
        ctl in any::<u32>().prop_map(|c| c & !1u32),
    ) {
        let mut a_mut: Fp = a;
        let mut b_mut: Fp = b;
        fp_cswap(&mut a_mut, &mut b_mut, ctl);
        prop_assert_eq!(a_mut, a);
        prop_assert_eq!(b_mut, b);
        for i in 0..NWORDS_FIELD {
            prop_assert_eq!(a_mut[i], a[i]);
            prop_assert_eq!(b_mut[i], b[i]);
        }
    }

    // (2) ctl with LSB set swaps a and b.
    #[test]
    fn lsb_one_swaps(
        a in uniform5(),
        b in uniform5(),
        ctl in any::<u32>().prop_map(|c| c | 1u32),
    ) {
        let mut a_mut: Fp = a;
        let mut b_mut: Fp = b;
        fp_cswap(&mut a_mut, &mut b_mut, ctl);
        prop_assert_eq!(a_mut, b);
        prop_assert_eq!(b_mut, a);
        for i in 0..NWORDS_FIELD {
            prop_assert_eq!(a_mut[i], b[i]);
            prop_assert_eq!(b_mut[i], a[i]);
        }
    }

    // (3) Double cswap with the same ctl is the identity (involution).
    #[test]
    fn double_cswap_is_identity(
        a in uniform5(),
        b in uniform5(),
        ctl in any::<u32>(),
    ) {
        let mut a_mut: Fp = a;
        let mut b_mut: Fp = b;
        fp_cswap(&mut a_mut, &mut b_mut, ctl);
        fp_cswap(&mut a_mut, &mut b_mut, ctl);
        prop_assert_eq!(a_mut, a);
        prop_assert_eq!(b_mut, b);
    }

    // (4) Only the LSB matters: high bits are dropped.
    #[test]
    fn only_lsb_matters(
        a in uniform5(),
        b in uniform5(),
        high_bits in any::<u32>().prop_map(|c| c & !1u32),
    ) {
        // ctl_clear: LSB == 0 (high bits free), must be a no-op.
        let ctl_clear = high_bits;
        let mut a_clear: Fp = a;
        let mut b_clear: Fp = b;
        fp_cswap(&mut a_clear, &mut b_clear, ctl_clear);
        prop_assert_eq!(a_clear, a);
        prop_assert_eq!(b_clear, b);

        // ctl_set: LSB == 1 (high bits free), must swap.
        let ctl_set = high_bits | 1u32;
        let mut a_set: Fp = a;
        let mut b_set: Fp = b;
        fp_cswap(&mut a_set, &mut b_set, ctl_set);
        prop_assert_eq!(a_set, b);
        prop_assert_eq!(b_set, a);
    }

    // (4b) Pin the two literal cases the task spec calls out: 0x42 acts
    // like 0, 0x43 acts like 1. Kept as a fixed-input property so the
    // intent is visible in the suite without relying on the random
    // strategy to land on those exact values.
    #[test]
    fn ctl_0x42_is_noop_and_0x43_swaps(a in uniform5(), b in uniform5()) {
        let mut a_42: Fp = a;
        let mut b_42: Fp = b;
        fp_cswap(&mut a_42, &mut b_42, 0x42);
        prop_assert_eq!(a_42, a);
        prop_assert_eq!(b_42, b);

        let mut a_43: Fp = a;
        let mut b_43: Fp = b;
        fp_cswap(&mut a_43, &mut b_43, 0x43);
        prop_assert_eq!(a_43, b);
        prop_assert_eq!(b_43, a);
    }
}
