//! Property tests for `mp_shiftl`.
//!
//! The defining identity, confirmed to hold across all 1168 C-derived
//! vectors, is `mp_shiftl(x, s) == (x << s) mod 2^(64*len)`. Properties
//! check it against a `u128`-built reference value for small widths and
//! check shift composition.

use proptest::prelude::*;
use sqisign_mp::{mp_add, mp_shiftl};

proptest! {
    // Two-limb shift equals the truncated 128-bit shift.
    #[test]
    fn two_limb_is_u128_shift(lo in any::<u64>(), hi in any::<u64>(), s in 1u32..=63) {
        let mut x = [lo, hi];
        let val = (lo as u128) | ((hi as u128) << 64);
        mp_shiftl(&mut x, s);
        let got = (x[0] as u128) | ((x[1] as u128) << 64);
        prop_assert_eq!(got, val.wrapping_shl(s));
    }

    // Shifting by 1 is adding the value to itself (x*2 == x+x mod 2^64n).
    #[test]
    fn shift_one_equals_self_add(v in proptest::collection::vec(any::<u64>(), 1..32)) {
        let mut shifted = v.clone();
        mp_shiftl(&mut shifted, 1);
        let mut doubled = vec![0u64; v.len()];
        mp_add(&mut doubled, &v, &v);
        prop_assert_eq!(shifted, doubled);
    }

    // Composition: shifting by a then by b equals shifting by a+b, as long
    // as a+b stays in the single-call domain (<= 63).
    #[test]
    fn shifts_compose(v in proptest::collection::vec(any::<u64>(), 1..24),
                       a in 1u32..=31, b in 1u32..=31) {
        let mut step = v.clone();
        mp_shiftl(&mut step, a);
        mp_shiftl(&mut step, b);
        let mut once = v.clone();
        mp_shiftl(&mut once, a + b);
        prop_assert_eq!(step, once);
    }

    // The low `shift` bits of the least significant limb are always zero
    // after a left shift (nothing spills into them).
    #[test]
    fn low_bits_clear(v in proptest::collection::vec(any::<u64>(), 1..32), s in 1u32..=63) {
        let mut x = v.clone();
        mp_shiftl(&mut x, s);
        prop_assert_eq!(x[0] & ((1u64 << s) - 1), 0);
    }
}
