//! Property tests for `multiple_mp_shiftl`.
//!
//! It is the arbitrary-amount generalisation of `mp_shiftl`; all 1280
//! C-derived vectors satisfy `r == (x << shift) mod 2^(64*len)`. These
//! pin its relationship to `mp_shiftl` and the over-width behaviour.

use proptest::prelude::*;
use sqisign_mp::{mp_shiftl, multiple_mp_shiftl};

proptest! {
    // In the single-call domain (1..=63) it equals mp_shiftl exactly.
    #[test]
    fn agrees_with_mp_shiftl_in_domain(
        v in proptest::collection::vec(any::<u64>(), 1..32),
        s in 1u32..=63,
    ) {
        let mut a = v.clone();
        let mut b = v.clone();
        multiple_mp_shiftl(&mut a, s);
        mp_shiftl(&mut b, s);
        prop_assert_eq!(a, b);
    }

    // Shifting by a then by b equals shifting once by a+b (composition,
    // now unrestricted since multiple_mp_shiftl takes any amount).
    #[test]
    fn shifts_compose(v in proptest::collection::vec(any::<u64>(), 1..24),
                       a in 1u32..400, b in 1u32..400) {
        let mut step = v.clone();
        multiple_mp_shiftl(&mut step, a);
        multiple_mp_shiftl(&mut step, b);
        let mut once = v.clone();
        multiple_mp_shiftl(&mut once, a + b);
        prop_assert_eq!(step, once);
    }

    // A shift at or beyond the full bit width zeroes the value.
    #[test]
    fn over_width_is_zero(v in proptest::collection::vec(any::<u64>(), 1..16),
                          extra in 0u32..256) {
        let full = 64 * v.len() as u32;
        let mut x = v.clone();
        multiple_mp_shiftl(&mut x, full + extra);
        prop_assert!(x.iter().all(|&w| w == 0));
    }

    // Equivalence to a u128-built reference for two limbs. Restricted to
    // s < 128: u128::wrapping_shl masks the shift by 128, so it would
    // disagree with the true over-width-is-zero behaviour for s >= 128
    // (which `over_width_is_zero` covers separately).
    #[test]
    fn two_limb_is_u128_shift(lo in any::<u64>(), hi in any::<u64>(), s in 1u32..128) {
        let mut x = [lo, hi];
        let val = (lo as u128) | ((hi as u128) << 64);
        multiple_mp_shiftl(&mut x, s);
        let got = (x[0] as u128) | ((x[1] as u128) << 64);
        prop_assert_eq!(got, val.wrapping_shl(s));
    }
}
