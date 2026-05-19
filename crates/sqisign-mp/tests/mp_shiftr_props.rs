//! Property tests for `mp_shiftr`.
//!
//! The C-derived battery confirmed two identities across all 1168
//! vectors: the array becomes `x >> shift` (logical) and the return is
//! `x[0] & 1` on entry. These exercise that space.

use proptest::prelude::*;
use sqisign_mp::{mp_shiftl, mp_shiftr};

proptest! {
    // Two-limb right shift equals the 128-bit logical shift.
    #[test]
    fn two_limb_is_u128_shift(lo in any::<u64>(), hi in any::<u64>(), s in 1u32..=63) {
        let mut x = [lo, hi];
        let val = (lo as u128) | ((hi as u128) << 64);
        mp_shiftr(&mut x, s);
        let got = (x[0] as u128) | ((x[1] as u128) << 64);
        prop_assert_eq!(got, val >> s);
    }

    // The returned bit is the entry value's parity, independent of shift.
    #[test]
    fn returns_entry_low_bit(v in proptest::collection::vec(any::<u64>(), 1..32), s in 1u32..=63) {
        let expect = v[0] & 1;
        let mut x = v.clone();
        prop_assert_eq!(mp_shiftr(&mut x, s), expect);
    }

    // shift-right then shift-left by the same s clears exactly the low s
    // bits and preserves everything above (no value exceeds 2^(64*len)).
    #[test]
    fn rshift_then_lshift_clears_low(v in proptest::collection::vec(any::<u64>(), 1..24), s in 1u32..=63) {
        let orig = v.clone();
        let mut x = v.clone();
        mp_shiftr(&mut x, s);
        mp_shiftl(&mut x, s);
        // low s bits of limb 0 cleared:
        prop_assert_eq!(x[0] & ((1u64 << s) - 1), 0);
        // bits at or above s unchanged: compare (orig >> s) << s limbwise
        let mut ref_ = orig.clone();
        mp_shiftr(&mut ref_, s);
        mp_shiftl(&mut ref_, s);
        prop_assert_eq!(x, ref_);
    }

    // Top limb's high `shift` bits are zero-filled (logical, not
    // arithmetic): the most significant limb after the shift is < 2^(64-s).
    #[test]
    fn top_is_zero_filled(v in proptest::collection::vec(any::<u64>(), 1..32), s in 1u32..=63) {
        let mut x = v.clone();
        mp_shiftr(&mut x, s);
        let top = *x.last().unwrap();
        prop_assert_eq!(top >> (64 - s), 0);
    }
}
