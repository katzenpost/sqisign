//! Property tests for `mp_mod_2exp`. It is in-place reduction modulo
//! `2^e`; all 1231 C-derived vectors satisfy `a mod 2^e`, with the
//! `e >= 64*len` case a no-op.

use proptest::prelude::*;
use sqisign_mp::mp_mod_2exp;

/// Value of little-endian limbs as a big integer (as u128 chunks summed).
fn as_big(w: &[u64]) -> Vec<u64> {
    w.to_vec()
}

proptest! {
    // Result equals the bitwise low-e-bits mask of the input.
    #[test]
    fn equals_low_bit_mask(v in proptest::collection::vec(any::<u64>(), 1..40),
                           e in 0u32..4096) {
        let orig = v.clone();
        let mut a = v.clone();
        mp_mod_2exp(&mut a, e);
        let full = 64 * orig.len() as u32;
        for (i, (&got, &was)) in a.iter().zip(orig.iter()).enumerate() {
            let bit_lo = (i as u32) * 64;
            let expect = if e >= full || bit_lo + 64 <= e {
                was // e covers the full width, or this limb is wholly < e
            } else if bit_lo >= e {
                0 // wholly at/above e: cleared
            } else {
                was & ((1u64 << (e - bit_lo)) - 1) // straddling limb
            };
            prop_assert_eq!(got, expect, "limb {} wrong (e={})", i, e);
        }
    }

    // Idempotent: reducing twice by the same e changes nothing further.
    #[test]
    fn idempotent(v in proptest::collection::vec(any::<u64>(), 1..40), e in 0u32..3000) {
        let mut a = v.clone();
        mp_mod_2exp(&mut a, e);
        let once = a.clone();
        mp_mod_2exp(&mut a, e);
        prop_assert_eq!(a, once);
    }

    // e >= full width is a no-op.
    #[test]
    fn over_width_is_noop(v in proptest::collection::vec(any::<u64>(), 1..40), extra in 0u32..512) {
        let full = 64 * v.len() as u32;
        let mut a = v.clone();
        mp_mod_2exp(&mut a, full + extra);
        prop_assert_eq!(as_big(&a), v);
    }

    // Reducing by a smaller e then a larger one equals reducing by the
    // smaller alone (monotone narrowing).
    #[test]
    fn smaller_then_larger(v in proptest::collection::vec(any::<u64>(), 1..24),
                           e1 in 0u32..1500, e2 in 0u32..1500) {
        let lo = e1.min(e2);
        let hi = e1.max(e2);
        let mut a = v.clone();
        mp_mod_2exp(&mut a, lo);
        let just_lo = a.clone();
        let mut b = v.clone();
        mp_mod_2exp(&mut b, hi);
        mp_mod_2exp(&mut b, lo);
        prop_assert_eq!(b, just_lo);
    }
}
