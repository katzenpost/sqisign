//! SQIsign `mp`: multiprecision integer arithmetic.
//!
//! Mirrors `vendor/the-sqisign/src/mp/ref/generic/mp.c`. **Phase 1,
//! unit 2.** Unlike `common`, this is a genuine reimplementation, not a
//! standardized primitive wired in: the reference is a fixed-width
//! word-array library (`digit_t = uint64_t`, little-endian limbs, a
//! runtime `nwords`), with no GCD/XGCD and no arbitrary precision, so it
//! carries no external-crate semantic risk and needs no numeric
//! dependency. (The arbitrary-precision / `gcdext` concern is the
//! `quaternion` module's, Phase 2/3; see the plan's Open Question 4.)
//!
//! Ported so far:
//! - [`mp_add`] is `mp_add(c, a, b, nwords)`: little-endian multiprecision
//!   addition, the final carry discarded exactly as the reference does.
//! - [`mp_sub`] is `mp_sub(c, a, b, nwords)`: the borrow counterpart, the
//!   final borrow discarded (so `a < b` wraps mod `2^(64n)`, as in C).
//!
//! Correctness is established as for the whole port: every committed
//! C-derived vector is replayed and bit-compared (`tests/`). Equivalence
//! to the reference is proven, not presumed.
#![forbid(unsafe_code)]

/// Multiprecision addition `c = (a + b) mod 2^(64*n)`, where `n` is the
/// common limb count of the three little-endian slices.
///
/// Mirrors the reference's `void mp_add(digit_t *c, const digit_t *a,
/// const digit_t *b, unsigned int nwords)`: it computes the truncated
/// `nwords`-limb sum and **discards the final carry** (the reference
/// returns nothing and writes no `c[nwords]`). Callers that need the carry
/// must widen the operands themselves, as the reference's callers do.
///
/// # Panics
/// If `a`, `b` and `c` do not all have the same length (the reference's
/// implicit `nwords` contract, made explicit and checked here).
pub fn mp_add(c: &mut [u64], a: &[u64], b: &[u64]) {
    assert!(
        a.len() == b.len() && a.len() == c.len(),
        "mp_add: a, b, c must share one limb count (nwords)"
    );
    let mut carry = 0u64;
    for i in 0..a.len() {
        let (s1, c1) = a[i].overflowing_add(b[i]);
        let (s2, c2) = s1.overflowing_add(carry);
        c[i] = s2;
        // c1 and c2 cannot both be set: if a[i]+b[i] overflowed, s1 is
        // small and adding carry (<=1) cannot overflow again. The OR is a
        // branchless 0/1, matching the reference's ADDC carry-out.
        carry = (c1 as u64) | (c2 as u64);
    }
}

/// Multiprecision subtraction `c = (a - b) mod 2^(64*n)`, `n` the common
/// limb count of the three little-endian slices.
///
/// Mirrors the reference's `void mp_sub(digit_t *c, const digit_t *a,
/// const digit_t *b, unsigned int nwords)`. The reference's comment says
/// "assuming a > b", but that is a caller contract, not an enforced
/// precondition: the function itself runs a borrow chain and **discards
/// the final borrow**, so `a < b` wraps modulo `2^(64n)`. The port
/// reproduces that wrapping exactly (the C-derived vectors include `a < b`
/// cases and pin it).
///
/// # Panics
/// If `a`, `b` and `c` do not all have the same length.
pub fn mp_sub(c: &mut [u64], a: &[u64], b: &[u64]) {
    assert!(
        a.len() == b.len() && a.len() == c.len(),
        "mp_sub: a, b, c must share one limb count (nwords)"
    );
    let mut borrow = 0u64;
    for i in 0..a.len() {
        let (t, b1) = a[i].overflowing_sub(b[i]);
        let (d, b2) = t.overflowing_sub(borrow);
        c[i] = d;
        // Matches the reference SUBC: borrow_out = (a<b) | (borrow_in &
        // (a-b == 0)). b1 is (a[i] < b[i]); b2 is (t < borrow), i.e.
        // t == 0 && borrow == 1. They cannot both be set.
        borrow = (b1 as u64) | (b2 as u64);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_limb_wrap_discards_carry() {
        // (2^64 - 1) + 1 = 0, carry discarded, as the reference does.
        let mut c = [0u64; 1];
        mp_add(&mut c, &[u64::MAX], &[1]);
        assert_eq!(c, [0]);
    }

    #[test]
    fn carry_ripples_across_limbs() {
        // [MAX, MAX] + [1, 0] = [0, 0] (carry out of limb 0 into limb 1,
        // then out of limb 1 and discarded).
        let mut c = [0u64; 2];
        mp_add(&mut c, &[u64::MAX, u64::MAX], &[1, 0]);
        assert_eq!(c, [0, 0]);
    }

    #[test]
    fn sub_underflow_wraps_discarding_borrow() {
        // 0 - 1 = 2^64 - 1 in one limb; [0,0] - [1,0] = [MAX, MAX]
        // (borrow ripples the whole length, final borrow discarded).
        let mut c = [0u64; 1];
        mp_sub(&mut c, &[0], &[1]);
        assert_eq!(c, [u64::MAX]);
        let mut c = [0u64; 2];
        mp_sub(&mut c, &[0, 0], &[1, 0]);
        assert_eq!(c, [u64::MAX, u64::MAX]);
    }

    #[test]
    fn add_then_sub_roundtrips() {
        let a = [0x1111_2222_3333_4444u64, 0xaaaa_bbbb_cccc_ddddu64];
        let b = [0xffff_0000_ffff_0000u64, 0x0123_4567_89ab_cdefu64];
        let mut s = [0u64; 2];
        mp_add(&mut s, &a, &b);
        let mut back = [0u64; 2];
        mp_sub(&mut back, &s, &b);
        assert_eq!(back, a);
    }
}
