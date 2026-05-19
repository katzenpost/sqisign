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
//! - [`mp_shiftl`] is `mp_shiftl(x, shift, nwords)`: in-place left shift
//!   by `1..=63` bits, truncated to `nwords` (`x <- (x << shift) mod
//!   2^(64n)`), exactly as the reference's bit-spill loop.
//! - [`mp_shiftr`] is `mp_shiftr(x, shift, nwords)`: in-place logical
//!   right shift by `1..=63` bits, returning the original `x[0] & 1`
//!   (the low bit *before* shifting), as the reference does.
//! - [`multiple_mp_shiftl`] is `multiple_mp_shiftl(x, shift, nwords)`:
//!   left shift by an *arbitrary* amount `>= 1` (the reference loops
//!   `mp_shiftl` by `RADIX-1`), i.e. `x <- (x << shift) mod 2^(64n)`
//!   for any `shift`, including amounts past the full bit width
//!   (result `0`).
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

/// In-place multiprecision left shift, `x = (x << shift) mod 2^(64*n)`,
/// where `n = x.len()` and `shift` is in `1..=63`.
///
/// Mirrors the reference's `void mp_shiftl(digit_t *x, unsigned int
/// shift, unsigned int nwords)`: each limb takes the low `shift` bits
/// spilled out of the limb below, the top of the most-significant limb is
/// discarded (no growth). The reference is only defined for
/// `1 <= shift <= RADIX-1` (`shift = 0` or `>= 64` is a shift-by-width,
/// undefined in C); this port documents and checks that domain rather
/// than silently producing a platform-dependent value.
///
/// # Panics
/// If `shift` is `0` or `>= 64` (outside the reference's defined domain),
/// or if `x` is empty.
pub fn mp_shiftl(x: &mut [u64], shift: u32) {
    assert!(
        (1..=63).contains(&shift),
        "mp_shiftl: shift must be in 1..=63 (the reference's domain)"
    );
    assert!(!x.is_empty(), "mp_shiftl: nwords must be >= 1");
    for i in (1..x.len()).rev() {
        x[i] = (x[i] << shift) ^ (x[i - 1] >> (64 - shift));
    }
    x[0] <<= shift;
}

/// In-place multiprecision logical right shift, `x = x >> shift`, where
/// `shift` is in `1..=63`. Returns the original least-significant bit
/// (`x[0] & 1`, captured *before* the shift), exactly as the reference's
/// `digit_t mp_shiftr(digit_t *x, unsigned int shift, unsigned int
/// nwords)`.
///
/// Each limb takes the high `shift` bits of the limb above; the top is
/// zero-filled (logical, not arithmetic). The returned bit is *not* the
/// last bit shifted out in general, only the value's bit 0 on entry, as
/// the reference defines it (its callers use it as an "is odd" probe).
///
/// # Panics
/// If `shift` is `0` or `>= 64`, or if `x` is empty.
pub fn mp_shiftr(x: &mut [u64], shift: u32) -> u64 {
    assert!(
        (1..=63).contains(&shift),
        "mp_shiftr: shift must be in 1..=63 (the reference's domain)"
    );
    assert!(!x.is_empty(), "mp_shiftr: nwords must be >= 1");
    let bit_out = x[0] & 1;
    let n = x.len();
    for i in 0..n - 1 {
        x[i] = (x[i] >> shift) ^ (x[i + 1] << (64 - shift));
    }
    x[n - 1] >>= shift;
    bit_out
}

/// In-place left shift by an arbitrary amount `shift >= 1`, equal to
/// `x = (x << shift) mod 2^(64*n)` for any `shift` (including amounts at
/// or beyond the full bit width, which leave `x` all zero).
///
/// Mirrors the reference's `multiple_mp_shiftl`, which composes
/// [`mp_shiftl`] in steps of `RADIX-1` (63) bits then a final step of the
/// remainder. The decomposition keeps every individual `mp_shiftl` within
/// its defined `1..=63` domain (for `shift >= 1` the remainder is always
/// in `1..=63`, never `0`), so the reference is well-defined for all
/// `shift >= 1`; the port reproduces that exactly. `shift == 0` is the
/// reference's one undefined input (it would call `mp_shiftl(x, 0)`); the
/// port rejects it rather than reproduce the C undefined behaviour.
///
/// # Panics
/// If `shift == 0`, or if `x` is empty.
pub fn multiple_mp_shiftl(x: &mut [u64], shift: u32) {
    assert!(shift >= 1, "multiple_mp_shiftl: shift must be >= 1");
    assert!(!x.is_empty(), "multiple_mp_shiftl: nwords must be >= 1");
    let mut t = shift;
    while t > 63 {
        mp_shiftl(x, 63);
        t -= 63;
    }
    mp_shiftl(x, t); // t is now in 1..=63
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
    fn shiftl_spills_across_limbs_and_truncates() {
        // [0x8000_0000_0000_0000, 0] << 1 = [0, 1] (top bit spills up).
        let mut x = [0x8000_0000_0000_0000u64, 0];
        mp_shiftl(&mut x, 1);
        assert_eq!(x, [0, 1]);
        // Single limb: top bit shifted out is discarded (truncation).
        let mut y = [0x8000_0000_0000_0000u64];
        mp_shiftl(&mut y, 1);
        assert_eq!(y, [0]);
    }

    #[test]
    fn shiftl_is_multiply_by_pow2_mod() {
        let mut x = [0x0123_4567_89ab_cdefu64, 0xfedc_ba98_7654_3210u64];
        let val = (x[0] as u128) | ((x[1] as u128) << 64);
        mp_shiftl(&mut x, 5);
        let got = (x[0] as u128) | ((x[1] as u128) << 64);
        assert_eq!(got, val.wrapping_shl(5));
    }

    #[test]
    fn shiftr_logical_and_returns_entry_low_bit() {
        // [0, 1] >> 1 = [0x8000_0000_0000_0000, 0]; returned bit is the
        // ORIGINAL x[0]&1 (0 here), not the bit that fell out.
        let mut x = [0u64, 1];
        let bit = mp_shiftr(&mut x, 1);
        assert_eq!(x, [0x8000_0000_0000_0000, 0]);
        assert_eq!(bit, 0);
        // Odd low limb: returned bit is 1, top zero-filled (logical).
        let mut y = [0xffff_ffff_ffff_ffffu64];
        let bit = mp_shiftr(&mut y, 4);
        assert_eq!(y, [0x0fff_ffff_ffff_ffff]);
        assert_eq!(bit, 1);
    }

    #[test]
    fn shiftr_then_shiftl_clears_low_bits() {
        let mut x = [0x0123_4567_89ab_cdefu64, 0xfedc_ba98_7654_3210u64];
        let orig = x;
        mp_shiftr(&mut x, 7);
        mp_shiftl(&mut x, 7);
        assert_eq!(x[0], orig[0] & !0x7f);
        assert_eq!(x[1], orig[1]);
    }

    #[test]
    fn multiple_shiftl_matches_single_in_domain_and_clears_past_width() {
        // For shift in 1..=63 it equals mp_shiftl.
        let mut a = [0x0123_4567_89ab_cdefu64, 0xfedc_ba98_7654_3210u64];
        let mut b = a;
        multiple_mp_shiftl(&mut a, 50);
        mp_shiftl(&mut b, 50);
        assert_eq!(a, b);
        // A shift at or beyond the full bit width zeroes the value.
        let mut c = [0xffff_ffff_ffff_ffffu64, 0x1];
        multiple_mp_shiftl(&mut c, 128);
        assert_eq!(c, [0, 0]);
        // A cross-RADIX shift: [1] << 64 in two limbs = [0, 1].
        let mut d = [1u64, 0];
        multiple_mp_shiftl(&mut d, 64);
        assert_eq!(d, [0, 1]);
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
