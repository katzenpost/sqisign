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
//! - [`mp_mul`] is `mp_mul(c, a, b, nwords)`: the low-half multiprecision
//!   product. **It faithfully reproduces an upstream defect**: for
//!   `nwords == 1` the reference double-counts column 0 and yields
//!   `2*(a*b) mod 2^64` rather than `a*b mod 2^64`. The port mirrors the
//!   reference's algorithm so this falls out by the same logic; per the
//!   plan the C reference is the oracle and divergence is never silently
//!   introduced. Upstream is being notified by a separate correction PR.
//!   See [`mp_mul`]'s own documentation.
//! - [`mp_mul2`] is `mp_mul2(c, a, b)`, the fixed two-digit multiply.
//!   **Another faithful reproduction**: the reference omits the `a1*b0`
//!   cross term, computing `a*b - (a1*b0)*2^64` rather than the full
//!   product. All 2296 vectors pin that identity. See its documentation.
//! - [`mp_mod_2exp`] is `mp_mod_2exp(a, e, nwords)`: in-place
//!   `a = a mod 2^e` (limb mask plus zero-fill, no-op when `e` covers
//!   the full width). Correct; 1231 vectors all satisfy `a mod 2^e`.
//! - [`mp_copy`] is `mp_copy(b, a, nwords)`: a plain limb-for-limb copy
//!   (`b == a`). No quirk; 1021 vectors.
//! - [`mp_compare`] is `mp_compare(a, b, nwords)`: three-way unsigned
//!   comparison (`1`/`0`/`-1`). No quirk; 1217 vectors == `sign(a-b)`.
//! - [`mp_is_zero`] is `mp_is_zero(a, nwords)`: all-limbs-zero predicate.
//!   No quirk; 1210 vectors.
//! - [`mp_neg`] is `mp_neg(a, nwords)`. **Third faithful reproduction**:
//!   the reference adds the two's-complement `+1` to limb 0 only with no
//!   carry propagation, so it equals `-a` iff `a[0] != 0`. 1042 vectors,
//!   339 exhibiting the `a[0] == 0` quirk. See its documentation.
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

/// Low-half multiprecision product: `c = (a * b)` keeping only the low
/// `n` limbs (`n = a.len()`), the high half discarded, mirroring the
/// reference's `void mp_mul(digit_t *c, const digit_t *a, const digit_t
/// *b, size_t nwords)` ("explicitly does not use the higher half of c, as
/// we do not need in our applications").
///
/// # Faithful reproduction of an upstream defect (`nwords == 1`)
///
/// For `nwords >= 2` this is exactly `(a*b) mod 2^(64n)`. For
/// `nwords == 1` the reference is **wrong**: its per-row code does the
/// first column with `MUL(t, a[i], b[0])` and then *unconditionally*
/// repeats the last column `j = nwords-1`; when `nwords == 1` those are
/// the same column 0, so `a[0]*b[0]` is added to the row twice and the
/// result is `2*(a*b) mod 2^64`, not `a*b mod 2^64`. (The reference also
/// writes two limbs into a one-limb VLA there, a separate latent
/// out-of-bounds store; it does not affect the truncated output and is
/// not reproduced.)
///
/// Per the plan the vendored C reference is the oracle and any divergence
/// is never silently introduced, so this port mirrors the reference's
/// *algorithm* (a scratch of `n + 1` limbs gives the two-limb `MUL`
/// writes a real slot, so the doubling arises from the same control flow,
/// not from undefined behaviour). All 1048 committed C-derived vectors,
/// including the 18 single-limb cases that exhibit the doubling, pass
/// bit-for-bit. `mp_mul` has no callers in the reference and SQIsign uses
/// widths in `{4,5,8,9}`, never 1, so the defect is latent there; a
/// correction PR has been opened upstream.
///
/// # Panics
/// If `a`, `b` and `c` do not all share one length, or it is zero.
pub fn mp_mul(c: &mut [u64], a: &[u64], b: &[u64]) {
    let n = a.len();
    assert!(
        n != 0 && b.len() == n && c.len() == n,
        "mp_mul: a, b, c must share one non-zero limb count (nwords)"
    );

    let mut cc = vec![0u64; n];
    // The reference's `digit_t t[nwords]`. MUL writes two limbs, so for
    // n == 1 the reference stores out of bounds; we size t at n + 1 to
    // give that store a real slot. The extra limb is never folded into
    // the truncated result, so the logic stays faithful.
    let mut t = vec![0u64; n + 1];

    for i in 0..n {
        // MUL(t, a[i], b[0]): t[0] = lo, t[1] = hi.
        let p = (a[i] as u128) * (b[0] as u128);
        t[0] = p as u64;
        t[1] = (p >> 64) as u64;

        // for j in 1 ..= nwords-2
        let mut j = 1usize;
        while j + 1 < n {
            let uv = (a[i] as u128) * (b[j] as u128);
            let uv0 = uv as u64;
            let uv1 = (uv >> 64) as u64;
            let (s, carry) = t[j].overflowing_add(uv0);
            t[j] = s;
            // hi(a*b) <= 2^64-2, so + carry (<=1) never wraps; wrapping
            // add documents the reference's `UV[1] + carry` intent.
            t[j + 1] = uv1.wrapping_add(carry as u64);
            j += 1;
        }

        // The unconditional last column j = nwords-1. For n == 1 this is
        // column 0 again: the defect that doubles the single-limb result.
        let jl = n - 1;
        let uv0 = ((a[i] as u128) * (b[jl] as u128)) as u64;
        t[jl] = t[jl].wrapping_add(uv0); // ADDC carry-out discarded

        // mp_add(&cc[i], &cc[i], t, nwords - i): add the low n-i limbs of
        // the row into cc starting at i, the final carry discarded.
        let len = n - i;
        let mut carry = 0u64;
        for k in 0..len {
            let (s1, c1) = cc[i + k].overflowing_add(t[k]);
            let (s2, c2) = s1.overflowing_add(carry);
            cc[i + k] = s2;
            carry = (c1 as u64) | (c2 as u64);
        }
    }

    c.copy_from_slice(&cc);
}

/// Fixed two-digit operand multiply: `a[0..2] * b[0..2]` into `c[0..4]`,
/// mirroring the reference's `void mp_mul2(digit_t *c, const digit_t *a,
/// const digit_t *b)`.
///
/// # Faithful reproduction of a partial-product reference
///
/// Despite the "multiplication" name and four-digit output, the reference
/// is **not** the full 2x2 product: its body forms `a0*b0`, `a0*b1` and
/// `a1*b1` but never the `a1*b0` cross term, so it computes
/// `c = a*b - (a1*b0) * 2^64`. Every one of the 2296 committed C-derived
/// vectors satisfies exactly that identity (and the 396 that also equal
/// the true product are precisely those with `a1 == 0` or `b0 == 0`).
///
/// Whether that omission is intentional specialisation or a defect is not
/// for this port to adjudicate: per the plan the C reference is the
/// oracle and divergence is never silently introduced. The port is a
/// structural transcription of the reference's `MUL`/`ADDC` sequence, so
/// the omission arises from the same code, and all vectors pass
/// bit-for-bit. (`mp_mul2`'s callers and status are noted for sir; a
/// separate upstream question may follow.)
///
/// # Panics
/// If `a` or `b` is not length 2, or `c` is not length 4.
pub fn mp_mul2(c: &mut [u64], a: &[u64], b: &[u64]) {
    assert!(
        a.len() == 2 && b.len() == 2 && c.len() == 4,
        "mp_mul2: a, b must be 2 limbs and c must be 4 (fixed shape)"
    );

    // 64x64 -> 128 split, the reference's MUL.
    fn mul(x: u64, y: u64) -> (u64, u64) {
        let p = (x as u128) * (y as u128);
        (p as u64, (p >> 64) as u64)
    }
    // Exactly the reference's ADDC macro: tempReg = addend1 + carryIn;
    // sumOut = addend2 + tempReg; carryOut = (tempReg < carryIn) |
    // (sumOut < tempReg).
    fn addc(addend1: u64, addend2: u64, cin: u64) -> (u64, u64) {
        let temp = addend1.wrapping_add(cin);
        let sum = addend2.wrapping_add(temp);
        let cout = ((temp < cin) as u64) | ((sum < temp) as u64);
        (sum, cout)
    }

    let (t0_0, t0_1) = mul(a[0], b[0]); // MUL(t0, a[0], b[0])
    let (t1_0, t1_1) = mul(a[0], b[1]); // MUL(t1, a[0], b[1])

    let (t0_1, carry) = addc(t0_1, t1_0, 0); // ADDC(t0[1],c,t0[1],t1[0],c)
    let (t1_1, carry) = addc(0, t1_1, carry); // ADDC(t1[1],c,0,t1[1],c)

    let (t2_0, t2_1) = mul(a[1], b[1]); // MUL(t2, a[1], b[1])

    let (t2_0, carry) = addc(t2_0, t1_1, carry); // ADDC(t2[0],c,t2[0],t1[1],c)
    let (t2_1, _carry) = addc(0, t2_1, carry); // ADDC(t2[1],c,0,t2[1],c)

    c[0] = t0_0;
    c[1] = t0_1;
    c[2] = t2_0;
    c[3] = t2_1;
}

/// In-place reduction modulo `2^e`: `a = a mod 2^e`, mirroring the
/// reference's `void mp_mod_2exp(digit_t *a, unsigned int e, unsigned int
/// nwords)`.
///
/// With `q = e / 64` and `r = e % 64`: if `q < a.len()` the limb at `q`
/// is masked to its low `r` bits and every higher limb is zeroed; if
/// `q >= a.len()` (i.e. `e` covers or exceeds the full width) it is a
/// no-op, exactly as the reference. `r == 0` masks limb `q` to zero,
/// which is correct for `e` a multiple of 64. No defect here: all 1231
/// committed vectors satisfy `r == a mod 2^e`.
pub fn mp_mod_2exp(a: &mut [u64], e: u32) {
    let q = (e >> 6) as usize;
    let r = e & 63;
    if q < a.len() {
        a[q] &= (1u64 << r) - 1; // r in 0..=63; r == 0 -> mask 0
        for limb in a.iter_mut().skip(q + 1) {
            *limb = 0;
        }
    }
}

/// In-place negation, mirroring the reference's `void mp_neg(digit_t *a,
/// unsigned int nwords)`.
///
/// # Faithful reproduction of a missing carry
///
/// True two's-complement negation is `~a + 1` with the `+ 1` carried
/// through all limbs. The reference complements every limb and then adds
/// `1` to **limb 0 only, with no carry propagation**:
///
/// ```c
/// for (i = 0; i < nwords; i++) a[i] ^= -1;
/// a[0] += 1;
/// ```
///
/// So it equals `(-a) mod 2^(64n)` exactly when `a[0] != 0` (then
/// `~a[0] != MAX`, the `+1` does not overflow, and there is no carry to
/// lose). When `a[0] == 0`, `~a[0] == MAX`, the `+1` wraps limb 0 to `0`
/// and the carry that true negation would propagate is dropped, so the
/// result differs from `-a` (e.g. `mp_neg([0,0]) == [0, MAX]`, not
/// `[0,0]`). All 1042 committed vectors match this exactly; 339 of them
/// exhibit the `a[0] == 0` quirk.
///
/// Per the plan the C reference is the oracle and divergence is never
/// silently introduced; this is a structural transcription, so the
/// missing carry arises from the same code. `mp_neg` has no callers in
/// the reference; flagged for sir alongside the `mp_mul`/`mp_mul2`
/// findings.
pub fn mp_neg(a: &mut [u64]) {
    for limb in a.iter_mut() {
        *limb = !*limb;
    }
    if let Some(lo) = a.first_mut() {
        *lo = lo.wrapping_add(1); // limb 0 only; no carry, exactly as C
    }
}

/// Copy `a` into `b`, limb for limb, mirroring the reference's
/// `void mp_copy(digit_t *b, const digit_t *a, size_t nwords)`. A plain
/// copy with no quirk; all 1021 committed vectors have `b == a`.
///
/// # Panics
/// If `a` and `b` do not share one length (the reference's implicit
/// `nwords` contract, made explicit and checked here).
pub fn mp_copy(b: &mut [u64], a: &[u64]) {
    assert!(
        a.len() == b.len(),
        "mp_copy: a and b must share one limb count (nwords)"
    );
    b.copy_from_slice(a);
}

/// Three-way unsigned comparison of two equal-length little-endian
/// values: `1` if `a > b`, `-1` if `a < b`, `0` if equal. Mirrors the
/// reference's `int mp_compare(const digit_t *a, const digit_t *b,
/// unsigned int nwords)`, scanning from the most-significant limb down.
/// No quirk; all 1217 vectors equal `sign(a - b)`. Not constant-time
/// (the reference's early return is data-dependent; constant-time
/// hardening is a separate, later phase per the plan).
///
/// # Panics
/// If `a` and `b` do not share one length.
pub fn mp_compare(a: &[u64], b: &[u64]) -> i32 {
    assert!(
        a.len() == b.len(),
        "mp_compare: a and b must share one limb count (nwords)"
    );
    for i in (0..a.len()).rev() {
        if a[i] > b[i] {
            return 1;
        } else if a[i] < b[i] {
            return -1;
        }
    }
    0
}

/// Whether every limb of `a` is zero, mirroring the reference's
/// `bool mp_is_zero(const digit_t *a, unsigned int nwords)` (an OR
/// reduction over the limbs). No quirk; all 1210 vectors equal
/// "all limbs zero". The reference's reduction is constant-time; this
/// port matches its *result* (constant-time hardening is a separate
/// later phase per the plan's non-goals, and `iter().all` early-exits).
pub fn mp_is_zero(a: &[u64]) -> bool {
    a.iter().all(|&x| x == 0)
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
    fn is_zero_iff_all_limbs_zero() {
        assert!(mp_is_zero(&[0, 0, 0]));
        assert!(mp_is_zero(&[]));
        assert!(!mp_is_zero(&[0, 0, 1]));
        assert!(!mp_is_zero(&[1, 0, 0]));
        assert!(!mp_is_zero(&[u64::MAX]));
    }

    #[test]
    fn compare_is_three_way_unsigned() {
        assert_eq!(mp_compare(&[5, 0], &[5, 0]), 0);
        assert_eq!(mp_compare(&[5, 0], &[4, 0]), 1);
        assert_eq!(mp_compare(&[4, 0], &[5, 0]), -1);
        // High limb decides regardless of low limb.
        assert_eq!(mp_compare(&[u64::MAX, 1], &[0, 2]), -1);
        assert_eq!(mp_compare(&[0, 2], &[u64::MAX, 1]), 1);
        // Equal high, low decides.
        assert_eq!(mp_compare(&[7, 9], &[8, 9]), -1);
        // Empty compares equal (loop body never runs).
        assert_eq!(mp_compare(&[], &[]), 0);
    }

    #[test]
    fn copy_is_a_plain_copy() {
        let a = [1u64, 2, 3, u64::MAX];
        let mut b = [0xdeadu64; 4];
        mp_copy(&mut b, &a);
        assert_eq!(b, a);
        // Empty is a no-op (not a panic): equal (zero) lengths.
        let a: [u64; 0] = [];
        let mut b: [u64; 0] = [];
        mp_copy(&mut b, &a);
        assert_eq!(b, a);
    }

    #[test]
    fn neg_is_true_negation_iff_low_limb_nonzero() {
        // a[0] != 0: exact two's-complement negation.
        let mut a = [5u64, 0];
        mp_neg(&mut a);
        // -(5) mod 2^128 = 2^128 - 5 = [MAX-4, MAX]
        assert_eq!(a, [u64::MAX - 4, u64::MAX]);
        // a[0] == 0: the dropped carry quirk. ~[0,0]=[MAX,MAX];
        // limb0 +1 wraps to 0, NO carry -> [0, MAX], not [0,0].
        let mut a = [0u64, 0];
        mp_neg(&mut a);
        assert_eq!(a, [0, u64::MAX]);
        // a[0] == 0, high nonzero: [0, 7] -> ~ = [MAX, MAX-7];
        // limb0 +1 -> 0 (no carry) -> [0, MAX-7].
        let mut a = [0u64, 7];
        mp_neg(&mut a);
        assert_eq!(a, [0, u64::MAX - 7]);
        // Single limb with a[0] != 0 is ordinary wrapping negation.
        let mut a = [12345u64];
        mp_neg(&mut a);
        assert_eq!(a, [0u64.wrapping_sub(12345)]);
    }

    #[test]
    fn mod_2exp_masks_and_zero_fills() {
        // e = 0: everything cleared.
        let mut a = [u64::MAX, u64::MAX];
        mp_mod_2exp(&mut a, 0);
        assert_eq!(a, [0, 0]);
        // e = 64: keep limb 0, zero the rest.
        let mut a = [0xdead_beef_dead_beef, u64::MAX];
        mp_mod_2exp(&mut a, 64);
        assert_eq!(a, [0xdead_beef_dead_beef, 0]);
        // e = 68: limb1 masked to low 4 bits, limb0 intact.
        let mut a = [u64::MAX, u64::MAX, u64::MAX];
        mp_mod_2exp(&mut a, 68);
        assert_eq!(a, [u64::MAX, 0xf, 0]);
        // e at/over the full width: no-op.
        let mut a = [1u64, 2, 3];
        mp_mod_2exp(&mut a, 192);
        assert_eq!(a, [1, 2, 3]);
        let mut a = [1u64, 2, 3];
        mp_mod_2exp(&mut a, 9999);
        assert_eq!(a, [1, 2, 3]);
    }

    #[test]
    fn mul2_omits_a1_b0_cross_term() {
        // The reference computes a*b - (a1*b0)*2^64, dropping that one
        // cross term. Verify against an explicit u256 product minus it.
        let a = [3u64, 5u64];
        let b = [7u64, 11u64];
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &a, &b);
        let mut expect = mul_u256(&a, &b);
        sub_at_limb(&mut expect, 1, (a[1] as u128) * (b[0] as u128));
        assert_eq!(c, expect, "mp_mul2 must equal a*b - a1*b0*2^64");
    }

    // Helpers for the 4-limb (256-bit) reference arithmetic in tests.
    fn mul_u256(a: &[u64; 2], b: &[u64; 2]) -> [u64; 4] {
        let mut r = [0u128; 5];
        for (i, &ai) in a.iter().enumerate() {
            for (j, &bj) in b.iter().enumerate() {
                let p = (ai as u128) * (bj as u128);
                r[i + j] += p & 0xffff_ffff_ffff_ffff;
                r[i + j + 1] += p >> 64;
            }
        }
        let mut out = [0u64; 4];
        let mut carry = 0u128;
        for k in 0..4 {
            let v = r[k] + carry;
            out[k] = v as u64;
            carry = v >> 64;
        }
        out
    }
    fn sub_at_limb(x: &mut [u64; 4], limb: usize, mut amount: u128) {
        let mut k = limb;
        while amount != 0 && k < 4 {
            let cur = x[k] as u128;
            let sub = amount & 0xffff_ffff_ffff_ffff;
            if cur >= sub {
                x[k] = (cur - sub) as u64;
                amount >>= 64;
            } else {
                x[k] = (cur + (1u128 << 64) - sub) as u64;
                amount = (amount >> 64) + 1; // borrow
            }
            k += 1;
        }
    }

    #[test]
    fn mul2_equals_full_product_when_a1_or_b0_zero() {
        // a1 == 0: a1*b0 == 0, so the omission vanishes -> full product.
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &[0x1234_5678, 0], &[0x9abc_def0, 0x11]);
        assert_eq!(c, mul_u256(&[0x1234_5678, 0], &[0x9abc_def0, 0x11]));
        // b0 == 0 likewise.
        let mut c = [0u64; 4];
        mp_mul2(&mut c, &[7, 9], &[0, 13]);
        assert_eq!(c, mul_u256(&[7, 9], &[0, 13]));
    }

    #[test]
    fn mul_single_limb_reproduces_upstream_doubling() {
        // The reference double-counts column 0 at nwords==1: c = 2*(a*b).
        let mut c = [0u64; 1];
        mp_mul(&mut c, &[7], &[9]);
        assert_eq!(c, [2u64.wrapping_mul(7 * 9)]); // 126, not 63
                                                   // (2^64-1)^2 mod 2^64 = 1, doubled = 2.
        let mut c = [0u64; 1];
        mp_mul(&mut c, &[u64::MAX], &[u64::MAX]);
        assert_eq!(c, [2]);
        // Zero product is unaffected by doubling.
        let mut c = [0u64; 1];
        mp_mul(&mut c, &[0], &[12345]);
        assert_eq!(c, [0]);
    }

    #[test]
    fn mul_multilimb_is_true_low_half() {
        // nwords>=2 is correct: 2^64 * 1 in two limbs = [0, ...]; use a
        // known small case. [3,0] * [5,0] = [15,0].
        let mut c = [0u64; 2];
        mp_mul(&mut c, &[3, 0], &[5, 0]);
        assert_eq!(c, [15, 0]);
        // Carry into the high limb: [2^64-1,0] * [2,0] = lo=2^64-2,
        // hi=1.
        let mut c = [0u64; 2];
        mp_mul(&mut c, &[u64::MAX, 0], &[2, 0]);
        assert_eq!(c, [u64::MAX - 1, 1]);
        // Cross term: [a0,a1]*[b0,b1] low half =
        // [lo(a0 b0), hi(a0 b0)+lo(a0 b1)+lo(a1 b0)].
        let a0 = 0x1111_1111_1111_1111u64;
        let b0 = 0x0000_0000_0000_0010u64;
        let mut c = [0u64; 2];
        mp_mul(&mut c, &[a0, 1], &[b0, 1]);
        let full = (a0 as u128) * (b0 as u128) + (((a0 as u128) + (b0 as u128)) << 64);
        assert_eq!(c[0], full as u64);
        assert_eq!(c[1], (full >> 64) as u64);
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
