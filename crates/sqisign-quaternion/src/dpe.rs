//! Double-precision-exponent arithmetic (`dpe.h`).
//!
//! Mirrors the header-only INRIA/LORIA library
//! `vendor/the-sqisign/src/quaternion/ref/generic/internal_quaternion_headers/dpe.h`
//! at the precision required by the L2 LLL reduction. A `Dpe` value is a
//! pair `(mantissa: f64, exp: i32)` representing the real number
//! `mantissa * 2^exp` with the invariant, post-normalisation, that the
//! mantissa lies in `[1/2, 1)` (or is exactly `0.0`, in which case `exp` is
//! `i32::MIN`).
//!
//! The C reference compiles in `DPE_USE_DOUBLE` mode and consumes the
//! header through `vendor/the-sqisign/src/quaternion/ref/generic/lll/l2.c`
//! and friends. We mirror that mode only: long-double and quad variants are
//! not in scope.
//!
//! ## Bit-exact intent
//!
//! All arithmetic is expressed in IEEE 754 `f64` ops in round-to-nearest
//! mode (the only mode the C reference uses), so the per-step result is
//! reproducible across host CPUs running the standard rounding mode. The
//! single externally-visible non-bit-exact spot is the integer-to-dpe
//! conversion [`dpe_set_z`], which mirrors the **mini-gmp** quirk in
//! `mpz_get_d_2exp` (see `vendor/the-sqisign/src/mini-gmp/mini-gmp-extra.c`)
//! rather than upstream GMP: it returns the bitsize of the integer as the
//! exponent and a mantissa from `frexp(mpz_get_d(top-shifted))`. Both the C
//! cdump harness and the Rust port use this exact convention so the
//! generated vectors are reproducible from either side.

use num_bigint::Sign;

use crate::ibz::Ibz;

/// `dpe_t`: a mantissa-exponent pair representing `mant * 2^exp`.
///
/// The invariant after [`Dpe::normalize`] is `|mant| in [0.5, 1.0)` or
/// `mant == 0.0 && exp == i32::MIN`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Dpe {
    pub mant: f64,
    pub exp: i32,
}

/// Mantissa precision matching `DPE_BITSIZE` in `dpe.h`. f64 has a 52-bit
/// significand plus the implicit leading bit, so 53.
pub const DPE_BITSIZE: i32 = 53;

impl Dpe {
    /// Construct the canonical zero (mantissa `0.0`, exponent `i32::MIN`).
    pub fn zero() -> Self {
        Self {
            mant: 0.0,
            exp: i32::MIN,
        }
    }

    /// `dpe_init`: in C this is a no-op (storage is on the stack); in Rust
    /// we return the canonical zero so callers can build their own scratch
    /// arrays mechanically.
    pub fn init() -> Self {
        Self::zero()
    }

    /// `dpe_clear`: in C this is a no-op; we keep it for shape parity.
    pub fn clear(&mut self) {
        // intentionally empty
    }

    /// `dpe_normalize(x)`: rewrite `x` so that the mantissa lies in
    /// `[0.5, 1.0)` (with `exp = i32::MIN` reserved for `0.0`). Equivalent
    /// to `frexp` of the mantissa, folding the returned exponent into
    /// `self.exp`.
    pub fn normalize(&mut self) {
        if self.mant == 0.0 || !self.mant.is_finite() {
            if self.mant == 0.0 {
                self.exp = i32::MIN;
            }
            // NaN/Inf: leave exp unchanged.
            return;
        }
        let (m, e) = frexp(self.mant);
        self.mant = m;
        self.exp = self.exp.saturating_add(e);
    }

    /// `DPE_SIGN(x)`: -1, 0, or 1.
    pub fn sign(&self) -> i32 {
        if self.mant < 0.0 {
            -1
        } else if self.mant > 0.0 {
            1
        } else {
            0
        }
    }
}

/// `frexp(d) -> (m, e)` such that `d = m * 2^e` with `|m| in [0.5, 1)` (or
/// `m == 0` and `e == 0`). Pure-Rust transcription of the IEEE 754 bit
/// manipulation; matches the libm semantics the C reference uses.
fn frexp(d: f64) -> (f64, i32) {
    if d == 0.0 || !d.is_finite() {
        return (d, 0);
    }
    let bits = d.to_bits();
    let exp_field = ((bits >> 52) & 0x7FF) as i32;
    if exp_field == 0 {
        // Subnormal: normalise by multiplying up.
        let (m, e) = frexp(d * f64::from_bits(0x4350000000000000)); // 2^54
        return (m, e - 54);
    }
    // Set the exponent field to 1022 (bias-1) so the result lies in
    // [0.5, 1) with the original sign and mantissa bits.
    let new_bits = (bits & 0x800FFFFFFFFFFFFFu64) | (1022u64 << 52);
    let m = f64::from_bits(new_bits);
    let e = exp_field - 1022;
    (m, e)
}

/// `ldexp(d, e) = d * 2^e`. Pure-Rust transcription via bit manipulation.
fn ldexp(d: f64, e: i32) -> f64 {
    if d == 0.0 || !d.is_finite() || e == 0 {
        return d;
    }
    // For sub/overflow safety, scale in steps when |e| is huge. The L2
    // routine never moves beyond a few hundred at a time, so this is a
    // safety net rather than a hot path.
    let mut result = d;
    let mut remaining = e;
    while remaining > 1023 {
        result *= f64::from_bits(0x7FE0000000000000); // 2^1023
        remaining -= 1023;
    }
    while remaining < -1022 {
        result *= f64::from_bits(0x0010000000000000); // 2^-1022
        remaining += 1022;
    }
    if remaining != 0 {
        let biased = (1023i32 + remaining) as u64;
        let scale = f64::from_bits(biased << 52);
        result *= scale;
    }
    result
}

/// `dpe_scale(d, s)`: for `-DPE_BITSIZE < s <= 0` and `d in [1/2, 1)`,
/// returns `d * 2^s`. Matches the C `dpe_scale_tab` table lookup with the
/// same numeric values.
fn dpe_scale(d: f64, s: i32) -> f64 {
    // The C code asserts -DPE_BITSIZE < s <= 0; the table covers indices
    // 0..=53. We use a direct multiply by 2^s via ldexp, which is
    // bit-equivalent on IEEE 754.
    ldexp(d, s)
}

/// `dpe_set(x, y)`: `x <- y`.
pub fn dpe_set(x: &mut Dpe, y: &Dpe) {
    *x = *y;
}

/// `dpe_neg(x, y)`: `x <- -y`.
pub fn dpe_neg(x: &mut Dpe, y: &Dpe) {
    x.mant = -y.mant;
    x.exp = y.exp;
}

/// `dpe_abs(x, y)`: `x <- |y|`.
pub fn dpe_abs(x: &mut Dpe, y: &Dpe) {
    x.mant = if y.mant >= 0.0 { y.mant } else { -y.mant };
    x.exp = y.exp;
}

/// `dpe_set_d(x, y)`.
pub fn dpe_set_d(x: &mut Dpe, y: f64) {
    x.mant = y;
    x.exp = 0;
    x.normalize();
}

/// `dpe_set_ui(x, y)`.
pub fn dpe_set_ui(x: &mut Dpe, y: u64) {
    x.mant = y as f64;
    x.exp = 0;
    x.normalize();
}

/// `dpe_set_si(x, y)`.
pub fn dpe_set_si(x: &mut Dpe, y: i64) {
    x.mant = y as f64;
    x.exp = 0;
    x.normalize();
}

/// `dpe_get_d(x)`: convert back to a plain `f64`.
pub fn dpe_get_d(x: &Dpe) -> f64 {
    ldexp(x.mant, x.exp)
}

/// `dpe_set_z(x, y)`: mirror of `mini-gmp`'s `mpz_get_d_2exp` semantics.
///
/// Specifically:
///  * if `y == 0`, set `x = 0` with the canonical `i32::MIN` exponent;
///  * else, the C path does:
///    ```text
///    *exp = mpz_sizeinbase(op, 2);          // i.e. the bit length
///    if (*exp > DBL_MAX_EXP) mpz_fdiv_q_2exp(tmp, tmp, *exp - DBL_MAX_EXP);
///    ret = frexp(mpz_get_d(tmp), &tmp_exp); // tmp_exp is discarded
///    ```
///    so the returned mantissa is in `[1/2, 1)` and the exponent is the
///    bit length, ignoring the residual `tmp_exp` from `frexp`. We replicate
///    that quirk here so the bit-exact differential against the C reference
///    holds.
///
/// Note that real GMP would return the IEEE-correct `(mant, exp)` pair with
/// `op = mant * 2^exp`. The mini-gmp version, and therefore the C dump,
/// can disagree from real GMP by at most one bit in the exponent. The
/// vectors are generated against mini-gmp; this routine matches mini-gmp.
pub fn dpe_set_z(x: &mut Dpe, y: &Ibz) {
    let sign = y.0.sign();
    if matches!(sign, Sign::NoSign) {
        x.mant = 0.0;
        x.exp = i32::MIN;
        return;
    }
    let bits = y.0.bits() as usize; // bit length of magnitude
    let mag = y.0.magnitude();
    // mpz_get_d_2exp equivalent: extract top 53 bits of the magnitude
    // directly into a [0.5, 1) mantissa, without ever building the
    // full-precision double (which would overflow f64 for >1024-bit
    // magnitudes). The C reference uses `mpz_get_d_2exp(&e, y)` for the
    // same purpose; mini-gmp behaves identically.
    let mant_unsigned = top_53_bits_as_normalized_f64(
        mag.iter_u64_digits().collect::<Vec<u64>>().as_slice(),
    );
    let signed = if matches!(sign, Sign::Minus) {
        -mant_unsigned
    } else {
        mant_unsigned
    };
    x.mant = signed;
    x.exp = bits as i32;
}

/// Extract the top 53 significant bits of a positive integer (given as
/// little-endian u64 limbs) and return them as an f64 in `[0.5, 1)`.
/// Matches the result of `mpz_get_d_2exp` in GMP: the magnitude is
/// scaled to land in the standard normalized range, with the implicit
/// exponent recorded separately by the caller (here, `y.0.bits()`).
fn top_53_bits_as_normalized_f64(limbs_le: &[u64]) -> f64 {
    let n = limbs_le.len();
    if n == 0 {
        return 0.0;
    }
    let high = limbs_le[n - 1];
    if high == 0 {
        return top_53_bits_as_normalized_f64(&limbs_le[..n - 1]);
    }
    // Top 53 bits in [2^52, 2^53), then divide by 2^53 to land in
    // [0.5, 1). For a multi-limb magnitude, the top bits straddle the
    // boundary between the high limb and the limb below it.
    let clz = high.leading_zeros() as u32;
    let high_bits = (64 - clz) as u32; // count of significant bits in high limb
    let target_bits: u32 = 53;
    let mantissa_u: u64;
    if high_bits >= target_bits {
        let shift = high_bits - target_bits;
        mantissa_u = high >> shift;
    } else if n >= 2 {
        let need = target_bits - high_bits;
        let lo_shift = 64 - need;
        // high<<need carries the top high_bits bits; (next>>lo_shift) brings
        // in the topmost `need` bits of the limb below.
        mantissa_u = (high << need) | (limbs_le[n - 2] >> lo_shift);
    } else {
        // Single limb with fewer than 53 significant bits.
        mantissa_u = high << (target_bits - high_bits);
    }
    // mantissa_u now sits in [2^52, 2^53). Divide by 2^53 to land in
    // [0.5, 1).
    debug_assert!(
        (1u64 << 52) <= mantissa_u && mantissa_u < (1u64 << 53),
        "mantissa out of [2^52, 2^53)"
    );
    (mantissa_u as f64) / (1u64 << 53) as f64
}

/// Compute `mpz_get_d` for a magnitude represented as little-endian
/// `u64` limbs. Returns the (non-negative) double-precision approximation.
/// Mirrors the limb-walk in `mini-gmp.c::mpz_get_d`.
///
/// Retained but unused: superseded by [`top_53_bits_as_normalized_f64`]
/// which never builds the full-precision double and therefore handles
/// magnitudes whose bit count exceeds `f64::MAX_EXP` (1024). Kept here
/// for reference parity with the C `mpz_get_d` walk.
#[allow(dead_code)]
fn top_bits_as_f64(limbs_le: &[u64]) -> f64 {
    let n = limbs_le.len();
    if n == 0 {
        return 0.0;
    }
    let high = limbs_le[n - 1];
    if high == 0 {
        // shouldn't happen for a normalised BigInt but be defensive
        return top_bits_as_f64(&limbs_le[..n - 1]);
    }
    let clz = high.leading_zeros() as i32;
    // GMP_DBL_MANT_BITS is 53 (f64 significand), GMP_LIMB_BITS is 64.
    let mut m: i32 = clz + 53 - 64;
    let masked_high = if m < 0 {
        high & (u64::MAX << (-m))
    } else {
        high
    };
    let mut x = masked_high as f64;
    let big_b = 4.0_f64 * (1u64 << 62) as f64; // 2^64
    for i in (0..n - 1).rev() {
        x *= big_b;
        if m > 0 {
            let mut l = limbs_le[i];
            m -= 64;
            if m < 0 {
                l &= u64::MAX << (-m);
            }
            x += l as f64;
        }
    }
    x
}

/// `dpe_get_z(x, y)`: round-to-nearest integer conversion. Mirrors the
/// C `dpe_get_z`: large positives shift the mantissa up by 53 and then up
/// by `(exp - 53)` more bits; small values round to zero.
pub fn dpe_get_z(out: &mut Ibz, y: &Dpe) {
    if y.exp >= DPE_BITSIZE {
        // y is an integer. d = mantissa * 2^53, then shift by exp - 53.
        let d = y.mant * (1u64 << 53) as f64;
        // mpz_set_d on an integer-valued double is exact.
        let bi = bigint_from_double_trunc(d);
        let shift = (y.exp - DPE_BITSIZE) as u32;
        out.0 = bi << shift;
    } else if y.exp < 0 {
        // |y| < 1/2, rounds to zero.
        out.0 = num_bigint::BigInt::from(0);
    } else {
        let d = ldexp(y.mant, y.exp);
        let rounded = round_ties_to_even(d);
        out.0 = bigint_from_double_trunc(rounded);
    }
}

/// `round` (round-half-to-even is the IEEE 754 default; C's `round`
/// rounds half-away-from-zero, but the dpe code only ever rounds values
/// that came out of the size-reduction loop and treats subsequent ones
/// integrally, so either convention gives an integer at the boundary).
///
/// `round(x)` here matches C's `round()` (ties-away-from-zero) to be
/// faithful to the C reference.
fn round_ties_to_even(d: f64) -> f64 {
    // C's `round`/`__builtin_round` rounds halves away from zero, so we
    // mirror that: floor(|d| + 0.5) with the sign of d. The fplll/L2
    // algorithm tolerates either convention; matching the reference
    // simplifies bit-exact comparison.
    if d == 0.0 || !d.is_finite() {
        return d;
    }
    if d > 0.0 {
        (d + 0.5).floor()
    } else {
        -((-d + 0.5).floor())
    }
}

fn bigint_from_double_trunc(d: f64) -> num_bigint::BigInt {
    // Convert an integer-valued f64 to a BigInt. We use the bit
    // representation so that values up to 2^1023 round trip exactly.
    if d == 0.0 {
        return num_bigint::BigInt::from(0);
    }
    let negative = d < 0.0;
    let absd = if negative { -d } else { d };
    let bits = absd.to_bits();
    let exp_field = ((bits >> 52) & 0x7FF) as i32;
    let mantissa_bits = bits & 0x000FFFFFFFFFFFFF;
    if exp_field == 0 {
        // subnormal -> truncates to zero
        return num_bigint::BigInt::from(0);
    }
    let unbiased = exp_field - 1023;
    // The implicit leading 1 plus the 52-bit fraction give a 53-bit integer
    // value `mant = 1.fraction * 2^52`. The actual integer absd is
    // mant * 2^(unbiased - 52).
    let mant_int: u64 = (1u64 << 52) | mantissa_bits;
    let bi = if unbiased >= 52 {
        let shift = (unbiased - 52) as u32;
        num_bigint::BigInt::from(mant_int) << shift
    } else if unbiased < 0 {
        // |d| < 1, truncates to zero
        num_bigint::BigInt::from(0)
    } else {
        // 0 <= unbiased < 52
        num_bigint::BigInt::from(mant_int >> (52 - unbiased) as u32)
    };
    if negative {
        -bi
    } else {
        bi
    }
}

/// `dpe_add(x, y, z)`: `x <- y + z`.
pub fn dpe_add(x: &mut Dpe, y: &Dpe, z: &Dpe) {
    if y.exp > z.exp.saturating_add(DPE_BITSIZE) {
        // |z| << ulp(y), result is y.
        *x = *y;
    } else if z.exp > y.exp.saturating_add(DPE_BITSIZE) {
        *x = *z;
    } else {
        let d = y.exp - z.exp;
        if d >= 0 {
            x.mant = y.mant + dpe_scale(z.mant, -d);
            x.exp = y.exp;
        } else {
            x.mant = z.mant + dpe_scale(y.mant, d);
            x.exp = z.exp;
        }
        x.normalize();
    }
}

/// `dpe_sub(x, y, z)`: `x <- y - z`.
pub fn dpe_sub(x: &mut Dpe, y: &Dpe, z: &Dpe) {
    if y.exp > z.exp.saturating_add(DPE_BITSIZE) {
        *x = *y;
    } else if z.exp > y.exp.saturating_add(DPE_BITSIZE) {
        dpe_neg(x, z);
    } else {
        let d = y.exp - z.exp;
        if d >= 0 {
            x.mant = y.mant - dpe_scale(z.mant, -d);
            x.exp = y.exp;
        } else {
            x.mant = dpe_scale(y.mant, d) - z.mant;
            x.exp = z.exp;
        }
        x.normalize();
    }
}

/// `dpe_mul(x, y, z)`: `x <- y * z`.
pub fn dpe_mul(x: &mut Dpe, y: &Dpe, z: &Dpe) {
    x.mant = y.mant * z.mant;
    x.exp = y.exp.saturating_add(z.exp);
    x.normalize();
}

/// `dpe_div(x, y, z)`: `x <- y / z`, assuming `z != 0`.
pub fn dpe_div(x: &mut Dpe, y: &Dpe, z: &Dpe) {
    x.mant = y.mant / z.mant;
    x.exp = y.exp.saturating_sub(z.exp);
    x.normalize();
}

/// `dpe_sqrt(x, y)`: `x <- sqrt(y)`, assuming `y >= 0`.
pub fn dpe_sqrt(x: &mut Dpe, y: &Dpe) {
    if y.exp % 2 != 0 {
        x.mant = (0.5 * y.mant).sqrt();
        x.exp = (y.exp + 1) / 2;
    } else {
        x.mant = y.mant.sqrt();
        x.exp = y.exp / 2;
    }
}

/// `dpe_round(x, y)`: round to the nearest integer (ties-away-from-zero).
pub fn dpe_round(x: &mut Dpe, y: &Dpe) {
    if y.exp < 0 {
        // |y| < 1/2
        dpe_set_ui(x, 0);
    } else if y.exp >= DPE_BITSIZE {
        // already an integer
        *x = *y;
    } else {
        let d = ldexp(y.mant, y.exp);
        dpe_set_d(x, round_ties_to_even(d));
    }
}

/// `dpe_zero_p(x)`.
pub fn dpe_zero_p(x: &Dpe) -> bool {
    x.mant == 0.0
}

/// `dpe_cmp(x, y)`: positive if `x > y`, negative if `x < y`, 0 if equal.
pub fn dpe_cmp(x: &Dpe, y: &Dpe) -> i32 {
    let sx = x.sign();
    let sy = y.sign();
    let d = sx - sy;
    if d != 0 {
        return d;
    }
    if x.exp > y.exp {
        return if sx > 0 { 1 } else { -1 };
    }
    if y.exp > x.exp {
        return if sx > 0 { -1 } else { 1 };
    }
    if x.mant < y.mant {
        -1
    } else if x.mant > y.mant {
        1
    } else {
        0
    }
}

/// `dpe_cmp_d(x, d)`: compare against a plain `f64`.
pub fn dpe_cmp_d(x: &Dpe, d: f64) -> i32 {
    let mut y = Dpe::zero();
    dpe_set_d(&mut y, d);
    dpe_cmp(x, &y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn zero_normalizes() {
        let mut x = Dpe::zero();
        dpe_set_d(&mut x, 0.0);
        assert_eq!(x.mant, 0.0);
        assert_eq!(x.exp, i32::MIN);
    }

    #[test]
    fn set_d_normalizes() {
        let mut x = Dpe::zero();
        dpe_set_d(&mut x, 3.0);
        // 3.0 = 0.75 * 2^2
        assert_eq!(x.mant, 0.75);
        assert_eq!(x.exp, 2);
    }

    #[test]
    fn add_simple() {
        let (mut x, mut y, mut z) = (Dpe::zero(), Dpe::zero(), Dpe::zero());
        dpe_set_d(&mut y, 3.0);
        dpe_set_d(&mut z, 4.0);
        dpe_add(&mut x, &y, &z);
        assert_eq!(dpe_get_d(&x), 7.0);
    }

    #[test]
    fn mul_simple() {
        let (mut x, mut y, mut z) = (Dpe::zero(), Dpe::zero(), Dpe::zero());
        dpe_set_d(&mut y, 6.0);
        dpe_set_d(&mut z, 7.0);
        dpe_mul(&mut x, &y, &z);
        assert_eq!(dpe_get_d(&x), 42.0);
    }

    #[test]
    fn sub_simple() {
        let (mut x, mut y, mut z) = (Dpe::zero(), Dpe::zero(), Dpe::zero());
        dpe_set_d(&mut y, 10.0);
        dpe_set_d(&mut z, 3.0);
        dpe_sub(&mut x, &y, &z);
        assert_eq!(dpe_get_d(&x), 7.0);
    }

    #[test]
    fn div_simple() {
        let (mut x, mut y, mut z) = (Dpe::zero(), Dpe::zero(), Dpe::zero());
        dpe_set_d(&mut y, 22.0);
        dpe_set_d(&mut z, 11.0);
        dpe_div(&mut x, &y, &z);
        assert_eq!(dpe_get_d(&x), 2.0);
    }

    #[test]
    fn round_basic() {
        let (mut x, mut y) = (Dpe::zero(), Dpe::zero());
        dpe_set_d(&mut y, 2.5);
        dpe_round(&mut x, &y);
        assert_eq!(dpe_get_d(&x), 3.0);
        dpe_set_d(&mut y, -2.5);
        dpe_round(&mut x, &y);
        assert_eq!(dpe_get_d(&x), -3.0);
    }

    #[test]
    fn set_z_zero() {
        let mut x = Dpe::zero();
        dpe_set_d(&mut x, 9.9);
        dpe_set_z(&mut x, &Ibz::zero());
        assert_eq!(x.mant, 0.0);
        assert_eq!(x.exp, i32::MIN);
    }

    #[test]
    fn set_z_small() {
        let mut x = Dpe::zero();
        let y = Ibz::from(num_bigint::BigInt::from(12345));
        dpe_set_z(&mut x, &y);
        // 12345 has bit length 14, so the C-side returns exp = 14 and
        // mantissa = frexp(12345.0).0 = 12345.0 / 2^14 = 0.7535...
        // We do not require this Rust-side test to match exactly to the
        // bit; the differential vectors will confirm the wider precision.
        assert!(x.exp == 14);
        let reconstructed = dpe_get_d(&x);
        assert!((reconstructed - 12345.0).abs() < 1e-6);
    }
}
