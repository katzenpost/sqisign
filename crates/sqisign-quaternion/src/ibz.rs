//! Big integer wrappers (`ibz_t`).
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/intbig.c`. The
//! C reference uses GMP's `mpz_t`; this port uses [`num_bigint::BigInt`].
//! See the crate root for the canonical-bytes boundary contract and the
//! cryptographic non-secrecy of these paths.

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// `ibz_t`: arbitrary-precision signed integer.
///
/// Newtype around [`num_bigint::BigInt`]. The wrapper exists so that we can
/// mirror the C reference's separate type identity (the same `mpz_t` is
/// used for many semantically distinct roles in the upstream code) and so
/// that we control the canonical-bytes (de)serialization in one place.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Ibz(pub BigInt);

impl Ibz {
    /// Construct from a [`BigInt`].
    pub fn new(n: BigInt) -> Self {
        Self(n)
    }

    /// Construct zero.
    pub fn zero() -> Self {
        Self(BigInt::zero())
    }

    /// Construct from a `i32`.
    pub fn from_i32(n: i32) -> Self {
        Self(BigInt::from(n))
    }

    /// Canonical bytes encoding used by the differential harness.
    ///
    /// Format: 1 byte sign tag (0x00 non-negative, 0x01 negative), then a
    /// little-endian `u32` magnitude length N, then N big-endian magnitude
    /// bytes. Zero is encoded as `[0x00, 0x00, 0x00, 0x00, 0x00]` (sign
    /// non-negative, length zero). Matches the natural shape of GMP's
    /// `mpz_export` (most-significant byte first, BIG endian) plus a sign
    /// tag carried separately because `mpz_export` ignores it.
    pub fn to_canonical_bytes(&self) -> Vec<u8> {
        let (sign, mag) = self.0.to_bytes_be();
        let sign_byte = if matches!(sign, Sign::Minus) {
            1u8
        } else {
            0u8
        };
        // Strip a single leading zero introduced by `to_bytes_be` for Sign::NoSign,
        // which already returns `vec![0]` for zero; we want zero-length here.
        let mag: &[u8] = if self.0.is_zero() { &[] } else { &mag };
        let mut out = Vec::with_capacity(1 + 4 + mag.len());
        out.push(sign_byte);
        out.extend_from_slice(&(mag.len() as u32).to_le_bytes());
        out.extend_from_slice(mag);
        out
    }

    /// Decode the canonical bytes encoding produced by [`Self::to_canonical_bytes`].
    pub fn from_canonical_bytes(bytes: &[u8]) -> Result<Self, &'static str> {
        if bytes.len() < 5 {
            return Err("ibz canonical bytes: header truncated");
        }
        let sign_byte = bytes[0];
        let len = u32::from_le_bytes([bytes[1], bytes[2], bytes[3], bytes[4]]) as usize;
        if bytes.len() != 5 + len {
            return Err("ibz canonical bytes: payload length mismatch");
        }
        let mag = &bytes[5..];
        let sign = match sign_byte {
            0 => {
                if len == 0 {
                    Sign::NoSign
                } else {
                    Sign::Plus
                }
            }
            1 => {
                if len == 0 {
                    return Err("ibz canonical bytes: negative zero is not canonical");
                }
                Sign::Minus
            }
            _ => return Err("ibz canonical bytes: invalid sign byte"),
        };
        Ok(Self(BigInt::from_bytes_be(sign, mag)))
    }
}

impl From<BigInt> for Ibz {
    fn from(n: BigInt) -> Self {
        Self(n)
    }
}

impl From<i32> for Ibz {
    fn from(n: i32) -> Self {
        Self(BigInt::from(n))
    }
}

/// `ibz_const_zero`.
pub fn ibz_const_zero() -> Ibz {
    Ibz(BigInt::zero())
}

/// `ibz_const_one`.
pub fn ibz_const_one() -> Ibz {
    Ibz(BigInt::one())
}

/// `ibz_const_two`.
pub fn ibz_const_two() -> Ibz {
    Ibz(BigInt::from(2u32))
}

/// `ibz_const_three`.
pub fn ibz_const_three() -> Ibz {
    Ibz(BigInt::from(3u32))
}

// The C reference exposes these as `extern const ibz_t`. We can't make those
// `const` in Rust (BigInt allocates), so the only public constructors are
// the snake-case functions above.

/// `ibz_add(sum, a, b)`: `sum = a + b`.
pub fn ibz_add(sum: &mut Ibz, a: &Ibz, b: &Ibz) {
    sum.0 = &a.0 + &b.0;
}

/// `ibz_sub(diff, a, b)`: `diff = a - b`.
pub fn ibz_sub(diff: &mut Ibz, a: &Ibz, b: &Ibz) {
    diff.0 = &a.0 - &b.0;
}

/// `ibz_mul(prod, a, b)`: `prod = a * b`.
pub fn ibz_mul(prod: &mut Ibz, a: &Ibz, b: &Ibz) {
    prod.0 = &a.0 * &b.0;
}

/// `ibz_neg(neg, a)`: `neg = -a`.
pub fn ibz_neg(neg: &mut Ibz, a: &Ibz) {
    neg.0 = -(&a.0);
}

/// `ibz_abs(abs, a)`: `abs = |a|`.
pub fn ibz_abs(abs: &mut Ibz, a: &Ibz) {
    abs.0 = a.0.abs();
}

/// `ibz_div(q, r, a, b)`: Euclidean division rounding the quotient
/// **towards zero** (truncated division, matching GMP's `mpz_tdiv_qr`).
///
/// Satisfies `r + q*b == a` with `0 <= |r| < |b|` and `sign(r) == sign(a)`
/// when `r != 0`. Caller must ensure `b != 0`.
pub fn ibz_div(q: &mut Ibz, r: &mut Ibz, a: &Ibz, b: &Ibz) {
    // num-bigint's Div / Rem for BigInt is truncated (round-toward-zero).
    let (qv, rv) = num_integer::Integer::div_rem(&a.0, &b.0);
    q.0 = qv;
    r.0 = rv;
}

/// `ibz_div_2exp(q, a, exp)`: `q = sign(a) * (|a| >> exp)`, i.e. truncated
/// division by `2^exp`.
pub fn ibz_div_2exp(q: &mut Ibz, a: &Ibz, exp: u32) {
    // GMP's mpz_tdiv_q_2exp shifts the magnitude and preserves sign. Match
    // by shifting the absolute value then re-applying the sign.
    let sign = a.0.sign();
    let mag = a.0.magnitude() >> exp;
    q.0 = BigInt::from_biguint(sign, mag);
}

/// `ibz_div_floor(q, r, n, d)`: floor-rounded Euclidean division. Matches
/// GMP's `mpz_fdiv_qr` (quotient rounded toward minus infinity).
pub fn ibz_div_floor(q: &mut Ibz, r: &mut Ibz, n: &Ibz, d: &Ibz) {
    let (qv, rv) = num_integer::Integer::div_mod_floor(&n.0, &d.0);
    q.0 = qv;
    r.0 = rv;
}

/// `ibz_mod(r, a, b)`: `r = a mod b`, with `r` always non-negative
/// regardless of the sign of `b`. Matches GMP's `mpz_mod`.
pub fn ibz_mod(r: &mut Ibz, a: &Ibz, b: &Ibz) {
    // GMP mpz_mod returns the non-negative remainder; the sign of `b` is
    // ignored. num-bigint provides `mod_floor` against |b| to match.
    let abs_b = b.0.abs();
    r.0 = num_integer::Integer::mod_floor(&a.0, &abs_b);
}

/// `ibz_mod_ui(n, d)`: `|n| mod d` as an unsigned long. Matches GMP's
/// `mpz_fdiv_ui` (returns the non-negative remainder when dividing by an
/// unsigned long).
///
/// Aborts if `d == 0` (the C reference relies on GMP's abort behaviour).
pub fn ibz_mod_ui(n: &Ibz, d: u64) -> u64 {
    assert!(d != 0, "ibz_mod_ui: division by zero");
    let d_big = BigInt::from(d);
    let r = num_integer::Integer::mod_floor(&n.0, &d_big);
    // r is in [0, d) and fits in u64.
    let (_, bytes) = r.to_bytes_be();
    let mut acc: u64 = 0;
    for b in bytes {
        acc = (acc << 8) | b as u64;
    }
    acc
}

/// `ibz_divides(a, b)`: 1 if `b | a` (i.e. `a` is divisible by `b`), 0
/// otherwise.
///
/// Note the C argument order matches the helper name (`a` is the
/// dividend); the divisibility *test* is `b | a`. This mirrors GMP's
/// `mpz_divisible_p(a, b)`.
pub fn ibz_divides(a: &Ibz, b: &Ibz) -> i32 {
    if b.0.is_zero() {
        // mpz_divisible_p: nonzero iff a is also zero.
        return if a.0.is_zero() { 1 } else { 0 };
    }
    if num_integer::Integer::is_multiple_of(&a.0, &b.0) {
        1
    } else {
        0
    }
}

/// `ibz_pow(pow, x, e)`: `pow = x^e`. `0^0` yields 1 (matching GMP).
pub fn ibz_pow(pow: &mut Ibz, x: &Ibz, e: u32) {
    pow.0 = x.0.pow(e);
}

/// `ibz_pow_mod(pow, x, e, m)`: `pow = (x^e) mod m`.
///
/// Matches GMP's `mpz_powm`: when `e < 0`, requires an inverse of `x` mod
/// `m` to exist; this port follows the same contract.
pub fn ibz_pow_mod(pow: &mut Ibz, x: &Ibz, e: &Ibz, m: &Ibz) {
    pow.0 = x.0.modpow(&e.0, &m.0);
}

/// `ibz_two_adic(pow)`: position of the first 1-bit in `|pow|`, i.e. the
/// 2-adic valuation. Returns -1 if `pow == 0` (matches GMP's `mpz_scan1`,
/// which returns ULONG_MAX; we return the signed sentinel for portability
/// in `i32`).
pub fn ibz_two_adic(pow: &Ibz) -> i32 {
    if pow.0.is_zero() {
        return -1;
    }
    let mag = pow.0.magnitude();
    // Scan limbs from low to high.
    let digits = mag.to_u32_digits();
    let mut bit = 0i32;
    for d in digits {
        if d == 0 {
            bit += 32;
            continue;
        }
        return bit + (d.trailing_zeros() as i32);
    }
    -1
}

/// `ibz_cmp(a, b)`: sign of `a - b`, normalized to `{-1, 0, +1}`.
pub fn ibz_cmp(a: &Ibz, b: &Ibz) -> i32 {
    match a.0.cmp(&b.0) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// `ibz_is_zero(x)`.
pub fn ibz_is_zero(x: &Ibz) -> i32 {
    if x.0.is_zero() {
        1
    } else {
        0
    }
}

/// `ibz_is_one(x)`.
pub fn ibz_is_one(x: &Ibz) -> i32 {
    if x.0.is_one() {
        1
    } else {
        0
    }
}

/// `ibz_cmp_int32(x, y)`: sign of `x - y`, normalized to `{-1, 0, +1}`.
pub fn ibz_cmp_int32(x: &Ibz, y: i32) -> i32 {
    let yb = BigInt::from(y);
    match x.0.cmp(&yb) {
        std::cmp::Ordering::Less => -1,
        std::cmp::Ordering::Equal => 0,
        std::cmp::Ordering::Greater => 1,
    }
}

/// `ibz_is_even(x)`.
pub fn ibz_is_even(x: &Ibz) -> i32 {
    if num_integer::Integer::is_even(&x.0) {
        1
    } else {
        0
    }
}

/// `ibz_is_odd(x)`.
pub fn ibz_is_odd(x: &Ibz) -> i32 {
    if num_integer::Integer::is_odd(&x.0) {
        1
    } else {
        0
    }
}

/// `ibz_set(i, x)`: set `i` to the signed `i32` `x`.
pub fn ibz_set(i: &mut Ibz, x: i32) {
    i.0 = BigInt::from(x);
}

/// `ibz_get(i)`: extract a 32-bit summary of `i` matching the C
/// reference's truncation rule.
///
/// The C reference does (on 64-bit `long`):
/// ```text
/// signed long t = mpz_get_si(i);
/// return ((int32_t)(((t >> 63) & 0x80000000) | (t & 0x7FFFFFFF));
/// ```
/// i.e. take the sign bit of the 64-bit `long` into bit 31 and the low 31
/// bits of the long verbatim. With the mini-gmp impl that the harness
/// links, `mpz_get_si(i)` returns `sign(i) * (|i| mod 2^63)` (the low 63
/// bits of the magnitude carrying the sign). The Rust port reproduces
/// this byte-for-byte.
pub fn ibz_get(i: &Ibz) -> i32 {
    // Compute mpz_get_si equivalent under mini-gmp semantics.
    let mag = i.0.magnitude();
    let low63: u64 = {
        let digits = mag.to_u64_digits();
        let low = *digits.first().unwrap_or(&0);
        low & 0x7FFF_FFFF_FFFF_FFFF
    };
    let t: i64 = if i.0.sign() == Sign::Minus {
        (low63 as i64).wrapping_neg()
    } else {
        low63 as i64
    };
    // The C reduction: top bit of t into bit 31, low 31 bits verbatim.
    let top_bit_in_31: u32 = (((t as u64) >> 63) as u32) << 31;
    let low31: u32 = (t as u64) as u32 & 0x7FFF_FFFF;
    (top_bit_in_31 | low31) as i32
}

/// `ibz_bitsize(a)`: number of bits needed to represent `|a|`, with the
/// special case `ibz_bitsize(0) == 1` (matching GMP's `mpz_sizeinbase(_, 2)`).
pub fn ibz_bitsize(a: &Ibz) -> i32 {
    if a.0.is_zero() {
        // mpz_sizeinbase(0, base) returns 1.
        1
    } else {
        a.0.bits() as i32
    }
}

/// `ibz_size_in_base(a, base)`: number of digits of `|a|` in `base`, with
/// the special case `size_in_base(0, base) == 1`.
///
/// Implemented for base 2, 10, and 16 (the bases the reference uses).
/// Other bases are not used by the upstream code paths in scope.
pub fn ibz_size_in_base(a: &Ibz, base: i32) -> i32 {
    assert!(
        base == 2 || base == 10 || base == 16,
        "ibz_size_in_base: only bases 2, 10, 16 supported (got {base})"
    );
    if a.0.is_zero() {
        return 1;
    }
    let mag = a.0.magnitude();
    match base {
        2 => mag.bits() as i32,
        16 => mag.bits().div_ceil(4) as i32,
        10 => {
            // mpz_sizeinbase(_, 10) returns the number of decimal digits or
            // one more. Match by formatting; the reference's only consumer
            // is buffer sizing for string conversion, where overshooting by
            // one is harmless. We choose the exact decimal-digit count to
            // make the boundary deterministic.
            mag.to_str_radix(10).len() as i32
        }
        _ => unreachable!(),
    }
}

/// `ibz_copy_digits(target, dig, dig_len)`: build a non-negative `Ibz`
/// from a little-endian array of `digit_t` (`u64`) words.
///
/// Matches GMP's `mpz_import(target, dig_len, -1, sizeof(digit_t), 0, 0,
/// dig)` (LE word order, native endianness within each word).
pub fn ibz_copy_digits(target: &mut Ibz, dig: &[u64]) {
    let mut bytes = Vec::with_capacity(dig.len() * 8);
    // mpz_import with order=-1 and endian=0 (native) consumes the words in
    // little-endian *word* order with each word in machine endianness.
    // We are on little-endian targets only (matching the C harness host),
    // so the byte serialization is just the LE bytes of each word in
    // order.
    for w in dig {
        bytes.extend_from_slice(&w.to_le_bytes());
    }
    // The reference treats the array as a big non-negative integer: words
    // are the magnitude in little-endian word order. Reverse to BE bytes.
    bytes.reverse();
    target.0 = BigInt::from_bytes_be(Sign::Plus, &bytes);
    // BigInt normalizes the zero case to Sign::NoSign automatically.
}

/// `ibz_to_digits(target, ibz)`: export an `Ibz` to a little-endian array
/// of `digit_t` (`u64`) words. The number must be non-negative. The
/// caller-supplied buffer must already be zeroed and long enough to hold
/// the value; this function only writes as many words as needed and
/// leaves the rest unchanged (matching the C reference which relies on
/// the `ibz_to_digit_array` macro to memset first).
pub fn ibz_to_digits(target: &mut [u64], ibz: &Ibz) {
    assert!(
        ibz.0.sign() != Sign::Minus,
        "ibz_to_digits: input must be non-negative"
    );
    // Match GMP behaviour: zero writes nothing.
    if ibz.0.is_zero() {
        // The C harness explicitly sets target[0] = 0 before calling; the
        // wrapper here mirrors that to keep the boundary deterministic for
        // the zero input.
        if let Some(w) = target.first_mut() {
            *w = 0;
        }
        return;
    }
    let digits = ibz.0.magnitude().to_u64_digits();
    assert!(
        digits.len() <= target.len(),
        "ibz_to_digits: target buffer too small ({} words for {} needed)",
        target.len(),
        digits.len()
    );
    for (slot, word) in target.iter_mut().zip(digits.iter()) {
        *slot = *word;
    }
}

/// `ibz_probab_prime(n, reps)`: probabilistic primality.
///
/// GMP returns 2 for "certainly prime", 1 for "probably prime", 0 for
/// "certainly composite". This port returns the same trichotomy.
/// Uses a Miller-Rabin test seeded deterministically to keep the
/// boundary reproducible; `reps` is interpreted as the witness count.
pub fn ibz_probab_prime(n: &Ibz, reps: i32) -> i32 {
    // Match GMP's semantics: 0 for compositeness, 2 for proven prime (we
    // can only prove primes up to a small bound here), 1 otherwise.
    if n.0.sign() != Sign::Plus {
        // Negative or zero: not prime.
        return 0;
    }
    if n.0.is_one() {
        return 0;
    }
    // Small primes up to a cap: deterministic certain answer.
    let small_primes: [u64; 25] = [
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89,
        97,
    ];
    // Use trial division by small primes to certify primality up to 97^2.
    // We must trial-divide before falling through to Miller-Rabin.
    let v_low: u64 = {
        let digs = n.0.to_u64_digits().1;
        *digs.first().unwrap_or(&0)
    };
    if n.0.bits() <= 7 && v_low <= 97 {
        return if small_primes.contains(&v_low) { 2 } else { 0 };
    }
    for p in &small_primes {
        if num_integer::Integer::is_multiple_of(&n.0, &BigInt::from(*p)) {
            return 0;
        }
    }
    // Miller-Rabin with deterministic bases. The differential boundary
    // pins specific witnesses so the count must be reproducible.
    miller_rabin(&n.0, reps.max(1) as usize)
}

fn miller_rabin(n: &BigInt, k: usize) -> i32 {
    let one: BigInt = BigInt::one();
    let two: BigInt = BigInt::from(2u32);
    let n_minus_1 = n - &one;
    let mut d = n_minus_1.clone();
    let mut r = 0u32;
    while d.is_even() {
        d >>= 1;
        r += 1;
    }
    // Deterministic witness set: the first `k` odd primes (skipping 2). We
    // already trial-divided by these, but they remain valid Miller-Rabin
    // witnesses.
    let witnesses: [u64; 13] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41];
    'outer: for &a_small in witnesses.iter().take(k.min(witnesses.len())) {
        let a = BigInt::from(a_small);
        if &a >= n {
            continue;
        }
        let mut x = a.modpow(&d, n);
        if x == one || x == n_minus_1 {
            continue;
        }
        for _ in 0..(r - 1) {
            x = x.modpow(&two, n);
            if x == n_minus_1 {
                continue 'outer;
            }
        }
        return 0;
    }
    // We cannot certify primality without a full proof; report 1 as GMP
    // would for "probably prime".
    1
}

/// `ibz_gcd(gcd, a, b)`: `gcd = gcd(a, b)` (non-negative).
pub fn ibz_gcd(gcd: &mut Ibz, a: &Ibz, b: &Ibz) {
    gcd.0 = num_integer::Integer::gcd(&a.0, &b.0);
}

/// `ibz_invmod(inv, a, mod_)`: compute `inv = a^{-1} (mod mod_)` if it
/// exists. Returns 1 on success, 0 if no inverse exists.
pub fn ibz_invmod(inv: &mut Ibz, a: &Ibz, mod_: &Ibz) -> i32 {
    use num_integer::Integer;
    let m = mod_.0.abs();
    if m.is_zero() {
        return 0;
    }
    let a_mod = a.0.mod_floor(&m);
    let ext = a_mod.extended_gcd(&m);
    if !ext.gcd.is_one() {
        return 0;
    }
    let mut x = ext.x.mod_floor(&m);
    if x.sign() == Sign::Minus {
        x += &m;
    }
    inv.0 = x;
    1
}

/// `ibz_legendre(a, p)`: Legendre symbol of `a` mod `p`. Returns -1, 0, or 1.
///
/// Assumes `p` is an odd prime (the reference asserts this in debug
/// builds).
pub fn ibz_legendre(a: &Ibz, p: &Ibz) -> i32 {
    // Reduce a mod p.
    let p_abs = p.0.abs();
    let a_mod = num_integer::Integer::mod_floor(&a.0, &p_abs);
    if a_mod.is_zero() {
        return 0;
    }
    // Euler criterion: a^((p-1)/2) mod p, lifted to {-1, 0, 1}.
    let exp = (&p_abs - 1u32) >> 1;
    let r = a_mod.modpow(&exp, &p_abs);
    if r.is_one() {
        1
    } else if r == &p_abs - 1u32 {
        -1
    } else {
        // Should not occur for prime p, but in the spirit of GMP's
        // mpz_legendre, treat any other result as 0.
        0
    }
}

/// `ibz_sqrt(sqrt, a)`: if `a` is a perfect square, set `sqrt` to its
/// non-negative square root and return 1; otherwise return 0.
pub fn ibz_sqrt(sqrt: &mut Ibz, a: &Ibz) -> i32 {
    if a.0.sign() == Sign::Minus {
        return 0;
    }
    let s = a.0.sqrt();
    if &s * &s == a.0 {
        sqrt.0 = s;
        1
    } else {
        0
    }
}

/// `ibz_sqrt_floor(sqrt, a)`: `sqrt = floor(sqrt(a))`. Caller must ensure
/// `a >= 0`.
pub fn ibz_sqrt_floor(sqrt: &mut Ibz, a: &Ibz) {
    assert!(
        a.0.sign() != Sign::Minus,
        "ibz_sqrt_floor: input must be non-negative"
    );
    sqrt.0 = a.0.sqrt();
}

/// `ibz_sqrt_mod_p(sqrt, a, p)`: square root of `a` mod `p` (assumed
/// prime). Returns 1 if a square root exists and was set; 0 otherwise.
///
/// Mirrors the C reference: dispatches by `p mod 8` to specialised
/// closed forms for `p == 3 (mod 4)` and `p == 5 (mod 8)`, falling back
/// to Tonelli-Shanks for the remaining `p == 1 (mod 8)` case.
pub fn ibz_sqrt_mod_p(sqrt: &mut Ibz, a: &Ibz, p: &Ibz) -> i32 {
    // Case a == 0 mod p: result is 0. The reference sets sqrt = 0 then
    // continues into the dispatch; we mirror that early-set behaviour.
    let p_abs = p.0.abs();
    let a_mod_zero_check = num_integer::Integer::mod_floor(&a.0, &p_abs);
    if a_mod_zero_check.is_zero() {
        sqrt.0 = BigInt::zero();
    }
    // amod = a mod p, lifted to [0, p)
    let mut amod = num_integer::Integer::mod_floor(&a.0, &p_abs);
    if amod.sign() == Sign::Minus {
        amod += &p_abs;
    }
    // Legendre symbol must be 1 for a square root to exist.
    if ibz_legendre(&Ibz(amod.clone()), p) != 1 {
        return 0;
    }
    let four = BigInt::from(4u32);
    let eight = BigInt::from(8u32);
    let p_mod_4 = num_integer::Integer::mod_floor(&p_abs, &four);
    let p_mod_8 = num_integer::Integer::mod_floor(&p_abs, &eight);
    if p_mod_4 == BigInt::from(3u32) {
        // p % 4 == 3: sqrt = a^((p+1)/4) mod p.
        let exp = (&p_abs + 1u32) >> 2;
        sqrt.0 = amod.modpow(&exp, &p_abs);
        return 1;
    }
    if p_mod_8 == BigInt::from(5u32) {
        let exp = (&p_abs - 1u32) >> 2; // (p-1)/4
        let t = amod.modpow(&exp, &p_abs);
        if t.is_one() {
            let exp2 = (&p_abs + 3u32) >> 3; // (p+3)/8
            sqrt.0 = amod.modpow(&exp2, &p_abs);
        } else {
            let exp2 = (&p_abs - 5u32) >> 3; // (p-5)/8
            let a4: BigInt = &amod << 2; // 4*a
            let t2 = a4.modpow(&exp2, &p_abs);
            let a2: BigInt = &amod << 1; // 2*a
            let prod = &a2 * &t2;
            sqrt.0 = num_integer::Integer::mod_floor(&prod, &p_abs);
        }
        return 1;
    }
    // p % 8 == 1: Tonelli-Shanks.
    let pm1 = &p_abs - 1u32;
    let mut e = 0u32;
    let mut q = pm1.clone();
    while !num_integer::Integer::is_odd(&q) {
        q >>= 1;
        e += 1;
    }
    // Find a non-residue.
    let mut qnr = BigInt::from(2u32);
    while ibz_legendre(&Ibz(qnr.clone()), p) != -1 {
        qnr += 1u32;
    }
    let z_init = qnr.modpow(&q, &p_abs);
    let mut z = z_init.clone();
    let mut y = amod.modpow(&q, &p_abs);
    let _y_init = y.clone(); // mirror the C variable shape

    let exp_half = (&q + 1u32) >> 1;
    let mut x = amod.modpow(&exp_half, &p_abs);

    let mut exp = BigInt::one() << (e.saturating_sub(2));
    let two = BigInt::from(2u32);
    for _ in 0..e {
        let b = y.modpow(&exp, &p_abs);
        if b == pm1 {
            x = num_integer::Integer::mod_floor(&(&x * &z), &p_abs);
            let zz = num_integer::Integer::mod_floor(&(&z * &z), &p_abs);
            y = num_integer::Integer::mod_floor(&(&y * &zz), &p_abs);
        }
        z = z.modpow(&two, &p_abs);
        exp >>= 1;
    }
    sqrt.0 = x;
    1
}
