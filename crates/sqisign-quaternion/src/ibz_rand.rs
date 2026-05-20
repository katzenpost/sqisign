//! RNG-driven `ibz_t` helpers.
//!
//! Mirrors the `ibz_rand_*` family in
//! `the-sqisign/src/quaternion/ref/generic/intbig.c` and the
//! random-prime helper in `integers.c`. Every entry point in this module
//! takes `&mut impl RngSource` and draws bytes through it; the C
//! reference reaches a thread-local global instead, but the byte
//! sequence the algorithm consumes is identical. The differential
//! boundary is at the *value* level: for a given seeded
//! [`sqisign_common::CtrDrbg`], the C reference and this port emit the
//! same `ibz_t` value, encoded with [`Ibz::to_canonical_bytes`].
//!
//! ## The `mp_limb_t` wire format
//!
//! `ibz_rand_interval` in the reference allocates an `mp_limb_t[len_limbs]`
//! VLA, then writes exactly `len_bytes` bytes into it via `randombytes`. On
//! the 64-bit little-endian host the dump harness targets (and that the
//! Rust port targets), `mp_limb_t` is a 64-bit little-endian word. The
//! bytes `randombytes` writes therefore fill the buffer in little-endian
//! byte order, low limb first. The top limb has its high bits masked off
//! with a precomputed `mask = ((mp_limb_t)-1) >> ((sizeof_limb_bits -
//! len_bits) % sizeof_limb_bits)`. The masked limb array is then
//! interpreted as a non-negative integer via `mpz_roinit_n`, which on
//! little-endian platforms is identical to reading the byte buffer as a
//! little-endian unsigned integer of `len_limbs * 8` bytes.
//!
//! We mirror that: read `len_bytes` random bytes, zero-extend to the next
//! multiple of 8, mask the top 8 bytes per the reference's `mask`
//! formula, then call [`BigInt::from_bytes_le`] (with `Sign::Plus`) to
//! obtain the candidate.
//!
//! ## Why `ibz_get` is fine in `ibz_rand_interval_i`
//!
//! `ibz_rand_interval_i` operates strictly on `i32`-range values: the
//! reference's `ibz_set(rand, rand32)` call writes the low 31 bits plus
//! sign. We mirror that with [`crate::ibz::ibz_set`] directly, no
//! big-int detour needed.

use num_bigint::{BigInt, Sign};

use sqisign_common::RngSource;

use crate::ibz::{
    ibz_add, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero, ibz_pow, ibz_probab_prime,
    ibz_set, ibz_sub, Ibz,
};

/// Bits per `mp_limb_t` on the target the cdump harness builds against.
/// See `tools/cdump/CMakeLists.txt`: `GMP_LIMB_BITS=64` is added to the
/// compile defines, matching the 64-bit little-endian host. Changing
/// this constant in isolation would silently desync the differential
/// boundary.
const LIMB_BITS: usize = 64;
const LIMB_BYTES: usize = LIMB_BITS / 8;

/// `ibz_rand_interval(rng, rand, a, b)`: sample a uniformly random
/// integer in the closed interval `[a, b]`. Returns 1 on success, 0 on
/// failure (the reference can only fail when `randombytes` itself fails,
/// which the trait contract treats as a panic, so this path always
/// returns 1 in practice).
///
/// Mirrors the C reference verbatim: same rejection loop, same wire
/// format for the candidate (see module docs). The `rng` argument
/// replaces the reference's thread-local DRBG global; bytes drawn match
/// what the reference would have drawn if seeded identically.
pub fn ibz_rand_interval<R: RngSource>(rng: &mut R, rand: &mut Ibz, a: &Ibz, b: &Ibz) -> i32 {
    // bmina = b - a
    let mut bmina = Ibz::zero();
    ibz_sub(&mut bmina, b, a);

    // Empty interval (a == b): write a, done.
    if bmina.0.sign() == Sign::NoSign {
        *rand = a.clone();
        return 1;
    }

    // len_bits = mpz_sizeinbase(bmina, 2): the position of the leading
    // 1-bit (one more than the index of the most significant set bit).
    // `num_bigint`'s `bits()` returns this same value for non-zero
    // unsigned magnitudes.
    let len_bits = bmina.0.bits() as usize;
    let len_bytes = len_bits.div_ceil(8);
    let len_limbs = len_bytes.div_ceil(LIMB_BYTES);
    let buf_bytes = len_limbs * LIMB_BYTES;

    // Top-limb mask, matching the reference's
    // `mask = ((mp_limb_t)-1) >> (sizeof_limb_bits - len_bits) % sizeof_limb_bits`.
    // On a 64-bit limb the reduction modulo 64 collapses the
    // unsigned-wrap case (`len_bits` exceeding the limb width) into a
    // small right-shift; for the canonical case `len_bits == 64` this
    // gives 0 and the mask becomes all-ones.
    let shift = (LIMB_BITS.wrapping_sub(len_bits)) % LIMB_BITS;
    let mask: u64 = u64::MAX >> shift;

    let mut buf = vec![0u8; buf_bytes];

    loop {
        // Bytes after `len_bytes` are left at zero (the reference's VLA
        // is uninitialised garbage there, but `mask` zeroes the same
        // bits anyway, so the canonical value is identical).
        for byte in &mut buf[..len_bytes] {
            *byte = 0;
        }
        rng.fill(&mut buf[..len_bytes]);

        // Apply the top-limb mask. On the 64-bit little-endian target
        // the top limb is bytes `buf[buf_bytes-8..buf_bytes]` in LE
        // order. Read, mask, write back.
        let top_lo = buf_bytes - LIMB_BYTES;
        let mut top = u64::from_le_bytes(buf[top_lo..buf_bytes].try_into().unwrap());
        top &= mask;
        buf[top_lo..buf_bytes].copy_from_slice(&top.to_le_bytes());

        let candidate = BigInt::from_bytes_le(Sign::Plus, &buf);
        if candidate <= bmina.0 {
            *rand = Ibz::new(candidate);
            ibz_add(&mut bmina, rand, a);
            *rand = bmina;
            return 1;
        }
    }
}

/// `ibz_rand_interval_i(rng, rand, a, b)`: sample uniformly from
/// `[a, b)` (an open upper bound, mirroring the reference signature
/// where `b > a`). Both bounds must be non-negative and fit in `i32`;
/// the result fits too.
pub fn ibz_rand_interval_i<R: RngSource>(rng: &mut R, rand: &mut Ibz, a: i32, b: i32) -> i32 {
    assert!(
        a >= 0 && b >= 0 && b > a,
        "ibz_rand_interval_i: a={a} b={b}"
    );

    let diff = (b - a) as u32;

    // mask = (1 << (32 - clz(diff))) - 1, the smallest power-of-two
    // mask whose width covers `diff`. `u32::leading_zeros` mirrors the
    // reference's `__builtin_clz` exactly.
    let mask: u32 = (1u32 << (32 - diff.leading_zeros())).wrapping_sub(1);

    debug_assert!(mask >= diff && mask < 2 * diff);

    let mut rand32: u32;
    loop {
        // The reference reads `sizeof(rand32) = 4` bytes into a signed
        // `int32_t` via `randombytes((unsigned char *)&rand32, sizeof(rand32))`
        // then bitwise-ANDs with the mask. On a little-endian host the
        // four bytes appear in low-to-high byte order, so a `from_le_bytes`
        // reinterpretation matches.
        let mut bytes = [0u8; 4];
        rng.fill(&mut bytes);
        rand32 = u32::from_le_bytes(bytes) & mask;
        if rand32 <= diff {
            break;
        }
    }

    let r = (rand32 as i32).wrapping_add(a);
    ibz_set(rand, r);
    1
}

/// `ibz_rand_interval_minm_m(rng, rand, m)`: sample uniformly from
/// `[-m, m]`, where `m >= 0`. Implemented as the reference does:
/// sample in `[0, 2m]` and shift by `m`.
pub fn ibz_rand_interval_minm_m<R: RngSource>(rng: &mut R, rand: &mut Ibz, m: i32) -> i32 {
    assert!(m >= 0, "ibz_rand_interval_minm_m: m={m}");
    let mut m_big = Ibz::zero();
    ibz_set(&mut m_big, m);
    let cur = m_big.clone();
    ibz_add(&mut m_big, &cur, &cur);

    let ret = ibz_rand_interval(rng, rand, &ibz_const_zero(), &m_big);

    // *rand -= m
    let mut shift = Ibz::zero();
    ibz_set(&mut shift, m);
    let saved = rand.clone();
    ibz_sub(rand, &saved, &shift);

    ret
}

/// `ibz_rand_interval_bits(rng, rand, m)`: sample uniformly from
/// `[-(2^m), 2^m]` then shift by `-m`, matching the reference. The
/// post-shift result lies in `[-(2^m + m), 2^m - m]`; this is a
/// faithful port of the C oddity (the reference subtracts `m` from
/// the sampled value, not `2^m`, so the output is biased toward
/// the negative side by `m`). See the C code in `intbig.c`.
pub fn ibz_rand_interval_bits<R: RngSource>(rng: &mut R, rand: &mut Ibz, m: u32) -> i32 {
    // high = 2^m
    let mut high = Ibz::zero();
    ibz_pow(&mut high, &ibz_const_two(), m);

    // low = -high
    let mut low = Ibz::zero();
    crate::ibz::ibz_neg(&mut low, &high);

    let ret = ibz_rand_interval(rng, rand, &low, &high);
    if ret != 1 {
        return ret;
    }

    // *rand -= m (treating m as an unsigned long, matching the C cast).
    let mut shift = Ibz::zero();
    // `m` is `u32`; the reference passes it via `mpz_sub_ui(rand, rand,
    // (unsigned long)m)`. `i32` covers `u32 <= i32::MAX`, but
    // `ibz_set(i32)` would truncate values above `2^31 - 1`. The C
    // reference is called with small `m` only (bitsizes), so this is
    // safe; assert nonetheless.
    assert!(m <= i32::MAX as u32, "ibz_rand_interval_bits: m too large");
    ibz_set(&mut shift, m as i32);
    let saved = rand.clone();
    ibz_sub(rand, &saved, &shift);

    1
}

/// `ibz_generate_random_prime(p, is3mod4, bitsize, probability_test_iterations)`.
///
/// Repeatedly sample a candidate uniformly from `[2^(bitsize-1-is3mod4),
/// 2^(bitsize-is3mod4))` and adjust to the correct residue (`p % 2 == 1`
/// always; `p % 4 == 3` if `is3mod4`), then test primality via the
/// reference's Miller-Rabin (`ibz_probab_prime`). Returns the same value
/// `ibz_probab_prime` reports on the accepted candidate (1 if certified
/// prime; the reference's loop only exits on a `found = 1`, so the
/// return is always 1 in practice).
pub fn ibz_generate_random_prime<R: RngSource>(
    rng: &mut R,
    p: &mut Ibz,
    is3mod4: i32,
    bitsize: i32,
    probability_test_iterations: i32,
) -> i32 {
    assert!(bitsize != 0, "ibz_generate_random_prime: bitsize == 0");
    let mut found = 0i32;
    let mut two_pow = Ibz::zero();
    let mut two_powp = Ibz::zero();
    let is3mod4_bit = (is3mod4 != 0) as u32;
    ibz_pow(
        &mut two_pow,
        &ibz_const_two(),
        ((bitsize - 1) as u32).wrapping_sub(is3mod4_bit),
    );
    ibz_pow(
        &mut two_powp,
        &ibz_const_two(),
        (bitsize as u32).wrapping_sub(is3mod4_bit),
    );

    while found == 0 {
        ibz_rand_interval(rng, p, &two_pow, &two_powp);
        let saved = p.clone();
        ibz_add(p, &saved, &saved);
        if is3mod4 != 0 {
            let saved = p.clone();
            ibz_add(p, &saved, &saved);
            let saved = p.clone();
            ibz_add(p, &ibz_const_two(), &saved);
        }
        let saved = p.clone();
        ibz_add(p, &ibz_const_one(), &saved);

        found = ibz_probab_prime(p, probability_test_iterations);
    }

    // Touch ibz_cmp to keep the symbol table mirror visible (mirrors C's
    // unused-but-exported helper paths).
    let _ = ibz_cmp;
    found
}
