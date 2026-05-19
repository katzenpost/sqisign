//! SQIsign `gf`: GF(p) and GF(p^2) arithmetic, the performance-critical
//! layer.
//!
//! Mirrors `vendor/the-sqisign/src/gf`. **Phase 1, unit 3.** This is a
//! genuine reimplementation, not a standardized primitive wired in: the
//! level-1 generic field is a self-contained `monty.py`-generated
//! word-array core (`vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c`,
//! `spint = uint64_t`, `Nlimbs 5`, `Radix 51`, prime
//! `p5248 = 5 * 2^248 - 1`), with no external-crate semantic risk and no
//! numeric dependency, exactly as `sqisign-mp` is for the `mp` word-array
//! core.
//!
//! ## Representation
//!
//! A field element is an `fp_t`: five `u64` limbs in an *unsaturated
//! radix-2^51* positional form (so `Nlimbs * Radix = 255 >= 251 = Nbits`).
//! The representation is deliberately **redundant**: a given residue class
//! has many limb encodings, and the reference's arithmetic neither expects
//! nor produces a unique canonical form (only `redc`/`modfsb`, ported
//! later, fully canonicalize). Equality of field elements is therefore
//! *not* raw-limb equality; it is the reference's own `modcmp` (also
//! ported later). The differential boundary here is the raw five-limb
//! representation: the port must be **bit-identical** to the reference at
//! the `fp_t` boundary, not merely congruent modulo `p`.
//!
//! Ported so far:
//! - [`fp_add`] is `fp_add(out, a, b)`, the thin wrapper the reference
//!   defines over `modadd`: modular addition reducing to *less than 2p*
//!   (not fully canonical). Faithfully reproduces the reference's
//!   `prop` carry-propagation helper, including its **signed** 51-bit
//!   arithmetic right shift of the running carry and its choice to leave
//!   limb 4 unmasked. No upstream defect observed; every committed
//!   C-derived vector replays bit-for-bit.
//! - [`fp_sub`] is `fp_sub(out, a, b)`, the thin wrapper the reference
//!   defines over `modsub`: modular subtraction reducing to *less than
//!   2p* (not fully canonical). It reuses the same `prop` helper. Unlike
//!   `modadd`, `modsub` does *not* pre-add `2p` before the first `prop`:
//!   it subtracts limbwise, runs `prop` (whose all-ones `carry` mask
//!   signals the difference went negative), and only then conditionally
//!   adds `2p` back in the redundant form. No upstream defect observed;
//!   every committed C-derived vector replays bit-for-bit.
//! - [`fp_neg`] is `fp_neg(out, a)`, the thin wrapper the reference
//!   defines over `modneg`: modular negation reducing to *less than 2p*
//!   (not fully canonical). It is the **unary analogue of `modsub`**:
//!   where `modsub` forms `a[i] - b[i]`, `modneg` forms `0 - b[i]`
//!   limbwise (the implicit minuend is zero), then runs the *identical*
//!   `prop` plus conditional `2p` correction tail `modsub` uses. It
//!   reuses the same `prop` helper and `TWO_P4` constant. No upstream
//!   defect observed; every committed C-derived vector replays
//!   bit-for-bit, and the all-zero input maps to the bit-exact canonical
//!   all-zero representative.
//! - [`fp_copy`] is `fp_copy(out, a)`, the thin wrapper the reference
//!   defines over `modcpy`: a plain five-limb assignment, `out[i] = a[i]`
//!   for `i` in `0..5`. No `prop`, no `2p` correction, no reduction. It
//!   is bit-exact on every input by construction, including non-canonical
//!   limbs: the recorded output equals the recorded input limb for limb
//!   on the full 1012-vector battery. No upstream defect observed.
//! - [`fp_set_zero`] is `fp_set_zero(x)`, the thin wrapper the reference
//!   defines over `modzer`: a plain five-limb zero-fill, `x[i] = 0` for
//!   `i` in `0..5`. No `prop`, no `2p` correction, no reduction; the
//!   recorded output is the bit-exact canonical all-zero representative
//!   regardless of the destination's prior contents. The differential
//!   boundary varies the destination *pre-fill* as the "input" (a setter
//!   has no field argument); a no-op or partial-write port is caught at
//!   that boundary because each non-zero pre-fill leaves a visible
//!   residue if any limb is forgotten. No upstream defect observed; every
//!   committed C-derived vector replays bit-for-bit.
//! - [`fp_set_one`] is `fp_set_one(x)`, the thin wrapper the reference
//!   defines over `modone`, which writes positional `1` and then calls
//!   `nres(a, a)` to convert it to its Montgomery representative. The
//!   on-the-wire output of `fp_set_one` is therefore the Montgomery
//!   `ONE`, `[0x19, 0, 0, 0, 0x300000000000]`, the same bit pattern
//!   the reference exposes as `extern const ONE` at lines 526..530 of
//!   `fp_p5248_64.c`. `nres`/`modmul`/`R2` are not yet ported, but the
//!   output of `fp_set_one` is a single fixed constant; the port
//!   writes it directly and bit-matches every recorded reference
//!   output. Like [`fp_set_zero`], the differential boundary varies
//!   the destination *pre-fill* as the "input"; any limb the port
//!   forgets to write leaves a visible non-zero residue that diverges
//!   from the recorded output. No upstream defect observed; every
//!   committed C-derived vector replays bit-for-bit.
//!
//! Correctness is established as for the whole port: every committed
//! C-derived vector is replayed and bit-compared (`tests/`). Equivalence
//! to the reference is proven, not presumed.
#![forbid(unsafe_code)]

/// Number of `u64` limbs in a level-1 `fp_t`.
///
/// `NWORDS_FIELD == 5` for the `RADIX_64`, non-Broadwell generic build
/// (`fp_constants.h`): `Nlimbs == 5`, `Radix == 51`.
pub const NWORDS_FIELD: usize = 5;

/// A level-1 GF(p) field element: five `u64` limbs in unsaturated
/// radix-2^51 positional form. See the crate-level note on the redundant,
/// non-canonical representation.
pub type Fp = [u64; NWORDS_FIELD];

/// The radix: each canonical limb carries 51 bits.
const RADIX: u32 = 51;

/// `(1 << 51) - 1`, the per-limb mask `prop` applies to limbs 0..=3.
const MASK51: u64 = (1u64 << RADIX) - 1;

/// `2 * p4`, where `p4 == 0x500000000000` is the prime's contribution at
/// limb 4. `modadd` adds `2p` in redundant form (`+2` at limb 0, this
/// subtracted at limb 4) and conditionally adds it back.
const TWO_P4: u64 = 0xa00000000000;

/// Propagate carries across a five-limb value in place; return an
/// all-ones mask iff the value is "negative" in the redundant form.
///
/// Mirrors the reference's
/// `inline static spint prop(spint *n)` exactly:
///
/// ```c
/// spint mask = ((spint)1 << 51u) - (spint)1;
/// sspint carry = (sspint)n[0];
/// carry >>= 51u;            // signed (arithmetic) shift
/// n[0] &= mask;
/// for (i = 1; i < 4; i++) {
///   carry += (sspint)n[i];
///   n[i] = (spint)carry & mask;
///   carry >>= 51u;          // signed (arithmetic) shift
/// }
/// n[4] += (spint)carry;
/// return -((n[4] >> 1) >> 62u);
/// ```
///
/// Two faithfully reproduced subtleties:
/// 1. `carry` is `sspint == int64_t`; `carry >>= 51` is therefore an
///    **arithmetic** right shift (sign-extending). The port uses `i64`
///    so Rust's `>>` on the signed type matches the reference's signed
///    shift bit-for-bit.
/// 2. Limb 4 is **not** masked. The reference leaves it holding the full
///    accumulated top word (possibly far above `2^51`); the port does the
///    same. Only `redc`/`modfsb` (ported later) canonicalize it.
///
/// The returned value is `-((n[4] >> 1) >> 62)` computed in unsigned
/// `u64` (the reference's `spint` shifts here are unsigned): it is the
/// all-ones mask `0xffff_ffff_ffff_ffff` when `n[4]`'s sign bit is set,
/// else `0`.
fn prop(n: &mut Fp) -> u64 {
    let mut carry = n[0] as i64;
    carry >>= RADIX;
    n[0] &= MASK51;
    // Explicit index `1..4` mirrors the reference's `for (i = 1; i < 4;
    // i++)` exactly: each step reads then overwrites n[i] through the
    // running signed carry. An iterator rewrite would obscure the
    // bit-for-bit correspondence with the oracle, so the lint is silenced
    // locally rather than the loop reshaped.
    #[allow(clippy::needless_range_loop)]
    for i in 1..4 {
        carry = carry.wrapping_add(n[i] as i64);
        n[i] = (carry as u64) & MASK51;
        carry >>= RADIX;
    }
    n[4] = n[4].wrapping_add(carry as u64);
    (0u64).wrapping_sub((n[4] >> 1) >> 62)
}

/// Modular addition, reducing to less than `2p` (not fully canonical).
///
/// Mirrors the reference's
/// `inline static void modadd(const spint *a, const spint *b, spint *n)`:
///
/// ```c
/// n[i] = a[i] + b[i];           // i = 0..4
/// n[0] += 2u;
/// n[4] -= 0xa00000000000u;      // == 2 * p4
/// carry = prop(n);
/// n[0] -= 2u & carry;
/// n[4] += 0xa00000000000u & carry;
/// (void)prop(n);
/// ```
///
/// The `+2` at limb 0 paired with `-2*p4` at limb 4 adds `2p` in the
/// redundant representation; `prop` then signals (via its all-ones
/// `carry` mask) whether the result went negative, in which case the
/// `2p` is masked back in and carries are propagated once more. The
/// reference accepts arbitrary, non-canonical limb inputs and the port
/// reproduces its output on those too (the C-derived vectors include
/// full-width 64-bit limb inputs and pin this).
///
/// Wrapping `u64` arithmetic throughout matches the reference's
/// `spint == uint64_t` two's-complement wraparound.
fn modadd(a: &Fp, b: &Fp, n: &mut Fp) {
    for i in 0..NWORDS_FIELD {
        n[i] = a[i].wrapping_add(b[i]);
    }
    n[0] = n[0].wrapping_add(2);
    n[4] = n[4].wrapping_sub(TWO_P4);
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[4] = n[4].wrapping_add(TWO_P4 & carry);
    let _ = prop(n);
}

/// GF(p) addition `out = a + b mod p`, in the redundant radix-2^51
/// representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_add(fp_t *out, const fp_t *a,
/// const fp_t *b)`, which is the thin wrapper `modadd(*a, *b, *out)`.
/// The output is *not* fully canonical: it is a valid representative of
/// `a + b` whose limbs 0..=3 are below `2^51` (the `prop` mask) but whose
/// limb 4 is left unmasked, exactly as the reference leaves it. Compare
/// field elements with the reference's equality (`modcmp`, ported later),
/// never by raw-limb equality.
pub fn fp_add(out: &mut Fp, a: &Fp, b: &Fp) {
    modadd(a, b, out);
}

/// Modular subtraction, reducing to less than `2p` (not fully canonical).
///
/// Mirrors the reference's
/// `inline static void modsub(const spint *a, const spint *b, spint *n)`:
///
/// ```c
/// n[i] = a[i] - b[i];           // i = 0..4
/// carry = prop(n);
/// n[0] -= 2u & carry;
/// n[4] += 0xa00000000000u & carry;
/// (void)prop(n);
/// ```
///
/// The structural contrast with `modadd` is faithfully reproduced: where
/// `modadd` *pre-adds* `2p` (`n[0] += 2`, `n[4] -= 2*p4`) before the first
/// `prop`, `modsub` does no such pre-add. It subtracts limbwise, then runs
/// `prop`; `prop`'s all-ones `carry` mask signals whether the difference
/// went negative, in which case `2p` is masked in (`n[0] -= 2`,
/// `n[4] += 2*p4`) and carries are propagated once more. The conditional
/// `2p` correction tail is identical to `modadd`'s; only the missing
/// pre-add distinguishes the two. As with `modadd`, the reference accepts
/// arbitrary non-canonical limb inputs and the port reproduces its output
/// on those too (the C-derived vectors include full-width 64-bit limb
/// inputs and pin this).
///
/// Wrapping `u64` arithmetic throughout matches the reference's
/// `spint == uint64_t` two's-complement wraparound (each `a[i] - b[i]`
/// borrows by wrapping, exactly as the unsigned C subtraction does).
fn modsub(a: &Fp, b: &Fp, n: &mut Fp) {
    for i in 0..NWORDS_FIELD {
        n[i] = a[i].wrapping_sub(b[i]);
    }
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[4] = n[4].wrapping_add(TWO_P4 & carry);
    let _ = prop(n);
}

/// GF(p) subtraction `out = a - b mod p`, in the redundant radix-2^51
/// representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_sub(fp_t *out, const fp_t *a,
/// const fp_t *b)`, which is the thin wrapper `modsub(*a, *b, *out)`.
/// As with [`fp_add`], the output is *not* fully canonical: limbs 0..=3
/// are below `2^51` (the `prop` mask) but limb 4 is left unmasked,
/// exactly as the reference leaves it. Compare field elements with the
/// reference's equality (`modcmp`, ported later), never by raw-limb
/// equality.
pub fn fp_sub(out: &mut Fp, a: &Fp, b: &Fp) {
    modsub(a, b, out);
}

/// Modular negation, reducing to less than `2p` (not fully canonical).
///
/// Mirrors the reference's
/// `inline static void modneg(const spint *b, spint *n)`:
///
/// ```c
/// n[i] = (spint)0 - b[i];       // i = 0..4
/// carry = prop(n);
/// n[0] -= 2u & carry;
/// n[4] += 0xa00000000000u & carry;
/// (void)prop(n);
/// ```
///
/// This is the **unary analogue of `modsub`**: the minuend is the
/// implicit constant `0`, so each limb is `0 - b[i]` rather than
/// `a[i] - b[i]`. The `prop` plus conditional `2p` correction tail is
/// byte-for-byte the one `modsub` uses; only the implicit-zero minuend
/// distinguishes the two, exactly as the missing pre-add distinguishes
/// `modsub` from `modadd`. `prop`'s all-ones `carry` mask signals
/// whether the negation went negative in the redundant form, in which
/// case `2p` is masked back in (`n[0] -= 2`, `n[4] += 2*p4`) and carries
/// are propagated once more.
///
/// `(spint)0 - b[i]` is an *unsigned* `uint64_t` subtraction in the
/// reference (`spint == uint64_t`): it is the two's-complement negation
/// `0 - b[i]` with wraparound, reproduced here as
/// `0u64.wrapping_sub(b[i])`. As with `modadd`/`modsub`, the reference
/// accepts arbitrary non-canonical limb inputs and the port reproduces
/// its output on those too (the C-derived vectors include full-width
/// 64-bit limb inputs and pin this). The all-zero input is a fixed
/// point: every limb stays `0`, `prop` returns a zero carry mask, no
/// correction fires, and the output is the bit-exact canonical all-zero
/// representative.
fn modneg(b: &Fp, n: &mut Fp) {
    for i in 0..NWORDS_FIELD {
        n[i] = 0u64.wrapping_sub(b[i]);
    }
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(2 & carry);
    n[4] = n[4].wrapping_add(TWO_P4 & carry);
    let _ = prop(n);
}

/// GF(p) negation `out = -a mod p`, in the redundant radix-2^51
/// representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_neg(fp_t *out, const fp_t *a)`,
/// which is the thin wrapper `modneg(*a, *out)`. As with [`fp_add`] and
/// [`fp_sub`], the output is *not* fully canonical: limbs 0..=3 are
/// below `2^51` (the `prop` mask) but limb 4 is left unmasked, exactly
/// as the reference leaves it. Compare field elements with the
/// reference's equality (`modcmp`, ported later), never by raw-limb
/// equality. The lone exception is the canonical zero: `fp_neg` of the
/// all-zero representative is the bit-exact all-zero representative.
pub fn fp_neg(out: &mut Fp, a: &Fp) {
    modneg(a, out);
}

/// GF(p) assignment `out = a`, in the redundant radix-2^51 representation.
///
/// Mirrors the reference's `void fp_copy(fp_t *out, const fp_t *a)`, which
/// is the thin wrapper `modcpy(*a, *out)`:
///
/// ```c
/// inline static void modcpy(const spint *a, spint *c) {
///   for (int i = 0; i < 5; i++) c[i] = a[i];
/// }
/// ```
///
/// A plain five-limb assignment: no `prop`, no `2p` correction, no
/// reduction. The output is bit-exact equal to the input by construction,
/// including for non-canonical limbs (the reference makes no assumption
/// about its argument's range and neither does the port). Unlike
/// [`fp_add`], [`fp_sub`] and [`fp_neg`], `fp_copy` is therefore exactly
/// raw-limb-equal on its output rather than merely a valid representative
/// of the same residue class.
pub fn fp_copy(out: &mut Fp, a: &Fp) {
    *out = *a;
}

/// GF(p) zero setter `*x = 0`, in the redundant radix-2^51 representation.
///
/// Mirrors the reference's `void fp_set_zero(fp_t *x)`, which is the thin
/// wrapper `modzer(*x)`:
///
/// ```c
/// static void modzer(spint *a) {
///   for (int i = 0; i < 5; i++) a[i] = 0;
/// }
/// ```
///
/// A plain five-limb zero-fill: no `prop`, no `2p` correction, no
/// reduction. The output is the bit-exact canonical all-zero
/// representative regardless of the destination's prior contents, and so
/// is unambiguously the canonical zero in *both* the redundant-mod-`p`
/// sense and the raw-limb sense (the lone exception, like the all-zero
/// fixed point of [`fp_neg`]: every other field result is congruent only
/// up to the redundant form).
pub fn fp_set_zero(out: &mut Fp) {
    *out = [0u64; NWORDS_FIELD];
}

/// GF(p) one setter: writes the Montgomery representative of `1`,
/// `[0x19, 0, 0, 0, 0x300000000000]`, regardless of the destination's
/// prior contents.
///
/// Mirrors the reference's `void fp_set_one(fp_t *x)`, which is the thin
/// wrapper `modone(*x)`. The reference's `modone` writes positional
/// `1 = 1 * 2^0` then calls `nres(a, a)` to convert it to its
/// Montgomery n-residue form:
///
/// ```c
/// static void modone(spint *a) {
///   int i;
///   a[0] = 1;
///   for (i = 1; i < 5; i++) a[i] = 0;
///   nres(a, a);
/// }
/// ```
///
/// `nres` (and the `modmul`/`R2` precomputation it relies on) is not
/// yet ported, but the differential-boundary output of `fp_set_one` is
/// a single fixed bit pattern: the Montgomery `ONE`, identical to the
/// public `extern const ONE` defined at lines 526..530 of
/// `fp_p5248_64.c`. The port writes that constant directly. At the
/// `fp_t` boundary the two are bit-equal; when `nres` is ported it can
/// (and should) be exercised at its own boundary, leaving `fp_set_one`
/// as the same constant-write the reference effectively evaluates to.
pub fn fp_set_one(out: &mut Fp) {
    *out = MONTGOMERY_ONE;
}

/// Montgomery representative of `1` on the level-1 generic field:
/// `1 * R mod p` in the unsaturated radix-2^51 limb layout, matching
/// the reference's public `extern const ONE`
/// (`fp_p5248_64.c:526..530`).
const MONTGOMERY_ONE: Fp = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
];
