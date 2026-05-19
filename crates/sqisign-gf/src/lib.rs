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
//! - [`fp_select`] is `fp_select(d, a0, a1, ctl)`, the branchless
//!   constant-time conditional select defined in
//!   `vendor/the-sqisign/src/gf/ref/lvlx/fp.c`. The contract is narrow:
//!   `ctl` is required to be either `0x00000000` or `0xFFFFFFFF`; on
//!   `0x00000000` the destination is set to `a0`, on `0xFFFFFFFF` to
//!   `a1`. The reference computes a per-limb bit blend
//!   `d[i] = a0[i] ^ (cw & (a0[i] ^ a1[i]))`, where `cw` is `ctl`
//!   *sign-extended* to `digit_t` (a `uint64_t`) via the C cast chain
//!   `digit_t cw = (int32_t)ctl;`: the implicit widening from `int32_t`
//!   to `uint64_t` is sign-extending in C, so `cw == 0` for `ctl == 0`
//!   and `cw == 0xFFFFFFFFFFFFFFFF` for `ctl == 0xFFFFFFFF`. The port
//!   reproduces that two-step cast explicitly as `(ctl as i32) as u64`
//!   (sign-extending `i32 -> i64` then bit-casting to `u64`) to keep
//!   the bit-for-bit oracle correspondence visible at the call site.
//!   Per the reference contract any other `ctl` value is undefined; the
//!   C-derived vector battery therefore exercises only the two declared
//!   endpoints. No upstream defect observed; every committed C-derived
//!   vector replays bit-for-bit.
//! - [`fp_mul`] is `fp_mul(out, a, b)`, the thin wrapper the reference
//!   defines over `modmul`: Montgomery modular multiplication on the
//!   level-1 generic field. `modmul` is a Granger-Scott style schoolbook
//!   multiplier with the Montgomery reduction folded **inline**, taking
//!   advantage of the special structure of `p5248` (only limb 4 of the
//!   prime, `p4 == 0x500000000000`, is non-zero). The diagonals
//!   `a[i] * b[j]` for each column `i + j == k` are summed into a 128-bit
//!   accumulator `t`; at columns 0..=4 the low 51 bits of `t` are stored
//!   as `v0..v4` and `t` is shifted right by 51; starting at column 4 each
//!   column also folds in `v_{k-4} * p4`, performing the Montgomery
//!   reduction inline with the multiplication. The output is the redundant
//!   reduced representative: limbs 0..=3 are below `2^51` (the column
//!   mask), limb 4 carries the full residual `t` *unmasked* (no `& mask`
//!   on the final write). The port transcribes the reference's column
//!   structure identifier-for-identifier and statement-for-statement;
//!   `(dpint)a[i] * b[j]` becomes `(a[i] as u128) * (b[j] as u128)`, the
//!   `(spint)t & mask` masks become `(t as u64) & MASK51`, and the
//!   `t >>= 51` shifts apply directly to `u128`. No upstream defect
//!   observed; every committed C-derived vector replays bit-for-bit. The
//!   commutativity hypothesis `fp_mul(a,b) == fp_mul(b,a)` was checked
//!   bit-exactly across the full 1144-vector battery before being pinned
//!   in the property suite (see `fp_mul_props.rs` for justification: the
//!   column sums are operand-symmetric over the unordered diagonal
//!   `a[i] * b[j]` set, and `u128` accumulation is associative and
//!   commutative).
//! - [`fp_cswap`] is `fp_cswap(a, b, ctl)`, the branchless constant-time
//!   conditional swap defined in `fp_p5248_64.c`. Its contract is
//!   wider than [`fp_select`]'s: only the LSB of `ctl` is consulted
//!   (the reference's wrapper narrows it via `(int)(ctl & 0x1)`), so
//!   any `ctl` with `ctl & 1 == 0` is a no-op and any with `ctl & 1
//!   == 1` swaps `a` and `b` limb for limb. The underlying `modcsw`
//!   computes the swap as a cross-multiplication of the two operands
//!   by `c0 = (1 - b) + r` and `c1 = b + r` (where `r` is the rotating
//!   constant `0x3cc3c33c5aa5a55a` and the arithmetic is unsigned
//!   wrapping on `spint == uint64_t`), subtracting `w = r * (t + s)`
//!   from each cross product so the `r` term cancels and the residue
//!   is `s` or `t` depending on `b`. The port reproduces the
//!   `wrapping_add`/`wrapping_sub`/`wrapping_mul` chain on `u64`
//!   bit-for-bit, and the `(ctl & 1) as u64` narrowing exactly mirrors
//!   the reference's `(int)(ctl & 0x1)`. No upstream defect observed;
//!   every committed C-derived vector replays bit-for-bit, both
//!   endpoints are pinned, and the LSB-only contract is pinned by
//!   recording two non-canonical ctl values (`0xfffffffe` acts as `0`,
//!   `0xffffffff` acts as `1`).
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

/// Branchless constant-time conditional select on `fp_t`.
///
/// Mirrors the reference's
/// `void fp_select(fp_t *d, const fp_t *a0, const fp_t *a1, uint32_t ctl)`
/// from `vendor/the-sqisign/src/gf/ref/lvlx/fp.c` exactly:
///
/// ```c
/// /*
///  * If ctl == 0x00000000, then *d is set to a0
///  * If ctl == 0xFFFFFFFF, then *d is set to a1
///  * ctl MUST be either 0x00000000 or 0xFFFFFFFF.
///  */
/// void fp_select(fp_t *d, const fp_t *a0, const fp_t *a1, uint32_t ctl) {
///     digit_t cw = (int32_t)ctl;
///     for (unsigned int i = 0; i < NWORDS_FIELD; i++) {
///         (*d)[i] = (*a0)[i] ^ (cw & ((*a0)[i] ^ (*a1)[i]));
///     }
/// }
/// ```
///
/// Per the reference's documented contract, `ctl` must be either
/// `0x00000000` or `0xFFFFFFFF`; any other value is undefined behaviour
/// at the reference, and the C-derived vector battery therefore only
/// exercises the two declared endpoints. The port does **not** narrow
/// the type to `bool` because the differential boundary records `ctl`
/// as the raw `u32` the reference takes.
///
/// The cast chain `digit_t cw = (int32_t)ctl;` is the load-bearing
/// subtlety. In C, the right-hand side casts `uint32_t` to `int32_t`
/// (a bit-preserving reinterpretation under two's complement) then the
/// implicit assignment to `digit_t == uint64_t` is a *widening* from a
/// signed type, which is sign-extending. So `ctl == 0` yields `cw == 0`
/// and `ctl == 0xFFFFFFFF` (interpreted as the signed `-1`) yields
/// `cw == 0xFFFFFFFFFFFFFFFF`. The port reproduces this as the
/// explicit two-step `(ctl as i32) as u64`: `as i32` is the
/// bit-preserving `uint32_t -> int32_t` reinterpret, and `as u64` on
/// a signed type sign-extends (`i32 -> i64`) before bit-casting to
/// `u64`. The resulting `cw` is bit-for-bit equal to the reference's,
/// preserving the oracle correspondence at the limb-XOR step.
pub fn fp_select(d: &mut Fp, a0: &Fp, a1: &Fp, ctl: u32) {
    let cw = (ctl as i32) as u64;
    for i in 0..NWORDS_FIELD {
        d[i] = a0[i] ^ (cw & (a0[i] ^ a1[i]));
    }
}

/// Rotating constant `r` used by `modcsw`'s cross-multiplication. The
/// reference picks this 64-bit pattern (alternating nibbles
/// `3c c3 c3 3c 5a a5 a5 5a`) so that `c0 = (1 - b) + r` and
/// `c1 = b + r` have the same algebraic role in the per-limb update and
/// the `w = r * (t + s)` subtraction cancels the `r` contribution on
/// both `f` and `g`, leaving only the `s`/`t` cross-residue.
const MODCSW_R: u64 = 0x3cc3_c33c_5aa5_a55a;

/// Branchless constant-time conditional swap on the LSB of `b`.
///
/// Mirrors the reference's
/// `static void modcsw(int b, volatile spint *g, volatile spint *f)`
/// from `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:409..424`
/// exactly:
///
/// ```c
/// static void modcsw(int b, volatile spint *g, volatile spint *f) {
///   int i;
///   spint c0, c1, s, t, w;
///   spint r = 0x3cc3c33c5aa5a55au;
///   c0 = (1 - b) + r;
///   c1 = b + r;
///   for (i = 0; i < 5; i++) {
///     s = g[i];
///     t = f[i];
///     w = r * (t + s);
///     f[i] = c0 * t + c1 * s;
///     f[i] -= w;
///     g[i] = c0 * s + c1 * t;
///     g[i] -= w;
///   }
/// }
/// ```
///
/// Two faithfully reproduced subtleties:
///
/// 1. All arithmetic is on `spint == uint64_t` with two's-complement
///    wraparound. The port uses `wrapping_add`/`wrapping_sub`/
///    `wrapping_mul` on `u64` throughout so the bit pattern in every
///    intermediate matches the reference's modulo-`2^64` value. Both
///    `(1 - b) + r` (where `1 - b` is the unsigned subtraction
///    `1u64.wrapping_sub(b)`) and `c0 * t + c1 * s - r * (t + s)` rely
///    on this exact wraparound to cancel the `r` contribution.
/// 2. `b` is `int` in the reference and the wrapper narrows
///    `ctl` to `(int)(ctl & 0x1)`. Only the LSB matters; the higher
///    bits of `ctl` are dropped before `b` is consumed. The port
///    therefore does `let b = (ctl & 1) as u64` and feeds that to the
///    `c0`/`c1` formulae directly: the resulting `u64` `b` is `0` or
///    `1`, and `1u64.wrapping_sub(b)` is `1` or `0` respectively. The
///    algebraic identity verified at both endpoints:
///    - `b == 0`: `c0 = (1 + r)`, `c1 = r`. Then
///      `f' = (1+r)*t + r*s - r*(t+s) = t + r*t + r*s - r*t - r*s = t`
///      and symmetrically `g' = s`; the pair is unchanged.
///    - `b == 1`: `c0 = r`, `c1 = (1 + r)`. Then
///      `f' = r*t + (1+r)*s - r*(t+s) = r*t + s + r*s - r*t - r*s = s`
///      and symmetrically `g' = t`; the pair is swapped.
///
/// The `volatile` qualifier on the reference's pointers is a side-effect
/// hint for the C compiler (preventing it from elideing the writes or
/// reordering them across the load) and has no observable effect at the
/// differential `fp_t` boundary; the port reproduces the limb-by-limb
/// reads and writes without any equivalent annotation (constant-time
/// guarantees in Rust are out of scope for the port-correctness
/// contract, exactly as for the existing constant-time ports
/// `select_ct`, `swap_ct`, and [`fp_select`]).
fn modcsw(b: u64, g: &mut Fp, f: &mut Fp) {
    let r = MODCSW_R;
    let c0 = 1u64.wrapping_sub(b).wrapping_add(r);
    let c1 = b.wrapping_add(r);
    // Explicit index `0..5` mirrors the reference's `for (i = 0; i < 5;
    // i++)` exactly: each iteration reads then overwrites g[i] and f[i]
    // through the running products. An iterator rewrite would obscure
    // the bit-for-bit correspondence with the oracle, so the lint is
    // silenced locally rather than the loop reshaped.
    #[allow(clippy::needless_range_loop)]
    for i in 0..NWORDS_FIELD {
        let s = g[i];
        let t = f[i];
        let w = r.wrapping_mul(t.wrapping_add(s));
        let new_f = c0.wrapping_mul(t).wrapping_add(c1.wrapping_mul(s));
        f[i] = new_f.wrapping_sub(w);
        let new_g = c0.wrapping_mul(s).wrapping_add(c1.wrapping_mul(t));
        g[i] = new_g.wrapping_sub(w);
    }
}

/// Branchless constant-time conditional swap on `fp_t`.
///
/// Mirrors the reference's
/// `void fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)` from
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:592..596`:
///
/// ```c
/// void
/// fp_cswap(fp_t *a, fp_t *b, uint32_t ctl)
/// {
///     modcsw((int)(ctl & 0x1), *a, *b);
/// }
/// ```
///
/// Only the **LSB** of `ctl` is consulted: the wrapper narrows it to
/// `(int)(ctl & 0x1)` before passing it to `modcsw` (see `modcsw` for
/// the cross-multiplication and the algebraic verification). When
/// `ctl & 1 == 0`, `a` and `b` are left unchanged limb for limb; when
/// `ctl & 1 == 1`, `a` and `b` are swapped limb for limb. The
/// higher bits of `ctl` are not consulted, so `ctl = 0xfffffffe` is a
/// no-op (LSB clear) and `ctl = 0xffffffff` is a swap (LSB set), in
/// contrast to [`fp_select`] which requires the full 32-bit mask.
///
/// The port performs the same `ctl & 1` narrowing as
/// `(ctl & 1) as u64`, matching the recorded boundary's encoding
/// (the C harness records the raw `uint32_t` and the port consumes
/// it identically).
pub fn fp_cswap(a: &mut Fp, b: &mut Fp, ctl: u32) {
    let bit = (ctl & 1) as u64;
    modcsw(bit, a, b);
}

/// Prime's contribution at limb 4, `p4 == 0x500000000000`. The
/// level-1 prime `p5248 = 5 * 2^248 - 1` has all-zero limbs at indices
/// 0..=3 in the radix-2^51 layout (`5 * 2^248 = 5 * 2^(51*4 + 44) =
/// (5 << 44) << (51 * 4) = 0x500000000000 << (51 * 4)`); `modmul`
/// exploits this by folding only `v_k * p4` into the running accumulator
/// at columns 4..=8 (the Montgomery reduction inline with the
/// multiplication).
const P4: u64 = 0x500000000000;

/// Montgomery modular multiplication, `c = a * b mod 2p` in the redundant
/// radix-2^51 representation.
///
/// Mirrors the reference's
/// `inline static void modmul(const spint *a, const spint *b, spint *c)`
/// from `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:98..153`
/// statement-for-statement; the column structure and the order of every
/// `t +=` are preserved so the bit-exact oracle correspondence at the
/// `fp_t` boundary is visible at the source level.
///
/// Granger-Scott structure for the special prime `p5248`:
///
/// - The product `a * b` is computed schoolbook, columns 0..=8, into a
///   128-bit accumulator `t` (the reference's `dpint == __uint128_t`,
///   ported here as `u128`). At each column `k` the low 51 bits of `t`
///   are stored (as `v0..v4` for columns 0..=4 and as the output limbs
///   `c[0..=3]` for columns 5..=8), then `t >>= 51`. The mask is
///   `MASK51 = (1 << 51) - 1`.
/// - The Montgomery reduction is folded inline by adding `v_{k-4} * p4`
///   to `t` at columns 4..=8, where `p4 == 0x500000000000` is the
///   prime's only non-zero limb in the radix-2^51 layout. Because the
///   reduction touches only column 4 of the prime, the per-column work
///   is a single extra `v_k * p4` add rather than a full
///   `q * p` Schoolbook reduction.
/// - The final write `c[4] = (spint)t` is a **full 64-bit truncation of
///   the residual `t`, NOT masked** (no `& mask`), matching the
///   reference's last line exactly: limb 4 of the output is left
///   unmasked, exactly as `modadd`/`modsub`/`modneg` leave it (only
///   `redc`/`modfsb`, ported later, fully canonicalise).
///
/// The transcription rules applied uniformly:
/// - `spint`/`dpint`/`udpint` -> `u64` / `u128` / `u128`.
/// - `(dpint)a[i] * b[j]` -> `(a[i] as u128) * (b[j] as u128)`. Each
///   product is at most `(2^64 - 1)^2 = 2^128 - 2^65 + 1 < 2^128`, so
///   the multiplication never overflows `u128` and plain `*` is the
///   bit-exact equivalent of the C `(__uint128_t)a * b`.
/// - `t += ...` -> `t = t.wrapping_add(...)`. The reference's `t +=` is
///   on `__uint128_t`, which is *unsigned wraparound* by the C standard.
///   For canonical inputs (limbs `< 2^51`) the reference's own analysis
///   shows `t` never exceeds `~2^104.4` so wraparound is invisible; the
///   differential battery, however, deliberately feeds full-width random
///   limbs (the cdump harness's stated intent: "the reference accepts
///   any limbs and the port must match it bit-for-bit on them too"), and
///   on those `t` can exceed `2^128` after summing five partial
///   products. The port therefore uses `wrapping_add` throughout to
///   reproduce the C wraparound semantics exactly, otherwise Rust's
///   debug-mode overflow checks would trap on inputs the C reference
///   handles silently. `wrapping_add` is associative and commutative on
///   `u128` modulo `2^128`, so the commutativity argument for property
///   (1) is unaffected.
/// - `(spint)t & mask` -> `(t as u64) & MASK51`.
/// - `(dpint)v_k * (dpint)p4` -> `(v_k as u128) * (P4 as u128)`. The
///   factors are at most `2^51 - 1` and `0x500000000000 < 2^47`, so the
///   product fits in 98 bits; plain `*` is safe and bit-exact.
/// - `t >>= 51` applies directly to `u128`.
///
/// Identifier mapping (each row preserves the reference's name verbatim):
/// `t`, `v0`, `v1`, `v2`, `v3`, `v4`, `c[0..=4]`, `p4` (constant inlined
/// as `P4 as u128` at each fold). The statement order and the per-column
/// number of `t += ...` lines are preserved, so the column-by-column
/// trace lines up one-to-one with `fp_p5248_64.c:100..152`.
fn modmul(a: &Fp, b: &Fp, c: &mut Fp) {
    let mut t: u128 = 0;
    t = t.wrapping_add((a[0] as u128) * (b[0] as u128));
    let v0: u64 = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[0] as u128) * (b[1] as u128));
    t = t.wrapping_add((a[1] as u128) * (b[0] as u128));
    let v1: u64 = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[0] as u128) * (b[2] as u128));
    t = t.wrapping_add((a[1] as u128) * (b[1] as u128));
    t = t.wrapping_add((a[2] as u128) * (b[0] as u128));
    let v2: u64 = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[0] as u128) * (b[3] as u128));
    t = t.wrapping_add((a[1] as u128) * (b[2] as u128));
    t = t.wrapping_add((a[2] as u128) * (b[1] as u128));
    t = t.wrapping_add((a[3] as u128) * (b[0] as u128));
    let v3: u64 = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[0] as u128) * (b[4] as u128));
    t = t.wrapping_add((a[1] as u128) * (b[3] as u128));
    t = t.wrapping_add((a[2] as u128) * (b[2] as u128));
    t = t.wrapping_add((a[3] as u128) * (b[1] as u128));
    t = t.wrapping_add((a[4] as u128) * (b[0] as u128));
    t = t.wrapping_add((v0 as u128) * (P4 as u128));
    let v4: u64 = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[1] as u128) * (b[4] as u128));
    t = t.wrapping_add((a[2] as u128) * (b[3] as u128));
    t = t.wrapping_add((a[3] as u128) * (b[2] as u128));
    t = t.wrapping_add((a[4] as u128) * (b[1] as u128));
    t = t.wrapping_add((v1 as u128) * (P4 as u128));
    c[0] = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[2] as u128) * (b[4] as u128));
    t = t.wrapping_add((a[3] as u128) * (b[3] as u128));
    t = t.wrapping_add((a[4] as u128) * (b[2] as u128));
    t = t.wrapping_add((v2 as u128) * (P4 as u128));
    c[1] = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[3] as u128) * (b[4] as u128));
    t = t.wrapping_add((a[4] as u128) * (b[3] as u128));
    t = t.wrapping_add((v3 as u128) * (P4 as u128));
    c[2] = (t as u64) & MASK51;
    t >>= 51;
    t = t.wrapping_add((a[4] as u128) * (b[4] as u128));
    t = t.wrapping_add((v4 as u128) * (P4 as u128));
    c[3] = (t as u64) & MASK51;
    t >>= 51;
    c[4] = t as u64;
}

/// GF(p) Montgomery multiplication `out = a * b * R^-1 mod p`, in the
/// redundant radix-2^51 representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_mul(fp_t *out, const fp_t *a,
/// const fp_t *b)`, which is the thin wrapper `modmul(*a, *b, *out)`.
/// As with [`fp_add`], [`fp_sub`] and [`fp_neg`], the output is *not*
/// fully canonical: limbs 0..=3 are below `2^51` (the column mask) but
/// limb 4 is left unmasked (the final `c[4] = (spint)t` is a full 64-bit
/// truncation of the residual accumulator, no `& mask`), exactly as the
/// reference leaves it. Compare field elements with the reference's
/// equality (`modcmp`, ported later), never by raw-limb equality.
///
/// Because Montgomery multiplication carries the inverse of `R` through,
/// `fp_mul(a, b)` is not the positional product `a * b mod p`; it is
/// `a * b * R^-1 mod p`. Two ways the Montgomery domain manifests:
/// - `fp_mul(a, MONTGOMERY_ONE)` is the canonical reduction (`redc`) of
///   `a`. `redc` is not yet ported, but this identity is the basis of
///   `nres`/`redc` once they land.
/// - The Montgomery `R^2 mod p` constant is what `nres` multiplies an
///   ordinary residue by to enter the Montgomery domain.
pub fn fp_mul(out: &mut Fp, a: &Fp, b: &Fp) {
    modmul(a, b, out);
}
