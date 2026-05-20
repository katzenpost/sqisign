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
//! - [`fp_is_zero`] is `fp_is_zero(a)`, the predicate boundary returning
//!   a `uint32_t` mask (`0xFFFFFFFF` for zero, `0` for nonzero). It is
//!   the thin wrapper the reference defines over `modis0`, which first
//!   `redc`s its argument to the canonical representative (an ordinary
//!   non-Montgomery limb vector with limbs `< 2^51` and positional value
//!   `< p`) and then OR-folds the five limbs into a single `spint d`;
//!   the bit-twiddle `(d - 1) >> 51 & 1` returns `1` when `d == 0`
//!   (`0 - 1` wraps to `0xFFFF_FFFF_FFFF_FFFF`, whose 51-bit-right shift
//!   has bit 0 set) and `0` otherwise (canonical `d` is bounded by the
//!   OR of canonical limbs each below `2^51`, so `d < 2^51` and
//!   `d - 1 < 2^51`, whose 51-bit shift is zero). The `-(uint32_t)`
//!   wrapper turns the `{0, 1}` result into the mask `{0, 0xFFFFFFFF}`
//!   the rest of the codebase consumes (e.g. `fp_select(d, _, _, ctl)`
//!   takes exactly such a mask). The Montgomery reduction path is the
//!   first appearance of `redc`/`modfsb`/`flatten`, which are ported as
//!   internal helpers here even though `redc`/`modfsb` themselves have
//!   no public boundary yet: `redc` is `modmul(n, [1, 0, 0, 0, 0]) +
//!   modfsb`, which inverts the Montgomery factor (and so reads off the
//!   ordinary residue of the n-residue), and `modfsb` (followed
//!   internally by `flatten`) performs the conditional final subtraction
//!   of `p` that brings the redundant `[2p)` representative back into
//!   `[0, p)`. The differential boundary records arbitrary five-limb
//!   inputs (the C-derived battery deliberately includes non-canonical
//!   limbs and several non-canonical encodings of the field zero); on
//!   those, the reference's `modis0` does fire `redc` first, so the port
//!   must reproduce the full Montgomery reduction chain bit-for-bit at
//!   the predicate level. No upstream defect observed; every committed
//!   C-derived vector replays bit-for-bit.
//! - [`fp_sqr`] is `fp_sqr(out, a)`, the thin wrapper the reference
//!   defines over `modsqr`: Montgomery modular squaring on the level-1
//!   generic field. `modsqr` is the squaring specialisation of `modmul`:
//!   it uses a per-column accumulator `tot`, exploits the symmetry
//!   `a[i] * a[j] == a[j] * a[i]` to compute each off-diagonal product
//!   *once* and double it (`tot *= 2;`) before folding into the running
//!   128-bit accumulator `t`, then adds the diagonal `a[i] * a[i]` (when
//!   present in the column) un-doubled. The Montgomery reduction
//!   (`v_{k-4} * p4` folds at columns 4..=8) and the unmasked limb-4
//!   final write are unchanged from `modmul`. The transcription preserves
//!   the reference's `tot = / tot += / tot *= 2 / t += tot / t += v_k*p4`
//!   per-column structure identifier-for-identifier; `tot *= 2;` is
//!   ported as `tot.wrapping_mul(2)` for the same reason `modmul` uses
//!   `wrapping_add` uniformly on `u128` (the C `__uint128_t` wraps
//!   silently, Rust panics in debug, the differential battery feeds
//!   non-canonical full-width limbs). The empirical equivalence
//!   `fp_sqr(a) == fp_mul(a, a)` is checked bit-exactly across the full
//!   1012-vector battery and pinned in `fp_sqr_props.rs`: algebraically
//!   the column sums are equal and `u128` wrapping_add is associative
//!   and commutative, so the running `t` is bit-equal at every masking
//!   point. No upstream defect observed; every committed C-derived vector
//!   replays bit-for-bit.
//! - [`fp_half`] is `fp_half(out, a)`, the thin wrapper the reference
//!   defines over `modmul`, multiplying by the precomputed Montgomery
//!   representative of `2^-1 mod p`. The reference literally is
//!   `modmul(TWO_INV, *a, *out);` (see `fp_p5248_64.c:646..650`); the
//!   port mirrors that one-line call site exactly, dispatching through
//!   the already-ported [`fp_mul`]'s `modmul` core. `TWO_INV` is the
//!   precomputed `extern const`-style limb table from
//!   `fp_p5248_64.c:532..537`, transcribed bit-for-bit into the internal
//!   constant [`TWO_INV`], the Montgomery representative of `2^-1 mod p`.
//!   The output is, as for [`fp_mul`], a redundant `[0, 2p)`
//!   representative (limbs 0..=3 below `2^51`, limb 4 unmasked); compare
//!   with [`fp_is_equal`] (the value-level equality), not raw-limb
//!   equality. The empirical equivalence `fp_half(a) == fp_mul(&TWO_INV,
//!   a)` is bit-exact by construction (`fp_half` *is* exactly that call)
//!   and is checked at the property boundary as a sanity oracle. The
//!   value-level identity `2 * fp_half(a) == a mod p` is also pinned via
//!   [`fp_is_equal`], turning the just-landed equality predicate into an
//!   independent value-level oracle for the halving operation. No
//!   upstream defect observed; every committed C-derived vector replays
//!   bit-for-bit.
//! - [`fp_div3`] is `fp_div3(out, a)`, the direct analogue of
//!   [`fp_half`] one slot over: the reference's
//!   `modmul(THREE_INV, *a, *out);` (see `fp_p5248_64.c:658..662`),
//!   multiplying by the precomputed Montgomery representative of
//!   `3^-1 mod p`. The port is the same one-line `modmul` call site
//!   substituting [`THREE_INV`] for [`TWO_INV`], with no new arithmetic
//!   introduced. `THREE_INV` is the precomputed `extern const`-style
//!   limb table from `fp_p5248_64.c:538..542`, transcribed bit-for-bit
//!   into the internal constant alongside [`TWO_INV`]. The output is,
//!   as for [`fp_mul`] and [`fp_half`], a redundant `[0, 2p)`
//!   representative (limbs 0..=3 below `2^51`, limb 4 unmasked); compare
//!   with [`fp_is_equal`] (the value-level equality), not raw-limb
//!   equality. The empirical equivalence `fp_div3(a) == fp_mul(&THREE_INV,
//!   a)` is bit-exact by construction (`fp_div3` *is* exactly that
//!   call) and is checked at the property boundary as a sanity oracle.
//!   The value-level identity `3 * fp_div3(a) == a mod p` is also
//!   pinned via [`fp_is_equal`] (tripling the third back; the redundant
//!   form is sound under the equality predicate because [`fp_is_equal`]
//!   `redc`s both sides to canonical), independently exercising
//!   `THREE_INV`'s defining property without leaning on raw-limb
//!   readings of the redundant Montgomery form. No upstream defect
//!   observed; every committed C-derived vector replays bit-for-bit.
//! - [`fp_set_small`] is `fp_set_small(x, val)`, the value setter the
//!   reference defines as `modint((int)val, *x)`. It is the third gf
//!   setter (after [`fp_set_zero`] and [`fp_set_one`]) and the first one
//!   to take a non-`fp_t` argument, which forces the harness's first
//!   setter-with-a-value record shape: `inputs := {prefill, val}`,
//!   `outputs := {out}`. The port introduces three internal helpers
//!   that the reference uses pervasively for Montgomery-domain
//!   construction: the precomputed n-residue conversion factor [`NRES_C`],
//!   the positional-to-Montgomery conversion `nres` (one `modmul` call
//!   with [`NRES_C`] as the second operand), and the int-to-Montgomery
//!   setter `modint`. The boundary itself is the one-liner
//!   `modint(val as i32, out)`: the C wrapper narrows its `digit_t`
//!   argument to `int` before calling `modint`, so only the low 32 bits
//!   of `val` are observable through the boundary, and the result
//!   `a[0] = (spint)(int32_t)val` sign-extends to `u64` for values above
//!   `2^31 - 1` (the port reproduces this as `val as i32` then
//!   `x as i64 as u64` for the limb-0 write, the explicit two-step cast
//!   that mirrors the C cast chain bit for bit). The fp_t output is the
//!   Montgomery representative of the narrowed and sign-extended
//!   integer; when `val` is `1` the output is the same
//!   [`MONTGOMERY_ONE`] [`fp_set_one`] writes directly, which exercises
//!   `nres` through the boundary and confirms [`NRES_C`] is correct
//!   algorithmically (an independent check beyond the bit-for-bit
//!   transcription) -- pinned by the unit test
//!   `nres_of_positional_one_is_montgomery_one` and re-pinned at the
//!   property level by `fp_set_small(1) == fp_set_one()` for every
//!   pre-fill. The high-bits-ignored narrowing (`val` differs from
//!   `val as i32 as u64` for any `val` outside `[-2^31, 2^31 - 1]`) is
//!   pinned in 1048 of the 1132 differential records and re-pinned at
//!   the property level by `fp_set_small(val) ==
//!   fp_set_small(val as i32 as u64)`. No upstream defect observed;
//!   every committed C-derived vector replays bit-for-bit.
//! - [`fp_mul_small`] is `fp_mul_small(out, a, val)`, the binary-mixed
//!   boundary (one `fp_t` input plus one `uint32_t` scalar) the reference
//!   defines over `modmli`. `modmli(a, b, c)` is itself the thin
//!   two-line combinator that builds the Montgomery representative of
//!   the integer `b` via [`modint`] (into a five-limb scratch) and then
//!   runs a single [`modmul`] of `a` against that scratch into `c`. No
//!   new arithmetic is introduced at any layer: [`modint`], [`nres`],
//!   [`NRES_C`] and [`modmul`] were all landed with the earlier ports;
//!   [`modmli`] is a pure combinator. The C wrapper narrows the
//!   `uint32_t val` to `int` before calling `modmli`, so only the low 32
//!   bits of `val` are observable through the boundary and the result
//!   sign-extends through the same chain `fp_set_small` documents (the
//!   port reproduces this as `val as i32` at the public-wrapper
//!   call-site). Because the boundary takes one `fp_t` input and one
//!   `u64`-encoded scalar (the cdump emitter widens the recorded `val`
//!   to 8 little-endian bytes for shape uniformity with
//!   `fp_set_small`'s recording, even though the C signature is
//!   `uint32_t`), the cdump harness adds a new emitter shape
//!   `emit_fp_mul_val` (`{a, val} -> {out}`) alongside the
//!   already-landed `emit_fp_value_setter`'s `{prefill, val} -> {out}`.
//!   The differential battery is 132 edges (12 fp patterns x 11 vals)
//!   plus 1000 sweep seeds, 1132 records total, the same shape
//!   `fp_set_small`'s battery uses. Three sound raw-limb properties
//!   are pinned: `val == 0` yields the canonical all-zero limb vector
//!   for arbitrary `a` (every cross-product column sum is zero, the
//!   final `c[4] = (spint)t` truncation reads off a zero `t`); `val ==
//!   1` yields the same redundant representative `fp_mul(a, &nres(1))`
//!   does (the boundary is literally that call, threaded through the
//!   `modmli -> modint(1) -> modmul` chain that produces
//!   `MONTGOMERY_ONE` at the scratch buffer); and cross-oracle
//!   `fp_mul_small(out, a, val) == fp_mul(a, &fp_set_small(val))`
//!   bit-exact for arbitrary inputs (the just-landed `fp_set_small +
//!   fp_mul` re-express the same `modint + modmul` chain `modmli`
//!   evaluates, so the equivalence is a structural identity verified
//!   empirically across all 1132 records before pinning). No upstream
//!   defect observed; every committed C-derived vector replays
//!   bit-for-bit.
//! - [`fp_exp3div4`] is `fp_exp3div4(out, a)`, the thin wrapper the
//!   reference defines over `modpro`:
//!
//!   ```c
//!   void fp_exp3div4(fp_t *out, const fp_t *a) {
//!       modpro(*a, *out);
//!   }
//!   ```
//!
//!   `modpro` is the level-1 field's hand-built fixed addition chain
//!   that computes the *progenitor* `a^((p-3)/4) mod p`, the Montgomery
//!   representative the rest of the chain ([`fp_inv`], [`fp_sqrt`],
//!   [`fp_is_square`]) builds on by squaring further and/or multiplying
//!   by `a`. The port adds two new internal helpers alongside the
//!   existing `modmul`/`modsqr`/`modcpy` infrastructure: [`modnsqr`],
//!   the trivial n-fold-squaring loop, and [`modpro`] itself, a
//!   step-for-step transcription of the reference's fixed chain of
//!   `modsqr`/`modmul`/`modnsqr` calls preserving the six scratch
//!   buffers (`x`, `t0..t4`) the reference uses. Bit-for-bit
//!   correspondence is established by the differential battery (every
//!   committed C-derived vector replays bit-for-bit at the redundant
//!   limb level); no upstream defect observed.
//! - [`fp_inv`] is `fp_inv(x)`, the in-place modular inverse the
//!   reference defines over `modinv`:
//!
//!   ```c
//!   void fp_inv(fp_t *x) {
//!       modinv(*x, NULL, *x);
//!   }
//!   ```
//!
//!   `modinv` builds the inverse via Fermat (`x^(p-2)`) by computing the
//!   progenitor `x^((p-3)/4)` (when `h == NULL`, as it is here), squaring
//!   it twice to obtain `x^(p-3)`, and multiplying by `x` to obtain
//!   `x^(p-2) == x^-1`. The port adds [`modinv`] as a new internal
//!   helper threading [`modpro`]/[`modnsqr`]/[`modmul`]/[`modcpy`]
//!   (already ported with the [`fp_exp3div4`] commit, modulo modinv
//!   itself); the in-place wrapper resolves Rust's borrow-checker
//!   conflict by snapshotting `*x` to a local before passing the
//!   snapshot as input and `x` as destination. The differential battery
//!   exercises the wrapper's `None` branch; the `Some(h)` branch is
//!   exercised by the cross-validating property test that recomputes
//!   the progenitor and inverts via the precomputed-progenitor path.
//!   The value-level identity `fp_mul(fp_inv(a), a) ==_field 1` is
//!   pinned as a property test (the named exception per sir's
//!   directive, sound on canonical-nonzero inputs). No upstream defect
//!   observed.
//! - [`fp_is_square`] is `fp_is_square(a)`, the predicate boundary
//!   returning a `uint32_t` mask (`0xFFFFFFFF` if `a` is a quadratic
//!   residue mod `p`, `0` otherwise). The reference's wrapper is
//!
//!   ```c
//!   uint32_t fp_is_square(const fp_t *a) {
//!       return -(uint32_t)modqr(NULL, *a);
//!   }
//!   ```
//!
//!   `modqr` evaluates the Euler criterion `x^((p-1)/2)`: when `h ==
//!   NULL` it computes the progenitor `x^((p-3)/4)`, squares it to get
//!   `x^((p-3)/2)`, multiplies by `x` to get `x^((p-1)/2)`, and tests
//!   for unity via [`modis1`] OR-d with [`modis0`] (so the field zero
//!   is treated as a square per the reference's convention). The
//!   `-(uint32_t)` cast turns the `{0, 1}` return into the
//!   `{0, 0xFFFFFFFF}` mask, the same shape as [`fp_is_zero`] and
//!   [`fp_is_equal`]. New internal helpers landed here: [`modis1`] (the
//!   is-Montgomery-one predicate, twin to [`modis0`] using the same
//!   `(d - 1) >> 51 & 1` zero-detect trick) and [`modqr`] itself. The
//!   differential battery includes the canonical zero (positive
//!   outcome by convention), the Montgomery [`MONTGOMERY_ONE`]
//!   (positive: `1` is a square), the radix-2^51 encoding of `p`
//!   (positive: reduces to zero), and arbitrary full-width limbs
//!   (mostly negative: roughly half the residues are quadratic
//!   non-residues). No upstream defect observed.
//! - [`fp_sqrt`] is `fp_sqrt(a)`, the in-place modular square root the
//!   reference defines over `modsqrt`:
//!
//!   ```c
//!   void fp_sqrt(fp_t *a) {
//!       modsqrt(*a, NULL, *a);
//!   }
//!   ```
//!
//!   `modsqrt` exploits `p == 3 mod 4` (which holds for the level-1
//!   `p5248`): `sqrt(x) == x^((p+1)/4) == x^((p-3)/4) * x ==
//!   progenitor * x mod p`. When `h == NULL` (as it is at this boundary)
//!   the progenitor is computed via [`modpro`] and the single [`modmul`]
//!   yields the root. On a non-residue, the returned value is
//!   meaningless garbage (the reference makes no defensive check; the
//!   port follows). [`modsqrt`] is added as a new internal helper. The
//!   in-place wrapper resolves the borrow-checker conflict as [`fp_inv`]
//!   does. The value-level identity `fp_sqr(fp_sqrt(a)) ==_field a`
//!   on quadratic residues is pinned as the named property exception
//!   (gated on [`fp_is_square`] returning the positive mask). No upstream
//!   defect observed.
//! - [`fp_encode`] is `fp_encode(dst, a)`, the 32-byte canonical
//!   little-endian serialization the reference defines as a modified
//!   `modexp`:
//!
//!   ```c
//!   void fp_encode(void *dst, const fp_t *a) {
//!       spint c[5];
//!       redc(*a, c);
//!       for (int i = 0; i < 32; i++) {
//!           ((char *)dst)[i] = c[0] & 0xff;
//!           (void)modshr(8, c);
//!       }
//!   }
//!   ```
//!
//!   [`redc`] canonicalises the Montgomery representative to the
//!   positional residue below `p`, then the byte loop peels off the low
//!   byte of limb 0 and shifts the five-limb value right by 8 bits
//!   ([`modshr`], a new internal helper) thirty-two times. The port
//!   takes `&mut [u8; 32]` for the destination so the 32-byte contract
//!   is encoded in the type. The boundary's record schema is the new
//!   `fp -> bytes` shape: `{"a": <5 u64 LE>} -> {"dst": <32 byte LE>}`.
//!   No upstream defect observed.
//! - [`fp_decode`] is `fp_decode(d, src)`, the canonical-range-checked
//!   32-byte deserialization the reference defines as a modified
//!   `modimp`:
//!
//!   ```c
//!   uint32_t fp_decode(fp_t *d, const void *src) {
//!       const unsigned char *b = src;
//!       for (int i = 0; i < 5; i++) (*d)[i] = 0;
//!       for (int i = 31; i >= 0; i--) {
//!           modshl(8, *d);
//!           (*d)[0] += (spint)b[i];
//!       }
//!       spint res = (spint)-modfsb(*d);
//!       nres(*d, *d);
//!       for (int i = 0; i < 5; i++) (*d)[i] &= res;
//!       return (uint32_t)res;
//!   }
//!   ```
//!
//!   The bytes are folded into `d` in descending address order via
//!   [`modshl`] (the new internal helper, in-place left shift by `n <
//!   51` bits) and a single-byte add into limb 0. [`modfsb`] returns
//!   `1` iff the decoded value is below `p`; the reference negates that
//!   to a full-width mask `res`, runs [`nres`] to convert to Montgomery
//!   form, and ANDs `res` into every limb so an out-of-range input is
//!   zeroed on the way out. The returned `uint32_t` mask is the same
//!   shape as the other predicate wrappers: `0xFFFFFFFF` on canonical
//!   in-range input, `0` on out-of-range. The boundary's record schema
//!   is the new `bytes -> fp + u32` shape: `{"src": <32 byte LE>} ->
//!   {"d": <5 u64 LE>, "result": <4 byte LE u32>}`. The differential
//!   battery deliberately partitions edge inputs into canonical
//!   (positive outcome, `d` set) and non-canonical (negative outcome,
//!   `d` zeroed) classes. No upstream defect observed.
//! - [`fp_decode_reduce`] is `fp_decode_reduce(d, src, len)`, the
//!   arbitrary-length-input reducer the reference defines as a two-phase
//!   fold:
//!
//!   ```c
//!   void fp_decode_reduce(fp_t *d, const void *src, size_t len);
//!   ```
//!
//!   The trailing partial block (`len % 32`) is decoded into `d` via
//!   [`fp_decode`] after zero-padding, then each preceding 32-byte
//!   block is partially-reduced via the level-1 prime's `5 * 2^248 ==
//!   1 mod p` identity ([`partial_reduce`], a new internal helper on a
//!   plain 4-limb 256-bit array), re-encoded, decoded via [`fp_decode`]
//!   again, and added to `d` *after* `d` has been multiplied by [`R2`]
//!   `== 2^256 mod p` (the Montgomery representative of the per-block
//!   shift). New internal helpers landed here: [`R2`], [`partial_reduce`],
//!   [`add_carry`], [`dec64le`], [`enc64le`]. The boundary's record
//!   schema is the variable-length `bytes -> fp` shape: `{"src": <hex
//!   of bytes>, "len": <8 byte LE u64>} -> {"d": <5 u64 LE>}`. The
//!   differential battery sweeps representative lengths including `0`
//!   (the empty-input early-return), `< 32` (partial-block-only),
//!   `32`/`33`/`63`/`64`/`100`/`200` (one full block crossed with
//!   partial-block variants and multi-block chains). No upstream defect
//!   observed.
//! - [`fp_is_equal`] is `fp_is_equal(a, b)`, the second predicate boundary
//!   in the gf battery: a *binary* predicate returning the same
//!   `uint32_t` mask shape as [`fp_is_zero`] (`0xFFFFFFFF` when `a` and
//!   `b` represent the same field element, `0` otherwise). It is the
//!   thin wrapper the reference defines over `modcmp`, which `redc`s
//!   **both** operands to their canonical representatives (each limb
//!   below `2^51`, positional value below `p`), then per limb applies
//!   the same `(x - 1) >> 51 & 1` zero-detect trick [`modis0`] uses to
//!   `c[i] ^ d[i]` and ANDs the five resulting bits into an `eq`
//!   accumulator initialised to `1`. The per-limb XOR is `0` iff the
//!   canonical limbs match, and the canonical-bounded property
//!   ([`redc`]'s post-condition: each limb below `2^51`) is exactly
//!   what makes the trick correct here, the same precondition [`modis0`]
//!   relies on for its single-OR-fold variant. The `-(uint32_t)`
//!   wrapper turns the inner `{0, 1}` result into the
//!   `{0, 0xFFFFFFFF}` mask the rest of the codebase consumes (the same
//!   negation used by [`fp_is_zero`]). The port reuses the
//!   `redc`/`modfsb`/`flatten` chain ported with [`fp_is_zero`] (commit
//!   `1fd71e360e6a42ee704f0436a63f91dc5916931e`) verbatim, and adds
//!   `modcmp` as a new internal helper alongside [`modis0`]; the
//!   per-limb zero-detect application is precisely the same bit
//!   pattern, just five times over `c[i] ^ d[i]` rather than once over
//!   the OR-fold. The differential battery exercises the cross product
//!   of the 12-pattern edge set (including the radix-2^51 encoding of
//!   `p` itself) plus the 1000-seed pseudo-random sweep, so multiple
//!   pairs reduce to the same canonical form via different
//!   non-canonical representatives (the `(canonical zero, p as limbs)`
//!   pair in particular: `redc` brings both to `[0, 0, 0, 0, 0]`,
//!   `fp_is_equal` returns the all-ones mask). No upstream defect
//!   observed; every committed C-derived vector replays bit-for-bit.
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

/// Montgomery modular squaring, `c = a * a mod 2p` in the redundant
/// radix-2^51 representation.
///
/// Mirrors the reference's
/// `inline static void modsqr(const spint *a, spint *c)` from
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:156..220`
/// statement-for-statement.
///
/// `modsqr` is the squaring specialisation of `modmul`. Where `modmul`
/// sums the full set of partial products `a[i] * b[j]` per column, `modsqr`
/// exploits `a[i] * a[j] == a[j] * a[i]` to compute each off-diagonal
/// product *once* and double it (in the reference, `tot *= 2;` on the
/// per-column `tot` accumulator before it is added to the running `t`); the
/// diagonal product `a[i] * a[i]` (when present in the column) is added
/// after the doubling, never doubled. The Montgomery reduction structure is
/// otherwise identical to `modmul`: at columns 4..=8 the running `t` also
/// absorbs `v_{k-4} * p4`, the inline fold of the prime's only non-zero
/// limb (`p4 == 0x500000000000`), and the final write `c[4] = (spint)t` is
/// a full 64-bit truncation of the residual accumulator (no `& mask`),
/// exactly as in `modmul`.
///
/// The transcription rules are the same as `modmul`'s:
/// - `spint`/`dpint`/`udpint` -> `u64` / `u128` / `u128`.
/// - `(udpint)a[i] * a[j]` -> `(a[i] as u128) * (a[j] as u128)`. Each
///   factor fits in `u64`, so the product fits in `u128` and plain `*` is
///   bit-exact to the reference's `(__uint128_t)a[i] * a[j]`.
/// - `tot *= 2;` is the reference's column-accumulator doubling. It is
///   *unsigned wraparound* on `__uint128_t` and the differential battery
///   feeds non-canonical full-width limbs, so the port uses
///   `tot.wrapping_mul(2)` to preserve the bit-exact wraparound rather
///   than the panicking `*= 2`. Equivalent to `tot.wrapping_shl(1)`; the
///   `wrapping_mul(2)` spelling is chosen so the identifier-by-identifier
///   correspondence with the reference's `tot *= 2;` is visible.
/// - `tot = ...`, `tot += ...`, `t = tot`, `t += tot` are all
///   `tot.wrapping_add(...)` / `t.wrapping_add(tot)` for the same reason
///   `modmul` uses `wrapping_add` uniformly on the u128 accumulator: the C
///   `__uint128_t` wraps silently, Rust `u128 += ...` panics in debug, and
///   the differential battery feeds non-canonical full-width limbs that
///   can push the accumulator past `2^128`. Wrapping_add is associative
///   and commutative on `u128` modulo `2^128`, so the column-by-column
///   `t` value at every masking point is the same as the reference's.
/// - `(spint)t & mask` -> `(t as u64) & MASK51`.
/// - `(udpint)v_k * p4` -> `(v_k as u128) * (P4 as u128)`, plain `*` (the
///   factors are bounded by `2^51` and `2^47` respectively, so no
///   overflow). Folded into `t` at columns 4..=8, exactly as in `modmul`.
/// - `t >>= 51` applies directly to `u128`.
///
/// Identifier mapping (each row preserves the reference's name verbatim):
/// `tot`, `t`, `v0`, `v1`, `v2`, `v3`, `v4`, `c[0..=4]`, `p4` (constant
/// inlined as `P4 as u128` at each fold). The statement order, the
/// per-column count of `tot = / tot += / tot *= 2 / t += tot /
/// t += v_k*p4 / t >>= 51 / mask write` constructs, and the precise
/// placement of every `tot *= 2;` (always *after* the off-diagonal `tot`
/// builds and *before* any same-index diagonal `tot += a[i]*a[i]` add)
/// are preserved so the column-by-column trace lines up one-to-one with
/// `fp_p5248_64.c:156..220`.
fn modsqr(a: &Fp, c: &mut Fp) {
    // The reference declares `udpint t = 0;` and overwrites it on the very
    // next statement (`t = tot;`); the initial zero is dead-store by
    // construction. The port collapses the dead-store into the
    // declaration's initial assignment from `tot` so clippy's
    // `unused_assignments` lint is honoured without altering the
    // statement-for-statement column-by-column trace below.
    let mut tot: u128 = (a[0] as u128) * (a[0] as u128);
    let mut t: u128 = tot;
    let v0: u64 = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[0] as u128) * (a[1] as u128);
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    let v1: u64 = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[0] as u128) * (a[2] as u128);
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[1] as u128) * (a[1] as u128));
    t = t.wrapping_add(tot);
    let v2: u64 = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[0] as u128) * (a[3] as u128);
    tot = tot.wrapping_add((a[1] as u128) * (a[2] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    let v3: u64 = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[0] as u128) * (a[4] as u128);
    tot = tot.wrapping_add((a[1] as u128) * (a[3] as u128));
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[2] as u128) * (a[2] as u128));
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v0 as u128) * (P4 as u128));
    let v4: u64 = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[1] as u128) * (a[4] as u128);
    tot = tot.wrapping_add((a[2] as u128) * (a[3] as u128));
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v1 as u128) * (P4 as u128));
    c[0] = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[2] as u128) * (a[4] as u128);
    tot = tot.wrapping_mul(2);
    tot = tot.wrapping_add((a[3] as u128) * (a[3] as u128));
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v2 as u128) * (P4 as u128));
    c[1] = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[3] as u128) * (a[4] as u128);
    tot = tot.wrapping_mul(2);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v3 as u128) * (P4 as u128));
    c[2] = (t as u64) & MASK51;
    t >>= 51;
    tot = (a[4] as u128) * (a[4] as u128);
    t = t.wrapping_add(tot);
    t = t.wrapping_add((v4 as u128) * (P4 as u128));
    c[3] = (t as u64) & MASK51;
    t >>= 51;
    c[4] = t as u64;
}

/// GF(p) Montgomery squaring `out = a * a * R^-1 mod p`, in the redundant
/// radix-2^51 representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_sqr(fp_t *out, const fp_t *a)`, which
/// is the thin wrapper `modsqr(*a, *out)`. As with [`fp_mul`], the output
/// is *not* fully canonical: limbs 0..=3 are below `2^51` (the column mask)
/// but limb 4 is left unmasked (the final `c[4] = (spint)t` is a full
/// 64-bit truncation of the residual accumulator, no `& mask`), exactly as
/// the reference leaves it. Compare field elements with the reference's
/// equality (`modcmp`, ported later), never by raw-limb equality.
///
/// `modsqr` is the squaring specialisation of `modmul`: the per-column set
/// of partial products `{ a[i] * a[j] : i + j == k }` is symmetric under
/// `(i, j) -> (j, i)`, so each off-diagonal product is computed once and
/// doubled (`tot *= 2;`) before being added to the running accumulator,
/// rather than computed twice; the diagonal product `a[i] * a[i]` (when
/// present in the column) is added after the doubling. The Montgomery
/// reduction (`v_{k-4} * p4` folds at columns 4..=8) and the unmasked
/// limb-4 write are unchanged. Because the column sum is algebraically the
/// same value as `modmul(a, a)`'s, `fp_sqr(a)` and `fp_mul(a, a)` produce
/// the bit-exact same `fp_t` output: the order of the `t.wrapping_add`
/// terms differs but `wrapping_add` is associative and commutative on
/// `u128` modulo `2^128`, so the running `t` is bit-equal at every
/// masking point. This equivalence was checked empirically across the
/// full differential battery before being pinned in
/// `fp_sqr_props.rs`.
pub fn fp_sqr(out: &mut Fp, a: &Fp) {
    modsqr(a, out);
}

/// Montgomery representative of `2^-1 mod p` on the level-1 generic
/// field: the precomputed constant `2^-1 * R mod p` in the unsaturated
/// radix-2^51 limb layout, transcribed verbatim from the reference's
/// `static const digit_t TWO_INV[NWORDS_FIELD]` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:532..536`:
///
/// ```c
/// // Montgomery representation of 2^-1
/// static const digit_t TWO_INV[NWORDS_FIELD] = { 0x000000000000000c,
///                                                0x0000000000000000,
///                                                0x0000000000000000,
///                                                0x0000000000000000,
///                                                0x0000400000000000 };
/// ```
///
/// Used by [`fp_half`], the thin wrapper `modmul(TWO_INV, *a, *out)`. The
/// constant is taken verbatim; no derivation is performed at the port
/// (deriving it would require the Montgomery `R^2` constant and `nres`,
/// neither of which is ported yet). The bit-for-bit correspondence with
/// the reference is the load-bearing property; the value-level identity
/// `2 * (TWO_INV * a * R^-1) == a mod p` is checked at the property
/// boundary via [`fp_is_equal`].
const TWO_INV: Fp = [
    0x0000_0000_0000_000c,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_4000_0000_0000,
];

/// GF(p) Montgomery halving `out = a / 2 mod p`, in the redundant
/// radix-2^51 representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_half(fp_t *out, const fp_t *a)`
/// (`vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:646..650`), which
/// is the one-liner `modmul(TWO_INV, *a, *out);`. The port is the same
/// one-liner: a single `modmul` call with the constant [`TWO_INV`] (the
/// Montgomery representative of `2^-1 mod p`) as the first operand and
/// `a` as the second, threading through the already-ported [`fp_mul`]'s
/// `modmul` core. No new arithmetic of any kind is introduced.
///
/// Two faithfully reproduced subtleties from the underlying [`fp_mul`]
/// path:
/// 1. The output is *not* fully canonical: limbs 0..=3 are below `2^51`
///    (the column mask) but limb 4 is left unmasked (the final
///    `c[4] = (spint)t` in `modmul` is a full 64-bit truncation of the
///    residual accumulator, no `& mask`), exactly as the reference
///    leaves it. Compare field elements with the value-level
///    [`fp_is_equal`], never by raw-limb equality.
/// 2. The Montgomery domain is preserved: if `a == A * R mod p`
///    positionally, then `modmul(TWO_INV, a)` computes
///    `(2^-1 * R) * (A * R) * R^-1 == (A / 2) * R mod p`. The output is
///    the Montgomery representative of `A / 2`, ready to be consumed by
///    further Montgomery-domain operations without re-conversion.
///
/// `TWO_INV` is taken **verbatim** from the reference's precomputed
/// limb table at `fp_p5248_64.c:532..536`; the port does not derive it
/// (deriving it would require `nres` and the Montgomery `R^2`
/// constant, neither of which is ported). The bit-for-bit
/// correspondence with the reference's constant is the load-bearing
/// property; the value-level identity `2 * fp_half(a) == a mod p` is
/// pinned at the property boundary via [`fp_is_equal`] (now that the
/// equality predicate is ported, the doubling-back oracle is sound on
/// the canonical-equality domain even though both sides are redundant
/// limb vectors).
pub fn fp_half(out: &mut Fp, a: &Fp) {
    modmul(&TWO_INV, a, out);
}

/// Montgomery representative of `3^-1 mod p` on the level-1 generic
/// field: the precomputed constant `3^-1 * R mod p` in the unsaturated
/// radix-2^51 limb layout, transcribed verbatim from the reference's
/// `static const digit_t THREE_INV[NWORDS_FIELD]` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:538..542`:
///
/// ```c
/// // Montgomery representation of 3^-1
/// static const digit_t THREE_INV[NWORDS_FIELD] = { 0x000555555555555d,
///                                                  0x0002aaaaaaaaaaaa,
///                                                  0x0005555555555555,
///                                                  0x0002aaaaaaaaaaaa,
///                                                  0x0000455555555555 };
/// ```
///
/// Used by [`fp_div3`], the thin wrapper `modmul(THREE_INV, *a, *out)`.
/// The constant is taken verbatim; no derivation is performed at the
/// port (deriving it would require the Montgomery `R^2` constant and
/// `nres`, neither of which is ported yet). The bit-for-bit
/// correspondence with the reference is the load-bearing property; the
/// value-level identity `3 * (THREE_INV * a * R^-1) == a mod p` is
/// checked at the property boundary via [`fp_is_equal`].
const THREE_INV: Fp = [
    0x0005_5555_5555_555d,
    0x0002_aaaa_aaaa_aaaa,
    0x0005_5555_5555_5555,
    0x0002_aaaa_aaaa_aaaa,
    0x0000_4555_5555_5555,
];

/// GF(p) Montgomery division-by-three `out = a / 3 mod p`, in the
/// redundant radix-2^51 representation, reduced to less than `2p`.
///
/// Mirrors the reference's `void fp_div3(fp_t *out, const fp_t *a)`
/// (`vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:658..662`), which
/// is the one-liner `modmul(THREE_INV, *a, *out);`. The port is the
/// same one-liner: a single `modmul` call with the constant
/// [`THREE_INV`] (the Montgomery representative of `3^-1 mod p`) as the
/// first operand and `a` as the second, threading through the
/// already-ported [`fp_mul`]'s `modmul` core. No new arithmetic of any
/// kind is introduced; the function is the direct analogue of
/// [`fp_half`] one slot over with [`THREE_INV`] substituted for
/// [`TWO_INV`].
///
/// Two faithfully reproduced subtleties from the underlying [`fp_mul`]
/// path:
/// 1. The output is *not* fully canonical: limbs 0..=3 are below `2^51`
///    (the column mask) but limb 4 is left unmasked (the final
///    `c[4] = (spint)t` in `modmul` is a full 64-bit truncation of the
///    residual accumulator, no `& mask`), exactly as the reference
///    leaves it. Compare field elements with the value-level
///    [`fp_is_equal`], never by raw-limb equality.
/// 2. The Montgomery domain is preserved: if `a == A * R mod p`
///    positionally, then `modmul(THREE_INV, a)` computes
///    `(3^-1 * R) * (A * R) * R^-1 == (A / 3) * R mod p`. The output
///    is the Montgomery representative of `A / 3`, ready to be consumed
///    by further Montgomery-domain operations without re-conversion.
///
/// `THREE_INV` is taken **verbatim** from the reference's precomputed
/// limb table at `fp_p5248_64.c:538..542`; the port does not derive it
/// (deriving it would require `nres` and the Montgomery `R^2`
/// constant, neither of which is ported). The bit-for-bit
/// correspondence with the reference's constant is the load-bearing
/// property; the value-level identity `3 * fp_div3(a) == a mod p` is
/// pinned at the property boundary via [`fp_is_equal`] (the equality
/// predicate's `redc` of both sides makes the tripling-back oracle
/// sound on the canonical-equality domain even though both sides are
/// redundant limb vectors).
pub fn fp_div3(out: &mut Fp, a: &Fp) {
    modmul(&THREE_INV, a, out);
}

/// Propagate carries and, if `prop` signalled the value went negative,
/// add `p` back limbwise; propagate carries once more. Returns `1` if the
/// correction fired, else `0`.
///
/// Mirrors the reference's
/// `inline static int flatten(spint *n)`:
///
/// ```c
/// spint carry = prop(n);
/// n[0] -= (spint)1u & carry;
/// n[4] += ((spint)0x500000000000u) & carry;
/// (void)prop(n);
/// return (int)(carry & 1);
/// ```
///
/// Composition with [`modfsb`] is the load-bearing thing to understand:
/// `modfsb` pre-adds `+1` at limb 0 and `-p4` at limb 4 (effectively a
/// trial `-p` in the redundant representation, since the level-1 prime
/// `p5248 = 5 * 2^248 - 1` has only `p4` non-zero in the radix-2^51
/// layout and the `-1` constant is the `-(-1)` from `p`'s low term), then
/// calls `flatten`. `flatten` then runs `prop`, which returns an all-ones
/// `carry` mask iff the value went negative; on that branch, `flatten`
/// adds the just-subtracted `p` back limbwise (the `+1 & carry` at limb 0
/// and `+p4 & carry` at limb 4) and propagates carries once more. The
/// returned `(int)(carry & 1)` is `1` exactly when the trial subtraction
/// was undone (i.e. the input was already below `p`); the reference's
/// `redc` discards this return value (`(void)modfsb(m);`), and the port
/// preserves the return for symmetry even though [`redc`] does the same.
///
/// All arithmetic is unsigned `u64` wraparound, exactly as in the
/// reference's `spint == uint64_t`: `n[0] -= 1 & carry` is
/// `n[0].wrapping_sub(1 & carry)` and so on. The `& carry` mask is `0`
/// or all-ones, so the corrections are either no-ops or full-width adds
/// of the prime's limb contributions.
fn flatten(n: &mut Fp) -> i32 {
    let carry = prop(n);
    n[0] = n[0].wrapping_sub(1 & carry);
    n[4] = n[4].wrapping_add(P4 & carry);
    let _ = prop(n);
    (carry & 1) as i32
}

/// Montgomery final subtract: trial-subtract `p` limbwise, then [`flatten`]
/// (which propagates and conditionally adds `p` back if the value went
/// negative). Returns `1` if the final value is below `p` already (the
/// correction undid the trial subtraction), else `0`.
///
/// Mirrors the reference's
/// `inline static int modfsb(spint *n)`:
///
/// ```c
/// n[0] += (spint)1u;
/// n[4] -= (spint)0x500000000000u;
/// return flatten(n);
/// ```
///
/// The `+1` at limb 0 and `-p4` at limb 4 together form the trial
/// `-p` in the redundant radix-2^51 form: the level-1 prime
/// `p5248 = 5 * 2^248 - 1` decomposes as `+p4 << (51*4) - 1` in this
/// layout, so subtracting it is `-(+1) + (-p4)` at limbs 0 and 4
/// respectively. [`flatten`] then runs `prop`; if the running value went
/// negative under the trial subtraction (the `prop` carry mask is all
/// ones), [`flatten`] adds `p` back limbwise (`+1` at limb 0, `+p4` at
/// limb 4) and propagates once more. The cumulative effect is: if the
/// input was `>= p`, leave the canonical reduced value below `p`; if it
/// was already below `p`, restore the original. This is exactly the
/// behaviour [`redc`] needs to canonicalise the Montgomery output of
/// `modmul`'s redundant `[0, 2p)` representative.
///
/// Wrapping `u64` arithmetic throughout matches the reference's
/// `spint == uint64_t` two's-complement wraparound (the `-p4` underflows
/// silently and is fixed up by the subsequent `+p4` on the
/// negative-branch).
fn modfsb(n: &mut Fp) -> i32 {
    n[0] = n[0].wrapping_add(1);
    n[4] = n[4].wrapping_sub(P4);
    flatten(n)
}

/// Convert an n-residue back to its ordinary residue: `m = n * R^-1 mod p`,
/// fully canonical (`m[0..=3] < 2^51`, positional value `< p`).
///
/// Mirrors the reference's
/// `static void redc(const spint *n, spint *m)`:
///
/// ```c
/// spint c[5] = {1, 0, 0, 0, 0};
/// modmul(n, c, m);
/// (void)modfsb(m);
/// ```
///
/// The trick: `modmul(n, [1, 0, 0, 0, 0])` is exactly the Montgomery
/// reduction `n * 1 * R^-1 mod p` (multiplication by positional `1`
/// followed by `modmul`'s inline Montgomery reduction), which leaves a
/// redundant `[0, 2p)` representative of the ordinary residue. [`modfsb`]
/// then performs the conditional final subtraction of `p` to bring the
/// result fully into `[0, p)`. The `(void)modfsb(m);` discards the
/// `flatten` return value, and the port follows: the return is ignored
/// here, used only for the return-value contract of [`modfsb`]/[`flatten`]
/// themselves (callers outside [`redc`] that need it can read it).
fn redc(n: &Fp, m: &mut Fp) {
    let c: Fp = [1, 0, 0, 0, 0];
    modmul(n, &c, m);
    let _ = modfsb(m);
}

/// Test whether the field element represented by `a` is zero. Returns `1`
/// if zero, else `0` (a plain integer; [`fp_is_zero`] negates it to the
/// `0xFFFFFFFF`/`0` mask the rest of the codebase consumes).
///
/// Mirrors the reference's
/// `static int modis0(const spint *a)`:
///
/// ```c
/// spint c[5];
/// spint d = 0;
/// redc(a, c);
/// for (i = 0; i < 5; i++) d |= c[i];
/// return ((spint)1 & ((d - (spint)1) >> 51u));
/// ```
///
/// Two faithfully reproduced subtleties:
/// 1. `redc(a, c)` is not a no-op: the differential boundary records
///    non-canonical limb encodings of the field zero as well as the
///    canonical one (the C-derived battery includes patterns like
///    `[mask51, mask51, mask51, mask51, p4]`, several representatives of
///    `0 mod p` in the redundant form, and arbitrary full-width limbs).
///    For any such representative, `redc` returns the canonical reduced
///    value (with `c[0..=3] < 2^51`, positional value `< p`), so the OR
///    of the canonical limbs is `0` iff the input was congruent to zero.
/// 2. The `((d - 1) >> 51) & 1` bit-twiddle: when `d == 0`,
///    `d - 1 == 0xFFFF_FFFF_FFFF_FFFF` (unsigned wraparound), whose
///    51-bit right shift `>> 51u` is `0x1FFF` (the top 13 bits packed),
///    whose low bit is `1`; the masked result is `1`. When `d != 0`,
///    `d` is the OR of canonical limbs each below `2^51`, so `d < 2^51`,
///    `d - 1 < 2^51 - 1` is non-negative, and `(d - 1) >> 51 == 0`, so
///    the masked result is `0`. The trick is faithful in `u64` wrapping
///    arithmetic for *canonical* `d`; this is exactly the domain the
///    `redc` precondition delivers, regardless of how non-canonical the
///    input was. (For arbitrary `u64 d`, the trick is not the same as
///    `d == 0`: e.g. `d = 1 << 51` gives `d - 1 = (1 << 51) - 1` whose
///    `>> 51` is `0`, then `& 1` is `0`, so it would report "zero" for a
///    nonzero `d`; this is fine because the canonical-bounded property
///    keeps us strictly below `2^51`.)
///
/// The return is a plain `u32` `{0, 1}` so the caller can either consume
/// the boolean directly (the way `modis1`'s wrapper ANDs two such bits)
/// or negate it to the `0xFFFFFFFF`/`0` mask the public boundary
/// returns.
fn modis0(a: &Fp) -> u32 {
    let mut c: Fp = [0u64; NWORDS_FIELD];
    redc(a, &mut c);
    let mut d: u64 = 0;
    for limb in &c {
        d |= *limb;
    }
    (1u64 & (d.wrapping_sub(1) >> 51)) as u32
}

/// GF(p) zero predicate: returns the constant-time mask `0xFFFFFFFF` if
/// `a` represents the field zero, else `0`.
///
/// Mirrors the reference's
/// `uint32_t fp_is_zero(const fp_t *a) { return -(uint32_t)modis0(*a); }`
/// from `fp_p5248_64.c:581..584`. [`modis0`] returns `0` or `1`; the
/// `-(uint32_t)` cast turns that into the `0xFFFFFFFF`/`0` mask
/// downstream consumers (`fp_select`, `fp_cswap`'s LSB, the
/// `fp2_is_zero` AND-of-two-fp_is_zeros) expect.
///
/// The cast chain `-(uint32_t)int01` in C is: widen the `int` `{0, 1}` to
/// `uint32_t` (no-op for non-negative values), then negate as `uint32_t`
/// (wraparound: `0 -> 0`, `1 -> 0xFFFFFFFF`). The port reproduces this as
/// `0u32.wrapping_sub(modis0(a))`, which is bit-for-bit the same:
/// `0u32 - 0 == 0`, `0u32 - 1 == 0xFFFFFFFF` under `u32` wraparound. The
/// alternative spelling `modis0(a).wrapping_neg()` is identical.
///
/// Operates on the redundant, non-canonical radix-2^51 form: the
/// reference's [`modis0`] first calls [`redc`] to canonicalise, so this
/// is the *value* zero predicate (i.e. `a == 0 mod p`), not raw-limb
/// equality to `[0, 0, 0, 0, 0]`. Any redundant representative of `0`
/// (e.g. `[1, 0, 0, 0, p4]`, which represents `1 + p4 << (51*4) = p`)
/// returns the all-ones mask.
pub fn fp_is_zero(a: &Fp) -> u32 {
    0u32.wrapping_sub(modis0(a))
}

/// Test whether the field elements represented by `a` and `b` are equal.
/// Returns `1` if they are, else `0` (a plain integer; [`fp_is_equal`]
/// negates it to the `0xFFFFFFFF`/`0` mask the rest of the codebase
/// consumes).
///
/// Mirrors the reference's
/// `static int modcmp(const spint *a, const spint *b)`:
///
/// ```c
/// spint c[5], d[5];
/// int i, eq = 1;
/// redc(a, c);
/// redc(b, d);
/// for (i = 0; i < 5; i++) {
///   eq &= (((c[i] ^ d[i]) - 1) >> 51) & 1;
/// }
/// return eq;
/// ```
///
/// Three faithfully reproduced subtleties:
/// 1. Both operands are reduced via [`redc`] (Montgomery to canonical) so
///    each limb of `c` and `d` is below `2^51`. The post-condition is
///    exactly what makes the `(x - 1) >> 51 & 1` zero-detect trick
///    correct in subtlety (2): the canonical-bounded property is the
///    same precondition [`modis0`] relies on for its single-OR-fold
///    variant.
/// 2. The per-limb application of the zero-detect trick. For each limb
///    `i`, let `x = c[i] ^ d[i]`: `x == 0` iff the canonical limbs
///    match. Then `(x - 1) >> 51 & 1` is `1` iff `x == 0` (`0 - 1`
///    wraps to `0xFFFF_FFFF_FFFF_FFFF` whose `>> 51` has bit 0 set;
///    canonical `x` is bounded above by the OR of canonical limbs each
///    below `2^51`, so `x < 2^51` and `x - 1 < 2^51` is non-negative,
///    `>> 51 == 0`, low bit `0`). This is precisely the bit pattern
///    [`modis0`] uses; the only difference is that [`modis0`] applies
///    it once to the OR-fold `d` of one operand's canonical limbs,
///    while `modcmp` applies it five times to per-limb XORs of two
///    operands' canonical limbs.
/// 3. `eq` is an `int` initialised to `1` and AND-ed with each per-limb
///    `{0, 1}` bit. The final value is `1` iff every limb XOR was zero,
///    else `0`. The port uses `u32` for `eq` (the inner-bit type, so
///    the AND is a `u32 & u32` matching the C `int & int` modulo the
///    `int01` invariant on the right operand: each per-limb bit is
///    `{0, 1}` regardless of the C type chosen). Equivalently `i32`
///    would work; `u32` makes the return-type widen at the call site
///    obvious. All arithmetic is unsigned `u64` wraparound (each
///    `(x - 1)` underflows silently when `x == 0`, exactly as the
///    reference's `spint == uint64_t`).
///
/// The return is a plain `u32` `{0, 1}` so the caller can either consume
/// the boolean directly (the way `modis1`'s wrapper ANDs two such bits)
/// or negate it to the `0xFFFFFFFF`/`0` mask the public boundary
/// returns.
fn modcmp(a: &Fp, b: &Fp) -> u32 {
    let mut c: Fp = [0u64; NWORDS_FIELD];
    let mut d: Fp = [0u64; NWORDS_FIELD];
    redc(a, &mut c);
    redc(b, &mut d);
    let mut eq: u32 = 1;
    for i in 0..NWORDS_FIELD {
        let x = c[i] ^ d[i];
        eq &= ((x.wrapping_sub(1) >> 51) & 1) as u32;
    }
    eq
}

/// GF(p) equality predicate: returns the constant-time mask `0xFFFFFFFF`
/// if `a` and `b` represent the same field element, else `0`.
///
/// Mirrors the reference's
/// `uint32_t fp_is_equal(const fp_t *a, const fp_t *b) { return
/// -(uint32_t)modcmp(*a, *b); }` from `fp_p5248_64.c:574..578`.
/// [`modcmp`] returns `0` or `1`; the `-(uint32_t)` cast turns that into
/// the `0xFFFFFFFF`/`0` mask downstream consumers ([`fp_select`],
/// [`fp_cswap`]'s LSB, the `fp2_is_equal` AND-of-two-fp_is_equals)
/// expect.
///
/// The cast chain `-(uint32_t)int01` in C is identical to the one
/// [`fp_is_zero`] uses: widen the `int` `{0, 1}` to `uint32_t` (no-op
/// for non-negative values), then negate as `uint32_t` (wraparound:
/// `0 -> 0`, `1 -> 0xFFFFFFFF`). The port reproduces this as
/// `0u32.wrapping_sub(modcmp(a, b))`, which is bit-for-bit the same:
/// `0u32 - 0 == 0`, `0u32 - 1 == 0xFFFFFFFF` under `u32` wraparound.
/// The alternative spelling `modcmp(a, b).wrapping_neg()` is identical.
///
/// Operates on the redundant, non-canonical radix-2^51 form: the
/// reference's [`modcmp`] first calls [`redc`] on **both** operands to
/// canonicalise them, so this is the *value* equality predicate
/// (i.e. `a == b mod p`), not raw-limb equality of the two limb
/// vectors. Two distinct redundant representatives of the same field
/// element (e.g. `[0, 0, 0, 0, 0]` and the radix-2^51 encoding of `p`
/// itself, `[MASK51, MASK51, MASK51, MASK51, P4 - 1]`, both
/// representing `0 mod p`) yield the all-ones mask.
///
/// Reflexivity (`fp_is_equal(a, a) == 0xFFFFFFFF`) holds for arbitrary
/// `a`: `redc(a) == redc(a)` bit-for-bit (the reference's [`redc`] is
/// deterministic with no global state), so every per-limb XOR is `0`,
/// every per-limb bit is `1`, the AND-fold is `1`, and the wrapper
/// returns the all-ones mask. Symmetry (`fp_is_equal(a, b) ==
/// fp_is_equal(b, a)`) holds bit-for-bit: XOR is symmetric, the per-limb
/// `{0, 1}` bits are independent of operand order, and the AND-fold is
/// commutative.
pub fn fp_is_equal(a: &Fp, b: &Fp) -> u32 {
    0u32.wrapping_sub(modcmp(a, b))
}

/// Montgomery n-residue conversion factor: the precomputed five-limb
/// constant `nres` multiplies a positional residue by to enter the
/// Montgomery domain. Transcribed verbatim from the reference's
/// `static void nres(const spint *m, spint *n)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:297..302`:
///
/// ```c
/// const spint c[5] = {0x4cccccccccf5cu, 0x1999999999999u,
///                     0x3333333333333u, 0x6666666666666u,
///                     0xcccccccccccu};
/// ```
///
/// Limb-by-limb correspondence (each constant copied byte-for-byte from
/// the C source's `0xNNNNNu` literal, with `u` stripped; the Rust
/// literal is zero-padded to a uniform 16-hex-digit grouping for
/// `clippy::unusual_byte_groupings` but the underlying value is the
/// same):
/// - limb 0: C `0x4cccccccccf5c` (13 hex digits, 51 bits), Rust
///   `0x0004_cccc_cccc_cf5c`.
/// - limb 1: C `0x1999999999999` (13 hex digits, 49 bits), Rust
///   `0x0001_9999_9999_9999`.
/// - limb 2: C `0x3333333333333` (13 hex digits, 50 bits), Rust
///   `0x0003_3333_3333_3333`.
/// - limb 3: C `0x6666666666666` (13 hex digits, 51 bits), Rust
///   `0x0006_6666_6666_6666`.
/// - limb 4: C `0xccccccccccc`   (11 hex digits, 44 bits), Rust
///   `0x0000_0ccc_cccc_cccc`. Note that this is **eleven** `c` digits,
///   not twelve; the C literal is `0xcccccccccccu` and dropping the
///   `0x` prefix and trailing `u` leaves an 11-character core.
///
/// Used by [`nres`]; the bit-for-bit correspondence with the reference's
/// constant is the load-bearing property. The independent value-level
/// oracle that the constant is the correct Montgomery `R^2 mod p` is
/// `nres(positional 1) == MONTGOMERY_ONE`, pinned in the unit test
/// `nres_of_positional_one_is_montgomery_one` below: that test exercises
/// the constant through the `modmul` core and confirms the result is the
/// already-cross-checked Montgomery representative of `1`.
const NRES_C: Fp = [
    0x0004_cccc_cccc_cf5c,
    0x0001_9999_9999_9999,
    0x0003_3333_3333_3333,
    0x0006_6666_6666_6666,
    0x0000_0ccc_cccc_cccc,
];

/// Convert a positional residue `m` to its Montgomery n-residue
/// representative `n = m * R^2 * R^-1 mod p == m * R mod p`.
///
/// Mirrors the reference's
/// `static void nres(const spint *m, spint *n)` at
/// `fp_p5248_64.c:297..302` exactly:
///
/// ```c
/// const spint c[5] = {0x4cccccccccf5cu, 0x1999999999999u,
///                     0x3333333333333u, 0x6666666666666u,
///                     0xcccccccccccu};
/// modmul(m, c, n);
/// ```
///
/// The single internal call is `modmul(m, c, n)` with the precomputed
/// constant [`NRES_C`] as the second operand. Operand order is preserved
/// for transcription clarity even though [`modmul`] is bit-exactly
/// commutative on the column sums (an established property pinned at the
/// `fp_mul` boundary); the reference puts `m` first so the port does too.
///
/// Note on aliasing: the reference passes `nres(a, a)` from inside
/// `modone` and `modint` (and elsewhere), reading and writing the same
/// buffer. The reference's [`modmul`] handles that aliasing correctly
/// because it writes `c[0]` only after the last column that needs `a[0]`
/// or `b[0]`, then `c[1]` only after the last column that needs `a[1]`
/// or `b[1]`, and so on. The Rust port's borrow checker will not permit
/// the equivalent `modmul(a, c, a)` call inside a Rust function,
/// however; callers that need `nres(a, a)` semantics must use a separate
/// destination buffer and copy back, exactly as [`modint`] does. The
/// `n` parameter is therefore documented as **distinct from `m`** at
/// the Rust call site, even though the underlying C is aliasing-safe.
fn nres(m: &Fp, n: &mut Fp) {
    modmul(m, &NRES_C, n);
}

/// Set `a` to the Montgomery n-residue representative of the integer
/// `x`. Mirrors the reference's
/// `static void modint(int x, spint *a)` at `fp_p5248_64.c:362..369`
/// exactly:
///
/// ```c
/// a[0] = (spint)x;
/// for (i = 1; i < 5; i++) a[i] = 0;
/// nres(a, a);
/// ```
///
/// Three faithfully reproduced subtleties:
///
/// 1. **Sign-extending cast at limb 0.** The reference writes
///    `a[0] = (spint)x` where `x` is `int` (`int32_t` on every
///    reasonable platform) and `spint` is `uint64_t`. The cast widens a
///    signed integer to an unsigned 64-bit one, which in C is
///    sign-extending for negative values: `(uint64_t)(int32_t)(-1)` is
///    `0xffff_ffff_ffff_ffff`. The port reproduces this as
///    `x as i64 as u64`: the `as i64` is the sign-extending widening
///    `i32 -> i64`, and the bit-cast `as u64` preserves the bit pattern.
///    Equivalently `x as u32 as u64 | sign_mask` would work but the
///    `x as i64 as u64` spelling makes the cast chain visible at the
///    source level.
/// 2. **Scratch destination for `nres(a, a)`.** The reference's `nres`
///    handles aliased in/out via the underlying [`modmul`]'s
///    column-ordering, but Rust's borrow checker does not permit a
///    direct `nres(a, a)`. The port computes into a fresh `tmp` then
///    moves `*a = tmp`; the resulting `a` is bit-equal to what the
///    reference would have produced in place (the [`modmul`] core
///    reads-then-writes column by column without revisiting an already-
///    written limb, so an in-place call would produce the same output).
/// 3. **Limbs 1..=4 are zeroed before the conversion.** The reference's
///    `modint` writes `a[0] = (spint)x` and then loops `a[i] = 0` for
///    `i in 1..5`, only then calls `nres`. The port matches: limb 0 is
///    set from `x` and limbs 1..=4 to `0` *before* the `nres` call, so
///    the input to `nres` is the positional residue with the same
///    pre-condition as the reference.
fn modint(x: i32, a: &mut Fp) {
    a[0] = x as i64 as u64;
    a[1] = 0;
    a[2] = 0;
    a[3] = 0;
    a[4] = 0;
    let mut tmp: Fp = [0u64; NWORDS_FIELD];
    nres(a, &mut tmp);
    *a = tmp;
}

/// GF(p) value setter: writes the Montgomery representative of the
/// integer `val`, narrowed first to `i32` and then sign-extended to
/// `u64` for the positional limb-0 write before `nres` converts the
/// result to Montgomery form.
///
/// Mirrors the reference's
/// `void fp_set_small(fp_t *x, const digit_t val)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:550..554`:
///
/// ```c
/// void fp_set_small(fp_t *x, const digit_t val) {
///   modint((int)val, *x);
/// }
/// ```
///
/// The full `digit_t` (`uint64_t` under `RADIX_64`) argument is taken
/// but immediately narrowed to `int` (`int32_t` on every reasonable
/// platform). This is the load-bearing observation: only the low 32
/// bits of `val` are observable through the boundary. The C cast chain
/// `(int)val` (where val is `uint64_t`) drops the high 32 bits then
/// reinterprets the low 32 bits as `int32_t`, so for `val ==
/// 0x80000000` the narrowed `int` is `INT32_MIN` (`-2_147_483_648`)
/// and the subsequent `a[0] = (spint)x` sign-extends to
/// `0xffff_ffff_8000_0000`. The port reproduces this as `val as i32`,
/// which is the bit-preserving `u64 -> i32` narrowing (Rust's `as` on
/// integer types performs the same truncation-and-reinterpret the C
/// cast chain does).
///
/// The intended argument domain per the function name is "small positive
/// integer fitting in int" (so the high 32 bits are zero and the low 32
/// bits are below `2^31`), but the C signature accepts any `digit_t` and
/// the port must reproduce its behaviour on the full battery, including
/// the high-bits-ignored narrowing and the sign-extension of values
/// above `2^31 - 1`. The differential battery exercises both: 12
/// edge prefills cross 11 representative `val` values including
/// `0x80000000`, `0xffffffff`, `0x100000000`, `0xffffffff00000000`,
/// `0xffffffffffffffff` (1132 records total, 1048 of which exercise the
/// high-bits-ignored narrowing).
///
/// Internal helpers landed here for the first time: [`NRES_C`] (the
/// Montgomery n-residue conversion factor), [`nres`] (the positional ->
/// Montgomery conversion threaded through the already-ported `modmul`
/// core), and [`modint`] (the int-to-Montgomery setter, the per-limb
/// shape `modone` and `modint` share). Going forward, anything in the
/// reference that calls `nres` or `modint` (or both, through `modmul`'s
/// `b * modint(_)` style fallthrough) can reuse these directly.
pub fn fp_set_small(out: &mut Fp, val: u64) {
    modint(val as i32, out);
}

/// Modular multiplication by an integer: `c = a * b mod 2p`, where `b` is
/// an `int` first widened through [`modint`] to its Montgomery
/// representative and then folded into the running multiplication via
/// [`modmul`].
///
/// Mirrors the reference's
/// `inline static void modmli(const spint *a, int b, spint *c)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:372..376`:
///
/// ```c
/// inline static void modmli(const spint *a, int b, spint *c) {
///   spint t[5];
///   modint(b, t);
///   modmul(a, t, c);
/// }
/// ```
///
/// It is the thin combinator the reference uses whenever a field element
/// must be multiplied by a small integer constant: build the Montgomery
/// representative of that integer with [`modint`] (so the cross product
/// stays in the Montgomery domain), then run a single [`modmul`]. Because
/// `modint`'s positional limb-0 write is sign-extending from `int32_t`,
/// the same high-bits-ignored narrowing [`fp_set_small`] documents applies
/// here too, only the low 32 bits of the caller's argument are observable
/// once the public wrapper [`fp_mul_small`] has performed its
/// `uint32_t -> int` cast.
///
/// No new arithmetic is introduced at this level: the port is exactly the
/// reference's two-line body, transcribed identifier-for-identifier with a
/// local five-limb scratch `t` for the `modint` output, then a single
/// `modmul(a, &t, c)` call into the already-ported core. The `t` scratch
/// is required for the same borrow-checker reason [`modint`] needs its
/// `tmp` buffer in the underlying `nres(a, a)` call: Rust will not permit
/// the equivalent aliased borrow, and a fresh stack-resident `Fp` is the
/// natural mirror of the reference's `spint t[5];`.
fn modmli(a: &Fp, b: i32, c: &mut Fp) {
    let mut t: Fp = [0u64; NWORDS_FIELD];
    modint(b, &mut t);
    modmul(a, &t, c);
}

/// GF(p) Montgomery multiplication by a small integer: `out = a * val * R^-1
/// mod p`, in the redundant radix-2^51 representation, reduced to less than
/// `2p`.
///
/// Mirrors the reference's `void fp_mul_small(fp_t *x, const fp_t *a,
/// const uint32_t val)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:557..560`:
///
/// ```c
/// void fp_mul_small(fp_t *x, const fp_t *a, const uint32_t val) {
///   modmli(*a, (int)val, *x);
/// }
/// ```
///
/// The full `uint32_t` argument is taken but immediately narrowed to
/// `int` (`int32_t` on every reasonable platform). The narrowing is the
/// load-bearing observation: the boundary's effective scalar domain is
/// `int32_t`, not `uint32_t`, so the upper half of values above
/// `0x7fffffff` is reinterpreted as negative and the subsequent
/// `(spint)x` cast inside [`modint`] sign-extends to `0xffffffff_xxxxxxxx`.
/// For example, `val == 0x80000000` becomes the `int32_t` `INT32_MIN`
/// (`-2_147_483_648`), and `modint` writes `0xffff_ffff_8000_0000` at
/// positional limb 0 before [`nres`] converts the result to its Montgomery
/// representative.
///
/// The port reproduces this as `val as i32`: Rust's `as` on integer types
/// performs the bit-preserving `u32 -> i32` reinterpret the C `(int)val`
/// cast does, and the subsequent `i32 -> u64` sign-extension is hidden
/// inside [`modint`]'s already-ported `x as i64 as u64` chain. No new
/// arithmetic is introduced at this boundary: the port is exactly the
/// reference's one-line body, dispatching through [`modmli`] (added here
/// alongside the existing [`modint`], [`nres`], [`redc`], [`modfsb`],
/// [`flatten`], [`modis0`], [`modcmp`] internal helpers) and reusing
/// [`modmul`] for the cross product.
///
/// Two faithfully reproduced subtleties from the underlying [`modmul`]
/// path:
/// 1. The output is *not* fully canonical: limbs 0..=3 are below `2^51`
///    (the column mask) but limb 4 is left unmasked (the final
///    `c[4] = (spint)t` in `modmul` is a full 64-bit truncation of the
///    residual accumulator, no `& mask`), exactly as the reference
///    leaves it. Compare field elements with the value-level
///    [`fp_is_equal`], never by raw-limb equality. The lone exception is
///    `val == 0`: the cross product is bit-exactly the canonical
///    all-zero representative for arbitrary `a` (every column sum is a
///    multiple of the all-zero Montgomery image of `0`).
/// 2. The Montgomery domain is preserved: if `a == A * R mod p`
///    positionally and `val` narrows to the signed integer `v`, then
///    `modmli(a, v)` computes `(A * R) * (v * R) * R^-1 == (A * v) * R
///    mod p`. The output is the Montgomery representative of `A * v`,
///    ready to be consumed by further Montgomery-domain operations
///    without re-conversion. This is exactly the same Montgomery-in,
///    Montgomery-out pattern [`fp_half`] and [`fp_div3`] exhibit, only
///    with the integer operand built dynamically through [`modint`]
///    rather than supplied as a precomputed constant.
pub fn fp_mul_small(out: &mut Fp, a: &Fp, val: u32) {
    modmli(a, val as i32, out);
}

/// In-place n-fold Montgomery squaring `a <- a^(2^n) mod 2p` in the
/// redundant radix-2^51 representation.
///
/// Mirrors the reference's
/// `static void modnsqr(spint *a, int n)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:228..233` exactly:
///
/// ```c
/// static void modnsqr(spint *a, int n) {
///   int i;
///   for (i = 0; i < n; i++) {
///     modsqr(a, a);
///   }
/// }
/// ```
///
/// A trivial loop calling [`modsqr`] in place `n` times. The reference
/// passes the same `a` as both source and destination of `modsqr`; the
/// port mirrors the in-place pattern by squaring into a scratch and
/// copying back limb for limb (Rust's borrow checker would reject the
/// aliased `modsqr(a, a)` directly even though the underlying [`modmul`]
/// is aliasing-safe by virtue of its column-ordered reads-then-writes).
/// The intermediate `tmp` is bit-equal to what the reference would
/// produce in place.
///
/// `n` is `int` in the reference; this port takes `i32` and treats a
/// zero or negative count as a no-op (the reference's `for` loop has the
/// same behaviour for `n <= 0`). The callers in [`modpro`] and
/// [`modinv`] pass small non-negative literals.
fn modnsqr(a: &mut Fp, n: i32) {
    for _ in 0..n {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modsqr(a, &mut tmp);
        *a = tmp;
    }
}

/// Compute the progenitor of `w`, `z = w^((p-3)/4) mod p`, the Montgomery
/// representative used by [`modsqrt`] and [`modinv`] to extract roots and
/// inverses via Fermat's little theorem.
///
/// Mirrors the reference's
/// `static void modpro(const spint *w, spint *z)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:236..281`
/// statement-for-statement. The reference is a hand-built fixed addition
/// chain of [`modsqr`]/[`modmul`]/[`modnsqr`] calls (see the upstream
/// comment "// Calculate progenitor"); the port preserves the exact
/// sequence and scratch-variable usage so the per-step intermediate
/// Montgomery representatives are bit-identical to the reference's.
///
/// Six five-limb scratches (`x`, `t0`, `t1`, `t2`, `t3`, `t4`) mirror the
/// reference's `spint x[5]; ... t4[5];` declarations. The reference uses
/// in-place [`modsqr`] on the same buffer (e.g. `modsqr(t0, z)`) and
/// in-place [`modmul`] with one operand aliased to the destination (e.g.
/// `modmul(x, z, z)`); the port uses fresh stack-resident scratches for
/// the [`modmul`] aliased case (the Rust borrow checker rejects the
/// equivalent aliased borrow even though the underlying column-ordered
/// writes are aliasing-safe), and likewise routes [`modsqr`] through a
/// `tmp` scratch followed by `*dst = tmp;`. The output `z` is bit-equal
/// to what the reference would produce.
///
/// The progenitor is the building block of:
/// - [`modinv`] when no precomputed `h` is supplied (computes `progenitor`
///   then squares twice and multiplies by `x` to obtain `x^-1`).
/// - [`modqr`] (squares the progenitor and multiplies by `x`, returning
///   `1` for a quadratic residue via [`modis1`]).
/// - [`modsqrt`] (multiplies the progenitor by `x` to obtain `sqrt(x)`).
fn modpro(w: &Fp, z: &mut Fp) {
    let x: Fp = *w;
    let mut t0: Fp = [0u64; NWORDS_FIELD];
    let mut t1: Fp = [0u64; NWORDS_FIELD];
    let mut t2: Fp = [0u64; NWORDS_FIELD];
    let mut t3: Fp = [0u64; NWORDS_FIELD];
    let mut t4: Fp;
    modsqr(&x, z);
    modmul(&x, z, &mut t0);
    {
        let src = t0;
        modsqr(&src, z);
    }
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&x, z, &mut tmp);
        *z = tmp;
    }
    modsqr(z, &mut t1);
    modsqr(&t1, &mut t3);
    modsqr(&t3, &mut t2);
    t4 = t2;
    modnsqr(&mut t4, 3);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t2, &t4, &mut tmp);
        t2 = tmp;
    }
    t4 = t2;
    modnsqr(&mut t4, 6);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t2, &t4, &mut tmp);
        t2 = tmp;
    }
    t4 = t2;
    modnsqr(&mut t4, 2);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t3, &t4, &mut tmp);
        t3 = tmp;
    }
    modnsqr(&mut t3, 13);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t2, &t3, &mut tmp);
        t2 = tmp;
    }
    t3 = t2;
    modnsqr(&mut t3, 27);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t2, &t3, &mut tmp);
        t2 = tmp;
    }
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(z, &t2, &mut tmp);
        *z = tmp;
    }
    t2 = *z;
    modnsqr(&mut t2, 4);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t1, &t2, &mut tmp);
        t1 = tmp;
    }
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t0, &t1, &mut tmp);
        t0 = tmp;
    }
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t1, &t0, &mut tmp);
        t1 = tmp;
    }
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t0, &t1, &mut tmp);
        t0 = tmp;
    }
    modmul(&t1, &t0, &mut t2);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t0, &t2, &mut tmp);
        t0 = tmp;
    }
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t1, &t0, &mut tmp);
        t1 = tmp;
    }
    modnsqr(&mut t1, 63);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t0, &t1, &mut tmp);
        t1 = tmp;
    }
    modnsqr(&mut t1, 64);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(&t0, &t1, &mut tmp);
        t0 = tmp;
    }
    modnsqr(&mut t0, 57);
    {
        let mut tmp: Fp = [0u64; NWORDS_FIELD];
        modmul(z, &t0, &mut tmp);
        *z = tmp;
    }
}

/// Compute the modular inverse `z = x^-1 mod p` via Fermat (`x * x^(p-2)
/// == 1 mod p`). Mirrors the reference's
/// `static void modinv(const spint *x, const spint *h, spint *z)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:284..295` exactly:
///
/// ```c
/// static void modinv(const spint *x, const spint *h, spint *z) {
///   spint s[5];
///   spint t[5];
///   if (h == NULL) {
///     modpro(x, t);
///   } else {
///     modcpy(h, t);
///   }
///   modcpy(x, s);
///   modnsqr(t, 2);
///   modmul(s, t, z);
/// }
/// ```
///
/// `h` is the optional precomputed progenitor of `x`; when `None`, this
/// port computes it via [`modpro`] (the reference does the same on the
/// `h == NULL` branch). The port's `Option<&Fp>` mirrors the reference's
/// nullable-pointer contract: `Some(h)` is the precomputed-progenitor
/// case, `None` triggers the [`modpro`] call. The only [`fp_inv`] caller
/// in this port passes `None`, so the precomputed-progenitor branch is
/// exercised by the differential battery only via the cross-validating
/// property tests that exist for symmetry with the reference.
///
/// The inversion identity is `x * x^((p-3)/4)^4 == x * x^(p-3) ==
/// x^(p-2) == x^-1 mod p` (using `progenitor == x^((p-3)/4)` and the
/// final two squarings raising to `^4`); the structure factors into a
/// progenitor evaluation, two further squarings, and one final multiply,
/// exactly as the reference encodes it.
fn modinv(x: &Fp, h: Option<&Fp>, z: &mut Fp) {
    let mut t: Fp = [0u64; NWORDS_FIELD];
    match h {
        None => modpro(x, &mut t),
        Some(h) => t = *h,
    }
    let s: Fp = *x;
    modnsqr(&mut t, 2);
    modmul(&s, &t, z);
}

/// Test whether `a` represents the field's multiplicative identity
/// (Montgomery `1`). Returns `1` if so, else `0` (a plain integer; the
/// `-(uint32_t)` wrapper turns it into the `0xFFFFFFFF`/`0` mask when
/// downstream callers need it; here the only consumer is [`modqr`], which
/// ORs the bit with [`modis0`]'s bit).
///
/// Mirrors the reference's
/// `static int modis1(const spint *a)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:317..329` exactly:
///
/// ```c
/// static int modis1(const spint *a) {
///   int i;
///   spint c[5];
///   spint c0;
///   spint d = 0;
///   redc(a, c);
///   for (i = 1; i < 5; i++) {
///     d |= c[i];
///   }
///   c0 = (spint)c[0];
///   return ((spint)1 & ((d - (spint)1) >> 51u) &
///           (((c0 ^ (spint)1) - (spint)1) >> 51u));
/// }
/// ```
///
/// Three faithfully reproduced subtleties:
/// 1. `redc(a, c)` canonicalises any redundant representative to limbs
///    below `2^51` with positional value below `p`. The canonical
///    representative of Montgomery `1` is the positional `[1, 0, 0, 0,
///    0]`, *not* the Montgomery [`MONTGOMERY_ONE`] (that's the
///    Montgomery domain representative pre-`redc`); the OR-fold of `c[1
///    ..= 4]` is `0` iff the canonical limbs 1..=4 are zero, the
///    canonical-bounded property that makes the `(d - 1) >> 51 & 1`
///    zero-detect trick correct.
/// 2. The same zero-detect trick is then applied to `c[0] ^ 1`: this
///    bit is `1` iff `c[0] == 1`. The combined AND is `1` iff every
///    canonical limb 1..=4 is `0` AND `c[0] == 1`, i.e. iff `a == 1
///    mod p`.
/// 3. All arithmetic is unsigned `u64` wraparound (each `x.wrapping_sub(1)`
///    underflows silently when `x == 0`, exactly as the reference's
///    `spint == uint64_t`).
fn modis1(a: &Fp) -> u32 {
    let mut c: Fp = [0u64; NWORDS_FIELD];
    redc(a, &mut c);
    let mut d: u64 = 0;
    for limb in c.iter().skip(1) {
        d |= *limb;
    }
    let c0 = c[0];
    let bit_d = (d.wrapping_sub(1) >> 51) & 1;
    let bit_c0 = ((c0 ^ 1).wrapping_sub(1) >> 51) & 1;
    (bit_d & bit_c0) as u32
}

/// Test whether `x` is a quadratic residue mod `p`. Returns `1` if so or
/// if `x` represents the field zero (per the reference's contract: zero
/// is conventionally treated as a square), else `0` (a plain integer;
/// the public [`fp_is_square`] wrapper negates it to the
/// `0xFFFFFFFF`/`0` mask).
///
/// Mirrors the reference's
/// `static int modqr(const spint *h, const spint *x)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:379..389` exactly:
///
/// ```c
/// static int modqr(const spint *h, const spint *x) {
///   spint r[5];
///   if (h == NULL) {
///     modpro(x, r);
///     modsqr(r, r);
///   } else {
///     modsqr(h, r);
///   }
///   modmul(r, x, r);
///   return modis1(r) | modis0(x);
/// }
/// ```
///
/// The Euler criterion in Montgomery form: `r = x^((p-1)/2) ==
/// x^((p-3)/4 * 2) * x^1 == progenitor^2 * x` is `1` iff `x` is a non-zero
/// square, `-1` iff `x` is a non-square, and `0` iff `x == 0`. The
/// `modis1(r) | modis0(x)` OR returns `1` for both the square (`r == 1`)
/// and the zero (`x == 0`) branches.
///
/// `h` is the optional precomputed progenitor of `x` (same role as in
/// [`modinv`]); when `None`, the reference computes it via [`modpro`]
/// and then squares it once, equivalent to `(x^((p-3)/4))^2 ==
/// x^((p-3)/2)`. When `Some(h)`, the caller has already supplied
/// `progenitor` and only one squaring is needed.
///
/// The [`modmul`] aliasing `modmul(r, x, r)` (output aliased with the
/// first input) and `modsqr`'s in-place form are routed through a `tmp`
/// scratch in the port for the same borrow-checker reason [`modpro`]
/// scratches its aliased modmul/modsqr calls.
fn modqr(h: Option<&Fp>, x: &Fp) -> u32 {
    let mut r: Fp = [0u64; NWORDS_FIELD];
    match h {
        None => {
            modpro(x, &mut r);
            let src = r;
            modsqr(&src, &mut r);
        }
        Some(h) => {
            modsqr(h, &mut r);
        }
    }
    {
        let src = r;
        modmul(&src, x, &mut r);
    }
    modis1(&r) | modis0(x)
}

/// Compute the square root `r = sqrt(x) mod p`, when one exists. Caller
/// must use [`modqr`] to determine whether `x` is in fact a quadratic
/// residue; on a non-residue, the returned value is meaningless garbage
/// (the reference makes no defensive check, the port follows).
///
/// Mirrors the reference's
/// `static void modsqrt(const spint *x, const spint *h, spint *r)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:427..437` exactly:
///
/// ```c
/// static void modsqrt(const spint *x, const spint *h, spint *r) {
///   spint s[5];
///   spint y[5];
///   if (h == NULL) {
///     modpro(x, y);
///   } else {
///     modcpy(h, y);
///   }
///   modmul(y, x, s);
///   modcpy(s, r);
/// }
/// ```
///
/// `h` is the optional precomputed progenitor; when `None`, the port
/// (mirroring the reference) computes it via [`modpro`]. The identity
/// underpinning the construction is `progenitor(x) * x ==
/// x^((p-3)/4) * x == x^((p+1)/4) == sqrt(x) mod p` (for `p == 3 mod 4`,
/// the case that holds for `p5248`); a single [`modmul`] of the
/// progenitor with `x` yields the root.
///
/// The `s`/`r` indirection in the reference (compute into `s`, then
/// `modcpy(s, r);`) is preserved at the source level even though it
/// reduces to a single buffer copy at the boundary; this is structurally
/// the same pattern [`modqr`] uses for its `r` write.
fn modsqrt(x: &Fp, h: Option<&Fp>, r: &mut Fp) {
    let mut s: Fp = [0u64; NWORDS_FIELD];
    let mut y: Fp = [0u64; NWORDS_FIELD];
    match h {
        None => modpro(x, &mut y),
        Some(h) => y = *h,
    }
    modmul(&y, x, &mut s);
    *r = s;
}

/// In-place left shift by `n < 51` bits across the five limbs of the
/// redundant radix-2^51 layout. Mirrors the reference's
/// `static void modshl(unsigned int n, spint *a)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:440..447` exactly:
///
/// ```c
/// static void modshl(unsigned int n, spint *a) {
///   int i;
///   a[4] = ((a[4] << n)) | (a[3] >> (51u - n));
///   for (i = 3; i > 0; i--) {
///     a[i] = ((a[i] << n) & (spint)0x7ffffffffffff) | (a[i - 1] >> (51u - n));
///   }
///   a[0] = (a[0] << n) & (spint)0x7ffffffffffff;
/// }
/// ```
///
/// Two faithfully reproduced subtleties:
/// 1. Limb 4 is shifted *without* the per-limb mask (so any bits shifted
///    above bit 50 stay in limb 4 unmasked); limbs 0..=3 are masked back
///    down to 51 bits each, the same `MASK51` the column writes use.
/// 2. The descending loop `for (i = 3; i > 0; i--)` ensures each `a[i]`
///    is read for its low bits *after* its high bits have already been
///    consumed by the previous iteration's write to `a[i+1]`; the port
///    preserves the iteration order so the bit-for-bit correspondence
///    is visible at the source level.
///
/// `n` is `unsigned int` in the reference (so `n == 0` produces an
/// undefined-behaviour `51 - 0 == 51` shift on `u64`, which is in-range
/// at exactly the type-width boundary; the reference's only caller
/// [`modimp`] passes `8`, so `n == 0` is out of the contract and the
/// port matches by passing `n` through to the same shifts without
/// guard). The port takes `u32` for the same reason [`modshr`] does.
fn modshl(n: u32, a: &mut Fp) {
    a[4] = (a[4] << n) | (a[3] >> (51 - n));
    for i in (1..=3).rev() {
        a[i] = ((a[i] << n) & MASK51) | (a[i - 1] >> (51 - n));
    }
    a[0] = (a[0] << n) & MASK51;
}

/// In-place right shift by `n < 51` bits across the five limbs of the
/// redundant radix-2^51 layout, returning the shifted-out low bits as a
/// plain integer (the `a[0] & ((1 << n) - 1)` mask captured *before* the
/// shift).
///
/// Mirrors the reference's
/// `static int modshr(unsigned int n, spint *a)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:450..458` exactly:
///
/// ```c
/// static int modshr(unsigned int n, spint *a) {
///   int i;
///   spint r = a[0] & (((spint)1 << n) - (spint)1);
///   for (i = 0; i < 4; i++) {
///     a[i] = (a[i] >> n) | ((a[i + 1] << (51u - n)) & (spint)0x7ffffffffffff);
///   }
///   a[4] = a[4] >> n;
///   return r;
/// }
/// ```
///
/// `n` is `unsigned int` in the reference; the port takes `u32`. The
/// return is `int` in the reference; the port returns `u64` so the
/// caller in [`fp_encode`] (which immediately downcasts to `u8`) can use
/// the value directly without losing the low byte. The reference's
/// `(int)r` truncation in the caller pattern is bit-equivalent to the
/// port's `r as u8`.
///
/// Limb 4 is shifted without the per-limb mask (consistent with
/// [`modshl`] and the limb-4-is-unmasked invariant `modmul`/`modsqr`
/// leave behind); limbs 0..=3 OR in the high-bit-borrow from `a[i+1]`
/// masked to 51 bits.
fn modshr(n: u32, a: &mut Fp) -> u64 {
    let r = a[0] & ((1u64 << n) - 1);
    for i in 0..4 {
        a[i] = (a[i] >> n) | ((a[i + 1] << (51 - n)) & MASK51);
    }
    a[4] >>= n;
    r
}

/// Montgomery representative of `2^256 mod p` on the level-1 generic
/// field, used by [`fp_decode_reduce`] to fold 32-byte blocks into a
/// running accumulator (each preceding block is multiplied by `R2` to
/// shift it `256` positional bits up before the next block is added).
///
/// Transcribed verbatim from the reference's
/// `static const digit_t R2[NWORDS_FIELD]` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:544..548`:
///
/// ```c
/// // Montgomery representation of 2^256
/// static const digit_t R2[NWORDS_FIELD] = { 0x0001999999999eb8,
///                                           0x0003333333333333,
///                                           0x0006666666666666,
///                                           0x0004cccccccccccc,
///                                           0x0000199999999999 };
/// ```
///
/// The constant is taken verbatim; no derivation is performed at the
/// port. The bit-for-bit correspondence with the reference is the
/// load-bearing property; the value-level identity that `R2 == 2^256 *
/// R mod p` (where `R == 2^256` for the level-1 `RADIX_64` build, so
/// `R2 == 2^512 mod p`, the Montgomery `R^2`) is what makes the chain
/// `d = d * R2 + decode(block)` in [`fp_decode_reduce`] equivalent to
/// the positional `d = d * 2^256 + value(block)` reduction.
const R2: Fp = [
    0x0001_9999_9999_9eb8,
    0x0003_3333_3333_3333,
    0x0006_6666_6666_6666,
    0x0004_cccc_cccc_cccc,
    0x0000_1999_9999_9999,
];

/// GF(p) modular exponentiation by `(p-3)/4` (the "progenitor" of `a`).
///
/// Mirrors the reference's `void fp_exp3div4(fp_t *out, const fp_t *a)`
/// at `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:652..656`:
///
/// ```c
/// void fp_exp3div4(fp_t *out, const fp_t *a) {
///   modpro(*a, *out);
/// }
/// ```
///
/// A thin wrapper around the internal [`modpro`] (the hand-built fixed
/// addition chain). Used by callers that need the progenitor explicitly
/// to feed [`fp_inv`] or [`fp_sqrt`] with a precomputed `h` argument
/// (the `Some(h)` branch of [`modinv`]/[`modsqrt`]) and amortise the
/// progenitor cost.
///
/// As with the rest of the Montgomery-domain gf ports, the output is in
/// the redundant `[0, 2p)` representation (limbs 0..=3 below `2^51`, limb
/// 4 unmasked). The differential boundary records the raw five-limb
/// output bit-for-bit against the reference's.
pub fn fp_exp3div4(out: &mut Fp, a: &Fp) {
    modpro(a, out);
}

/// GF(p) modular inverse, in place: `x <- x^-1 mod p`. On the field
/// zero, the reference returns whatever Fermat would (since the chain
/// computes `x^(p-2)`, an all-zero input squares-and-multiplies down to
/// the canonical zero); the port mirrors that behaviour.
///
/// Mirrors the reference's `void fp_inv(fp_t *x)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:628..632`:
///
/// ```c
/// void fp_inv(fp_t *x) {
///   modinv(*x, NULL, *x);
/// }
/// ```
///
/// The reference passes the same buffer as both input `x` and destination
/// `z` of [`modinv`]; the underlying [`modmul`] tolerates this aliasing,
/// but the Rust port resolves the borrow-checker conflict by snapshotting
/// `x` into a local `input` and computing through that, then writing back
/// to `x` from the [`modinv`] destination. The output is bit-equal to
/// what the reference would have produced in place.
///
/// `h` is hard-coded to `None` at this boundary (matching the
/// reference's `modinv(*x, NULL, *x)`), so [`modinv`] computes the
/// progenitor internally via [`modpro`]. The `Some(h)` branch of
/// [`modinv`] is exercised only by future callers that have a
/// precomputed progenitor.
pub fn fp_inv(x: &mut Fp) {
    let input = *x;
    modinv(&input, None, x);
}

/// GF(p) quadratic-residue predicate: returns the constant-time mask
/// `0xFFFFFFFF` if `a` is a quadratic residue mod `p` (or the field
/// zero), else `0`.
///
/// Mirrors the reference's `uint32_t fp_is_square(const fp_t *a)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:634..638`:
///
/// ```c
/// uint32_t fp_is_square(const fp_t *a) {
///   return -(uint32_t)modqr(NULL, *a);
/// }
/// ```
///
/// The `-(uint32_t)` cast turns [`modqr`]'s `{0, 1}` return into the
/// `0xFFFFFFFF`/`0` mask the rest of the codebase consumes (the same
/// negation [`fp_is_zero`] and [`fp_is_equal`] use). `h` is hard-coded
/// to `None` (the reference's `modqr(NULL, *a)`), so [`modqr`] computes
/// the progenitor internally.
///
/// Per the reference, the field zero is treated as a square (the
/// `modis0(x)` branch of [`modqr`]'s final OR returns `1` on zero); the
/// port matches.
pub fn fp_is_square(a: &Fp) -> u32 {
    0u32.wrapping_sub(modqr(None, a))
}

/// GF(p) modular square root, in place: `a <- sqrt(a) mod p` for the
/// caller's chosen branch. On a non-residue input, the returned value is
/// meaningless garbage (the reference makes no defensive check; the port
/// follows). Caller should first test with [`fp_is_square`].
///
/// Mirrors the reference's `void fp_sqrt(fp_t *a)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:640..644`:
///
/// ```c
/// void fp_sqrt(fp_t *a) {
///   modsqrt(*a, NULL, *a);
/// }
/// ```
///
/// The reference passes the same buffer as both input `x` and destination
/// `r` of [`modsqrt`]; the underlying [`modmul`] tolerates this aliasing,
/// but the Rust port resolves the borrow-checker conflict by
/// snapshotting `a` into a local `input` and writing the result back into
/// `a`. The output is bit-equal to what the reference would have produced
/// in place.
///
/// `h` is hard-coded to `None` (the reference's `modsqrt(*a, NULL,
/// *a)`), so [`modsqrt`] computes the progenitor internally via
/// [`modpro`].
pub fn fp_sqrt(a: &mut Fp) {
    let input = *a;
    modsqrt(&input, None, a);
}

/// GF(p) canonical serialization to 32 little-endian bytes.
///
/// Mirrors the reference's `void fp_encode(void *dst, const fp_t *a)`
/// at `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:664..675`:
///
/// ```c
/// void fp_encode(void *dst, const fp_t *a) {
///   int i;
///   spint c[5];
///   redc(*a, c);
///   for (i = 0; i < 32; i++) {
///     ((char *)dst)[i] = c[0] & (spint)0xff;
///     (void)modshr(8, c);
///   }
/// }
/// ```
///
/// The reference's `modexp`-derived path: [`redc`] the Montgomery
/// representative to its canonical positional form, then iteratively
/// peel off the low byte of limb 0 and right-shift the five-limb value
/// by 8 bits (the [`modshr`] call), 32 times in total. The output is the
/// 32-byte little-endian encoding of the canonical residue mod `p`.
///
/// The destination buffer must be exactly 32 bytes; the port's `&mut
/// [u8; 32]` is the Rust analogue of the reference's `void *dst` (the
/// reference does no bounds check, the port takes a fixed-size array
/// reference so the contract is encoded in the type).
pub fn fp_encode(dst: &mut [u8; 32], a: &Fp) {
    let mut c: Fp = [0u64; NWORDS_FIELD];
    redc(a, &mut c);
    for byte in dst.iter_mut() {
        *byte = (c[0] & 0xff) as u8;
        let _ = modshr(8, &mut c);
    }
}

/// GF(p) canonical deserialization from 32 little-endian bytes. Returns
/// the constant-time mask `0xFFFFFFFF` if the decoded value was in the
/// canonical range `[0, p)`, else `0` (in the non-canonical case, the
/// output `d` is zeroed). Mirrors the reference's `uint32_t
/// fp_decode(fp_t *d, const void *src)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:677..698`:
///
/// ```c
/// uint32_t fp_decode(fp_t *d, const void *src) {
///   int i;
///   spint res;
///   const unsigned char *b = src;
///   for (i = 0; i < 5; i++) {
///     (*d)[i] = 0;
///   }
///   for (i = 31; i >= 0; i--) {
///     modshl(8, *d);
///     (*d)[0] += (spint)b[i];
///   }
///   res = (spint)-modfsb(*d);
///   nres(*d, *d);
///   // If the value was canonical then res = -1; otherwise, res = 0
///   for (i = 0; i < 5; i++) {
///     (*d)[i] &= res;
///   }
///   return (uint32_t)res;
/// }
/// ```
///
/// Three faithfully reproduced subtleties:
/// 1. The byte-by-byte input is folded in descending address order
///    (`for (i = 31; i >= 0; i--)`): each iteration shifts the running
///    `d` left by 8 bits via [`modshl`] then adds the byte at the next
///    descending input position into `d[0]`. The cumulative effect is
///    `d = sum(b[i] << (8 * i)) for i in 0..32`, the little-endian
///    positional decoding.
/// 2. [`modfsb`] returns `1` if the decoded value is below `p` (the
///    trial subtraction was undone) and `0` otherwise. The reference
///    takes that as a `spint == uint64_t`, negates it
///    (`(spint)-modfsb(*d)`), so the result is `0xffff_ffff_ffff_ffff`
///    on the canonical-in-range branch and `0` on the out-of-range
///    branch. The port reproduces this via `0u64.wrapping_sub(...)`.
/// 3. [`nres`] is called *after* the canonical-range check but *before*
///    the per-limb `& res` mask: the limbs of `d` therefore hold the
///    Montgomery representative of the canonical value when `res ==
///    -1`, and arbitrary [`nres`] output (irrelevant) when `res == 0`,
///    which the subsequent `& res` zeroes out. The returned `uint32_t`
///    is the low 32 bits of the `spint` `res`: `0xffff_ffff` on canonical,
///    `0` on out-of-range.
pub fn fp_decode(d: &mut Fp, src: &[u8; 32]) -> u32 {
    for limb in d.iter_mut() {
        *limb = 0;
    }
    for i in (0..32).rev() {
        modshl(8, d);
        d[0] = d[0].wrapping_add(src[i] as u64);
    }
    let res = 0u64.wrapping_sub(modfsb(d) as u64);
    let input = *d;
    nres(&input, d);
    for limb in d.iter_mut() {
        *limb &= res;
    }
    res as u32
}

/// 8-bit add-with-carry on two `u64`s, producing one `u64` sum and the
/// outgoing carry bit. Mirrors the reference's
/// `static inline unsigned char add_carry(unsigned char cc, spint a,
/// spint b, spint *d)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:700..706` exactly:
///
/// ```c
/// static inline unsigned char
/// add_carry(unsigned char cc, spint a, spint b, spint *d) {
///   udpint t = (udpint)a + (udpint)b + cc;
///   *d = (spint)t;
///   return (unsigned char)(t >> Wordlength);
/// }
/// ```
///
/// `udpint == __uint128_t`, so the sum fits without overflow; the low 64
/// bits go to `d` and the high bit (one of `Wordlength == 64`) is the
/// outgoing carry. Used by [`partial_reduce`] to fold the high 8 bits of
/// a 256-bit value into the low 248 bits without losing precision.
fn add_carry(cc: u8, a: u64, b: u64, d: &mut u64) -> u8 {
    let t: u128 = (a as u128).wrapping_add(b as u128).wrapping_add(cc as u128);
    *d = t as u64;
    (t >> 64) as u8
}

/// Reduce a 256-bit value to a 248-bit value congruent mod `p == 5 *
/// 2^248 - 1`, in place on a four-`u64` array (NOT the Montgomery 5-limb
/// form). Mirrors the reference's
/// `static void partial_reduce(spint *out, const spint *src)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:708..726` exactly:
///
/// ```c
/// static void partial_reduce(spint *out, const spint *src) {
///   spint h, l, quo, rem;
///   unsigned char cc;
///   // Split value in high (8 bits) and low (248 bits) parts.
///   h = src[3] >> 56;
///   l = src[3] & 0x00FFFFFFFFFFFFFF;
///   // 5*2^248 = 1 mod q; hence, we add floor(h/5) + (h mod 5)*2^248
///   // to the low part.
///   quo = (h * 0xCD) >> 10;
///   rem = h - (5 * quo);
///   cc = add_carry(0, src[0], quo, &out[0]);
///   cc = add_carry(cc, src[1], 0, &out[1]);
///   cc = add_carry(cc, src[2], 0, &out[2]);
///   (void)add_carry(cc, l, rem << 56, &out[3]);
/// }
/// ```
///
/// The level-1 prime's defining identity `5 * 2^248 == 1 mod p` is the
/// crux: any value `(h << 248) + l` (where `h` is the top 8 bits and `l`
/// is the low 248 bits) reduces to `floor(h / 5) + (h mod 5) * 2^248 + l`
/// mod `p`. The `(h * 0xCD) >> 10` is a fixed-point reciprocal estimate
/// for `floor(h / 5)` exact for `h < 256` (`0xCD == 205 ≈ 1024 / 5`).
///
/// This helper operates on a 4-limb 256-bit plain integer (NOT the
/// 5-limb Montgomery `Fp`). It is internal to [`fp_decode_reduce`]'s
/// per-block fold: a 32-byte input block is decoded into four `u64`s
/// little-endian, [`partial_reduce`]'d to a 248-bit equivalent, and then
/// re-encoded to 32 bytes for [`fp_decode`] (which expects a canonical
/// value below `p`).
///
/// In-place is permitted: `out` may alias `src` (the reference's
/// `partial_reduce(t, t)` exercises this; the per-element add-carry
/// chain reads `src[i]` before writing `out[i]`, so aliasing is safe).
/// The port handles aliasing via a four-limb snapshot copy (Rust's
/// borrow checker would reject the literal aliased call).
fn partial_reduce(out: &mut [u64; 4], src: &[u64; 4]) {
    let snap = *src;
    let h: u64 = snap[3] >> 56;
    let l: u64 = snap[3] & 0x00FF_FFFF_FFFF_FFFF;
    let quo: u64 = (h.wrapping_mul(0xCD)) >> 10;
    let rem: u64 = h.wrapping_sub(5u64.wrapping_mul(quo));
    let mut cc = add_carry(0, snap[0], quo, &mut out[0]);
    cc = add_carry(cc, snap[1], 0, &mut out[1]);
    cc = add_carry(cc, snap[2], 0, &mut out[2]);
    let _ = add_carry(cc, l, rem << 56, &mut out[3]);
}

/// Little-endian decoding of an 8-byte slice as a `u64`. Mirrors the
/// reference's
/// `static inline uint64_t dec64le(const void *src)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:743..750` exactly.
/// `u64::from_le_bytes` is the same bit pattern; the port spells out
/// the shift-and-OR chain in the same form the reference does so the
/// source-level correspondence is visible at a glance.
fn dec64le(src: &[u8]) -> u64 {
    (src[0] as u64)
        | ((src[1] as u64) << 8)
        | ((src[2] as u64) << 16)
        | ((src[3] as u64) << 24)
        | ((src[4] as u64) << 32)
        | ((src[5] as u64) << 40)
        | ((src[6] as u64) << 48)
        | ((src[7] as u64) << 56)
}

/// Little-endian encoding of a `u64` as an 8-byte slice. Mirrors the
/// reference's
/// `static inline void enc64le(void *dst, uint64_t x)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:729..741` exactly.
/// `x.to_le_bytes()` is the same bit pattern; the port spells out the
/// individual byte writes in the same form the reference does so the
/// source-level correspondence is visible at a glance.
fn enc64le(dst: &mut [u8], x: u64) {
    dst[0] = x as u8;
    dst[1] = (x >> 8) as u8;
    dst[2] = (x >> 16) as u8;
    dst[3] = (x >> 24) as u8;
    dst[4] = (x >> 32) as u8;
    dst[5] = (x >> 40) as u8;
    dst[6] = (x >> 48) as u8;
    dst[7] = (x >> 56) as u8;
}

/// GF(p) deserialization-and-reduction from an arbitrary-length byte
/// slice. Mirrors the reference's
/// `void fp_decode_reduce(fp_t *d, const void *src, size_t len)` at
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:752..791` exactly:
///
/// ```c
/// void fp_decode_reduce(fp_t *d, const void *src, size_t len) {
///   uint64_t t[4];   // Stores Nbytes * 8 bits
///   uint8_t tmp[32]; // Nbytes
///   const uint8_t *b = src;
///   fp_set_zero(d);
///   if (len == 0) {
///     return;
///   }
///   size_t rem = len % 32;
///   if (rem != 0) {
///     // Input size is not a multiple of 32, we decode a partial
///     // block, which is already less than 2^248.
///     size_t k = len - rem;
///     memcpy(tmp, b + k, len - k);
///     memset(tmp + len - k, 0, (sizeof tmp) - (len - k));
///     fp_decode(d, tmp);
///     len = k;
///   }
///   // Process all remaining blocks, in descending address order.
///   while (len > 0) {
///     fp_mul(d, d, &R2);
///     len -= 32;
///     t[0] = dec64le(b + len);
///     t[1] = dec64le(b + len + 8);
///     t[2] = dec64le(b + len + 16);
///     t[3] = dec64le(b + len + 24);
///     partial_reduce(t, t);
///     enc64le(tmp, t[0]);
///     enc64le(tmp + 8, t[1]);
///     enc64le(tmp + 16, t[2]);
///     enc64le(tmp + 24, t[3]);
///     fp_t a;
///     fp_decode(&a, tmp);
///     fp_add(d, d, &a);
///   }
/// }
/// ```
///
/// Two-phase reduction:
/// 1. **Partial-block prefix.** If `len % 32 != 0`, the trailing partial
///    block (at most 31 bytes) is decoded into `d` via [`fp_decode`]
///    after zero-padding to 32 bytes. The partial block's value is below
///    `2^248 < p`, so the canonical-range check inside [`fp_decode`]
///    always succeeds and `d` is its Montgomery representative.
/// 2. **Full-block descending fold.** Each preceding 32-byte block is
///    decoded via [`partial_reduce`] (the level-1 prime's `5 * 2^248 ==
///    1 mod p` identity is the reduction kernel), re-encoded to 32 bytes
///    in canonical form, decoded via [`fp_decode`], and added to `d`
///    *after* `d` has been multiplied by `R2 == 2^256 mod p` (which
///    shifts the running accumulator up by one block's worth of
///    positional bits).
///
/// The `R2` constant carries the Montgomery factor: positionally, `d *
/// R2 == d * 2^256 mod p`, exactly the block-shift the descending fold
/// needs. Combined with the partial-block prefix, the final `d` is the
/// Montgomery representative of the full input's residue mod `p`. The
/// reference's [`fp_decode`] call discards its return value (the
/// canonical-range check is guaranteed by the partial-block bound and
/// by the [`partial_reduce`] guarantee that the re-encoded block is also
/// below `p`); the port follows by ignoring the return.
pub fn fp_decode_reduce(d: &mut Fp, src: &[u8]) {
    let mut tmp: [u8; 32] = [0u8; 32];
    fp_set_zero(d);
    let mut len = src.len();
    if len == 0 {
        return;
    }
    let rem = len % 32;
    if rem != 0 {
        let k = len - rem;
        tmp[..(len - k)].copy_from_slice(&src[k..len]);
        for byte in tmp.iter_mut().skip(len - k) {
            *byte = 0;
        }
        let _ = fp_decode(d, &tmp);
        len = k;
    }
    while len > 0 {
        let mut prod: Fp = [0u64; NWORDS_FIELD];
        fp_mul(&mut prod, d, &R2);
        *d = prod;
        len -= 32;
        let mut t: [u64; 4] = [0u64; 4];
        t[0] = dec64le(&src[len..(len + 8)]);
        t[1] = dec64le(&src[(len + 8)..(len + 16)]);
        t[2] = dec64le(&src[(len + 16)..(len + 24)]);
        t[3] = dec64le(&src[(len + 24)..(len + 32)]);
        let snap = t;
        partial_reduce(&mut t, &snap);
        enc64le(&mut tmp[0..8], t[0]);
        enc64le(&mut tmp[8..16], t[1]);
        enc64le(&mut tmp[16..24], t[2]);
        enc64le(&mut tmp[24..32], t[3]);
        let mut a: Fp = [0u64; NWORDS_FIELD];
        let _ = fp_decode(&mut a, &tmp);
        let mut sum: Fp = [0u64; NWORDS_FIELD];
        fp_add(&mut sum, d, &a);
        *d = sum;
    }
}

// ---------------------------------------------------------------------------
// GF(p^2) layer. Mirrors `vendor/the-sqisign/src/gf/ref/lvlx/fp2.c`.
//
// The level-1 quadratic extension is `Fp[X]/(X^2 + 1)`, the reference's own
// comment at the top of `fp2.c` ("Arithmetic modulo X^2 + 1"). An element is
// the C struct `fp2_t { fp_t re, im; }`; the Rust mirror keeps the same
// field naming (`re`, `im`) and the same memory order so the bit-for-bit
// boundary correspondence stays visible at the call site.
// ---------------------------------------------------------------------------

/// Byte length of the canonical `fp2` serialization, the
/// `FP2_ENCODED_BYTES` constant from
/// `vendor/the-sqisign/src/precomp/ref/lvl1/include/encoded_sizes.h`.
/// It is exactly `2 * FP_ENCODED_BYTES` (64 = 2 * 32) for this level: the
/// reference's `fp2_encode` writes `re` then `im` back to back.
pub const FP2_ENCODED_BYTES: usize = 64;

/// A level-1 GF(p^2) field element. Mirrors the C `struct fp2_t { fp_t
/// re, im; }`; the field naming (`re`, `im`) and the storage order are
/// preserved so the differential vector records can pin which `fp_t` half
/// is which. Like [`Fp`], each component is the redundant radix-2^51
/// representation: equality at this layer is therefore *not* raw struct
/// equality, it is [`fp2_is_equal`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Fp2 {
    pub re: Fp,
    pub im: Fp,
}

/// Alias matching the reference's typedef `fp2_t`; the C uses the
/// `_t`-suffixed spelling at every call site, so the alias is exposed for
/// symmetry with [`crate::Fp`] (which mirrors `fp_t`).
pub type Fp2t = Fp2;

/// `fp2_set_small(x, val)`: write the scalar `val` into the real part and
/// zero the imaginary part. Wraps the already-ported [`fp_set_small`] and
/// [`fp_set_zero`]; the same int32-narrowing the [`fp_set_small`] wrapper
/// applies to `val` is observable at this boundary (the imaginary part is
/// not affected by `val`).
pub fn fp2_set_small(x: &mut Fp2, val: u64) {
    fp_set_small(&mut x.re, val);
    fp_set_zero(&mut x.im);
}

/// `fp2_mul_small(x, y, n)`: multiply BOTH components of `y` by the
/// scalar `n` and write into `x`. Distinct shape from [`fp2_set_small`]:
/// `n` is consumed by both halves rather than just the real part. The
/// same int32-narrowing the per-half [`fp_mul_small`] applies fires
/// twice, once per component.
pub fn fp2_mul_small(x: &mut Fp2, y: &Fp2, n: u32) {
    fp_mul_small(&mut x.re, &y.re, n);
    fp_mul_small(&mut x.im, &y.im, n);
}

/// `fp2_set_one(x)`: write the canonical multiplicative identity, which
/// is `(1, 0)` in `Fp[X]/(X^2 + 1)`. The real part receives the
/// Montgomery representative of `1`; the imaginary part is zeroed.
pub fn fp2_set_one(x: &mut Fp2) {
    fp_set_one(&mut x.re);
    fp_set_zero(&mut x.im);
}

/// `fp2_set_zero(x)`: write the canonical additive identity. Both
/// components are zeroed.
pub fn fp2_set_zero(x: &mut Fp2) {
    fp_set_zero(&mut x.re);
    fp_set_zero(&mut x.im);
}

/// `fp2_is_zero(a)`: branchless `uint32_t` predicate returning the
/// all-ones mask iff both components reduce to the field zero. The
/// reference ANDs the two per-component masks, propagating the
/// `{0, 0xFFFFFFFF}` convention; the port mirrors that exactly.
pub fn fp2_is_zero(a: &Fp2) -> u32 {
    fp_is_zero(&a.re) & fp_is_zero(&a.im)
}

/// `fp2_is_equal(a, b)`: branchless `uint32_t` predicate returning the
/// all-ones mask iff both `re` and `im` are componentwise equal modulo
/// `p`. As for [`fp2_is_zero`], the two per-component masks are ANDed.
pub fn fp2_is_equal(a: &Fp2, b: &Fp2) -> u32 {
    fp_is_equal(&a.re, &b.re) & fp_is_equal(&a.im, &b.im)
}

/// `fp2_is_one(a)`: branchless `uint32_t` predicate returning the
/// all-ones mask iff `a == (1, 0)`. The reference compares the real
/// part against the `extern const ONE` Montgomery representative
/// directly via `fp_is_equal`, and ANDs that with the imaginary
/// part's `fp_is_zero`. The port mirrors that exactly using
/// [`MONTGOMERY_ONE`] (the same bit pattern the reference's `ONE`
/// holds, pinned by the `nres_of_positional_one_is_montgomery_one`
/// test in this crate).
pub fn fp2_is_one(a: &Fp2) -> u32 {
    fp_is_equal(&a.re, &MONTGOMERY_ONE) & fp_is_zero(&a.im)
}

/// `fp2_copy(x, y)`: componentwise copy. The same redundancy
/// considerations [`fp_copy`] documents apply: this is a plain
/// per-component limb assignment, NOT a canonicalisation.
pub fn fp2_copy(x: &mut Fp2, y: &Fp2) {
    fp_copy(&mut x.re, &y.re);
    fp_copy(&mut x.im, &y.im);
}

/// `fp2_add(x, y, z)`: componentwise addition.
pub fn fp2_add(x: &mut Fp2, y: &Fp2, z: &Fp2) {
    fp_add(&mut x.re, &y.re, &z.re);
    fp_add(&mut x.im, &y.im, &z.im);
}

/// `fp2_add_one(x, y)`: `x = y + 1`. The real part receives
/// `y.re + ONE` (where `ONE` is the Montgomery representative of `1`,
/// matching the reference's `&ONE` constant). The imaginary part is
/// copied through unchanged via [`fp_copy`], NOT re-added with zero;
/// this preserves the exact non-canonical limb pattern of `y.im`.
pub fn fp2_add_one(x: &mut Fp2, y: &Fp2) {
    fp_add(&mut x.re, &y.re, &MONTGOMERY_ONE);
    fp_copy(&mut x.im, &y.im);
}

/// `fp2_sub(x, y, z)`: componentwise subtraction.
pub fn fp2_sub(x: &mut Fp2, y: &Fp2, z: &Fp2) {
    fp_sub(&mut x.re, &y.re, &z.re);
    fp_sub(&mut x.im, &y.im, &z.im);
}

/// `fp2_neg(x, y)`: componentwise negation.
pub fn fp2_neg(x: &mut Fp2, y: &Fp2) {
    fp_neg(&mut x.re, &y.re);
    fp_neg(&mut x.im, &y.im);
}

/// `fp2_mul(x, y, z)`: multiplication in `Fp[X]/(X^2 + 1)`. The
/// reference uses the Karatsuba-style three-multiplication identity
///
/// ```text
/// (y.re + y.im * i) * (z.re + z.im * i)
///   = (y.re * z.re - y.im * z.im)
///   + ((y.re + y.im) * (z.re + z.im) - y.re * z.re - y.im * z.im) * i
/// ```
///
/// computed in the exact statement order of `fp2.c:95..107`:
///
/// 1. `t0 = y.re + y.im`
/// 2. `t1 = z.re + z.im`
/// 3. `t0 = t0 * t1` (sum of cross terms plus the two diagonal products)
/// 4. `t1 = y.im * z.im` (one diagonal)
/// 5. `x.re = y.re * z.re` (the other diagonal, stored into the
///    destination's real part)
/// 6. `x.im = t0 - t1` (sum of cross terms plus `y.re * z.re`)
/// 7. `x.im = x.im - x.re` (cross terms only: this is the imaginary part)
/// 8. `x.re = x.re - t1` (diagonal difference: this is the real part)
///
/// The port transcribes those eight statements identifier for
/// identifier; an in-place write to `x.re` would alias `y.re` or
/// `z.re` if `x` overlaps an input, so the reference's structure (write
/// `x.re` only at step 5, then read it again at step 7) is preserved.
pub fn fp2_mul(x: &mut Fp2, y: &Fp2, z: &Fp2) {
    let mut t0: Fp = [0u64; NWORDS_FIELD];
    let mut t1: Fp = [0u64; NWORDS_FIELD];

    fp_add(&mut t0, &y.re, &y.im);
    fp_add(&mut t1, &z.re, &z.im);
    let snap_t0 = t0;
    fp_mul(&mut t0, &snap_t0, &t1);
    fp_mul(&mut t1, &y.im, &z.im);
    fp_mul(&mut x.re, &y.re, &z.re);
    let snap_xre = x.re;
    fp_sub(&mut x.im, &t0, &t1);
    let snap_xim = x.im;
    fp_sub(&mut x.im, &snap_xim, &snap_xre);
    let snap_xre2 = x.re;
    fp_sub(&mut x.re, &snap_xre2, &t1);
}

/// `fp2_sqr(x, y)`: squaring in `Fp[X]/(X^2 + 1)`. The reference uses
/// the standard difference-of-squares identity
///
/// ```text
/// (y.re + y.im * i) ^ 2
///   = (y.re + y.im) * (y.re - y.im) + 2 * y.re * y.im * i
/// ```
///
/// computed in the exact statement order of `fp2.c:110..119`:
///
/// 1. `sum  = y.re + y.im`
/// 2. `diff = y.re - y.im`
/// 3. `x.im = y.re * y.im`
/// 4. `x.im = x.im + x.im`         (`2 * y.re * y.im`)
/// 5. `x.re = sum * diff`          (`y.re^2 - y.im^2`)
pub fn fp2_sqr(x: &mut Fp2, y: &Fp2) {
    let mut sum: Fp = [0u64; NWORDS_FIELD];
    let mut diff: Fp = [0u64; NWORDS_FIELD];

    fp_add(&mut sum, &y.re, &y.im);
    fp_sub(&mut diff, &y.re, &y.im);
    fp_mul(&mut x.im, &y.re, &y.im);
    let snap_xim = x.im;
    fp_add(&mut x.im, &snap_xim, &snap_xim);
    fp_mul(&mut x.re, &sum, &diff);
}

/// `fp2_inv(x)`: in-place multiplicative inverse. The reference computes
///
/// ```text
/// d   = re^2 + im^2
/// x   = (re / d, -im / d)
/// ```
///
/// in the exact statement order of `fp2.c:122..133`:
///
/// 1. `t0 = re^2`
/// 2. `t1 = im^2`
/// 3. `t0 = t0 + t1`              (norm)
/// 4. `t0 = t0^-1`                (inverse of the norm)
/// 5. `re  = re * t0`
/// 6. `im  = im * t0`
/// 7. `im  = -im`
///
/// On `x == (0, 0)` the norm is zero, [`fp_inv`] squares-and-multiplies
/// down to the canonical zero (see [`fp_inv`]'s comment), and the
/// resulting `x` is the canonical zero pair, matching the reference.
pub fn fp2_inv(x: &mut Fp2) {
    let mut t0: Fp = [0u64; NWORDS_FIELD];
    let mut t1: Fp = [0u64; NWORDS_FIELD];

    fp_sqr(&mut t0, &x.re);
    fp_sqr(&mut t1, &x.im);
    let snap_t0 = t0;
    fp_add(&mut t0, &snap_t0, &t1);
    fp_inv(&mut t0);
    let snap_xre = x.re;
    fp_mul(&mut x.re, &snap_xre, &t0);
    let snap_xim = x.im;
    fp_mul(&mut x.im, &snap_xim, &t0);
    let snap_xim2 = x.im;
    fp_neg(&mut x.im, &snap_xim2);
}

/// `fp2_is_square(x)`: branchless `uint32_t` predicate returning the
/// all-ones mask iff `x` is a quadratic residue in `Fp[X]/(X^2 + 1)`.
/// The reference uses the norm criterion: `x` is a square in `Fp^2` iff
/// `re^2 + im^2` is a square in `Fp` (the standard reduction via the
/// `Fp^2 -> Fp` norm). Three statements:
///
/// 1. `t0 = re^2`
/// 2. `t1 = im^2`
/// 3. `t0 = t0 + t1`
/// 4. return `fp_is_square(t0)`
pub fn fp2_is_square(x: &Fp2) -> u32 {
    let mut t0: Fp = [0u64; NWORDS_FIELD];
    let mut t1: Fp = [0u64; NWORDS_FIELD];

    fp_sqr(&mut t0, &x.re);
    fp_sqr(&mut t1, &x.im);
    let snap = t0;
    fp_add(&mut t0, &snap, &t1);
    fp_is_square(&t0)
}

/// `fp2_sqrt(a)`: in-place square root, from Aardal et al, eprint
/// 2024/1563 ("Optimized One-Dimensional SQIsign Verification on Intel
/// and Cortex-M4"). The transcription preserves the reference's
/// `fp2.c:148..202` statement order identifier-for-identifier; see the
/// reference comment block for the algebra. The choice of representative
/// (negate-when-odd-or-zero-and-other-odd) is canonicalised via
/// [`fp_encode`] on `t0` and `t1` and reading the parity off the
/// low byte, exactly as the reference does.
pub fn fp2_sqrt(a: &mut Fp2) {
    let mut x0: Fp = [0u64; NWORDS_FIELD];
    let mut x1: Fp = [0u64; NWORDS_FIELD];
    let mut t0: Fp = [0u64; NWORDS_FIELD];
    let mut t1: Fp = [0u64; NWORDS_FIELD];

    // x0 = delta = sqrt(a.re^2 + a.im^2)
    fp_sqr(&mut x0, &a.re);
    fp_sqr(&mut x1, &a.im);
    let snap_x0 = x0;
    fp_add(&mut x0, &snap_x0, &x1);
    fp_sqrt(&mut x0);
    // If a.im == 0, restore delta = a.re.
    let snap_x0b = x0;
    fp_select(&mut x0, &snap_x0b, &a.re, fp_is_zero(&a.im));
    // x0 = delta + a.re; t0 = 2 * x0.
    let snap_x0c = x0;
    fp_add(&mut x0, &snap_x0c, &a.re);
    fp_add(&mut t0, &x0, &x0);
    // x1 = t0^((p-3)/4)
    fp_exp3div4(&mut x1, &t0);
    // x0 = x0 * x1, x1 = x1 * a.im, t1 = (2*x0)^2.
    let snap_x0d = x0;
    fp_mul(&mut x0, &snap_x0d, &x1);
    let snap_x1 = x1;
    fp_mul(&mut x1, &snap_x1, &a.im);
    fp_add(&mut t1, &x0, &x0);
    let snap_t1 = t1;
    fp_sqr(&mut t1, &snap_t1);
    // If t1 == t0, return x0 + x1*i; otherwise x1 - x0*i.
    let snap_t0b = t0;
    fp_sub(&mut t0, &snap_t0b, &t1);
    let f = fp_is_zero(&t0);
    fp_neg(&mut t1, &x0);
    fp_copy(&mut t0, &x1);
    let snap_t0c = t0;
    fp_select(&mut t0, &snap_t0c, &x0, f);
    let snap_t1b = t1;
    fp_select(&mut t1, &snap_t1b, &x1, f);

    let t0_is_zero = fp_is_zero(&t0);

    let mut tmp_bytes = [0u8; FP_ENCODED_BYTES];
    fp_encode(&mut tmp_bytes, &t0);
    let t0_is_odd = 0u32.wrapping_sub((tmp_bytes[0] as u32) & 1);
    fp_encode(&mut tmp_bytes, &t1);
    let t1_is_odd = 0u32.wrapping_sub((tmp_bytes[0] as u32) & 1);

    let negate_output = t0_is_odd | (t0_is_zero & t1_is_odd);
    fp_neg(&mut x0, &t0);
    fp_select(&mut a.re, &t0, &x0, negate_output);
    fp_neg(&mut x0, &t1);
    fp_select(&mut a.im, &t1, &x0, negate_output);
}

/// `fp2_sqrt_verify(a)`: in-place square root with a verification
/// check. The reference saves `a` to `t0`, computes `fp2_sqrt(a)`,
/// squares the result into `t1`, and returns `fp2_is_equal(&t0, &t1)`.
/// The post-call value of `a` is therefore the *square root* (or
/// whatever the sqrt routine produced on a non-square, which is not
/// guaranteed to satisfy the verify), and the returned `uint32_t` is
/// `0xFFFFFFFF` exactly when the input was a square (so the squaring
/// reproduces the input). On a non-square input, the post-call `a` is
/// not generally `0`; it is whatever `fp2_sqrt` deposited, and the
/// returned mask is `0`.
pub fn fp2_sqrt_verify(a: &mut Fp2) -> u32 {
    let t0 = *a;
    fp2_sqrt(a);
    let mut t1 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_sqr(&mut t1, a);
    fp2_is_equal(&t0, &t1)
}

/// `fp2_half(x, y)`: componentwise halving (multiplication by
/// `2^-1 mod p`), routed through the already-ported [`fp_half`] twice.
pub fn fp2_half(x: &mut Fp2, y: &Fp2) {
    fp_half(&mut x.re, &y.re);
    fp_half(&mut x.im, &y.im);
}

/// `fp2_encode(dst, a)`: write 64 bytes, the canonical serialization of
/// `a`. The first 32 bytes are [`fp_encode`] of `a.re`, the second 32
/// are [`fp_encode`] of `a.im`. The reference uses `void *dst` with
/// pointer arithmetic; the port takes the same 64-byte buffer as a
/// fixed-size array.
pub fn fp2_encode(dst: &mut [u8; FP2_ENCODED_BYTES], a: &Fp2) {
    let (re_dst, im_dst) = dst.split_at_mut(FP_ENCODED_BYTES);
    let re_arr: &mut [u8; FP_ENCODED_BYTES] = re_dst.try_into().unwrap();
    let im_arr: &mut [u8; FP_ENCODED_BYTES] = im_dst.try_into().unwrap();
    fp_encode(re_arr, &a.re);
    fp_encode(im_arr, &a.im);
}

/// `fp2_decode(d, src)`: read 64 bytes, write `d` and return the
/// canonical-range mask `re_mask & im_mask`. Both per-half decodes must
/// be in range for the combined result to be `0xFFFFFFFF`; either out
/// of range, the combined mask is `0`. The reference returns `re & im`
/// (bitwise AND of the two `uint32_t` per-half masks), the port mirrors
/// that exactly. The decoded `d.re` and `d.im` may carry the
/// out-of-range zeroing applied per-half by [`fp_decode`] (see
/// [`fp_decode`]'s ANDing of `res` into every limb), so the output is
/// only meaningful on the in-range branch.
pub fn fp2_decode(d: &mut Fp2, src: &[u8; FP2_ENCODED_BYTES]) -> u32 {
    let (re_src, im_src) = src.split_at(FP_ENCODED_BYTES);
    let re_arr: &[u8; FP_ENCODED_BYTES] = re_src.try_into().unwrap();
    let im_arr: &[u8; FP_ENCODED_BYTES] = im_src.try_into().unwrap();
    let re = fp_decode(&mut d.re, re_arr);
    let im = fp_decode(&mut d.im, im_arr);
    re & im
}

/// `fp2_select(d, a0, a1, ctl)`: branchless componentwise conditional
/// select. The same `ctl` contract [`fp_select`] documents applies: the
/// reference restricts `ctl` to `0x00000000` (select `a0`) or
/// `0xFFFFFFFF` (select `a1`); any other `ctl` is undefined.
pub fn fp2_select(d: &mut Fp2, a0: &Fp2, a1: &Fp2, ctl: u32) {
    fp_select(&mut d.re, &a0.re, &a1.re, ctl);
    fp_select(&mut d.im, &a0.im, &a1.im, ctl);
}

/// `fp2_cswap(a, b, ctl)`: branchless componentwise conditional swap.
/// The same `ctl & 1` LSB-only contract [`fp_cswap`] documents applies.
pub fn fp2_cswap(a: &mut Fp2, b: &mut Fp2, ctl: u32) {
    fp_cswap(&mut a.re, &mut b.re, ctl);
    fp_cswap(&mut a.im, &mut b.im, ctl);
}

/// `FP_ENCODED_BYTES`: byte length of the canonical `fp` serialization,
/// matching the reference's per-level constant in
/// `encoded_sizes.h`. Hardcoded to 32 here because [`fp_encode`] /
/// [`fp_decode`] are typed against `[u8; 32]` directly; the `fp2`
/// boundaries above use this constant to slice the 64-byte fp2 buffer
/// into its two halves.
const FP_ENCODED_BYTES: usize = 32;

/// `fp2_batched_inv(x, len)`: Montgomery's batched inverse. Replaces
/// each element of `x` with its inverse in place; the cost is one
/// `fp2_inv` plus `3 * (len - 1)` `fp2_mul` calls (instead of the
/// `len` `fp2_inv` calls a naive loop would do). The reference uses
/// C99 VLAs `t1[len], t2[len]`; the port uses heap `Vec`s of the
/// same length for the same prefix-product / suffix-product
/// scratchpads. Statement order matches `fp2.c:223..251` exactly.
///
/// A differential boundary for this routine would require a new
/// emitter shape (variable-length list of fp2 in / variable-length
/// list of fp2 out); pinning equivalence to the reference for a
/// composite chain of [`fp2_mul`] and [`fp2_inv`] adds no new
/// algebraic content over those primitives' own per-call differential
/// vectors. Left out of the fp2 mega-batch differential gate
/// deliberately; the primitives it composes are pinned bit-for-bit.
pub fn fp2_batched_inv(x: &mut [Fp2]) {
    let len = x.len();
    if len == 0 {
        return;
    }
    let zero_fp2 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1: Vec<Fp2> = vec![zero_fp2; len];
    let mut t2: Vec<Fp2> = vec![zero_fp2; len];

    // t1 = x0, x0*x1, ..., x0 * x1 * ... * xn
    fp2_copy(&mut t1[0], &x[0]);
    for i in 1..len {
        let prev = t1[i - 1];
        fp2_mul(&mut t1[i], &prev, &x[i]);
    }

    // inverse = 1 / (x0 * x1 * ... * xn)
    let mut inverse = t1[len - 1];
    fp2_inv(&mut inverse);

    fp2_copy(&mut t2[0], &inverse);
    // t2 = 1 / (x0 * x1 * ... * xn), 1 / (x0 * x1 * ... * x(n-1)), ..., 1 / x0
    for i in 1..len {
        let prev = t2[i - 1];
        fp2_mul(&mut t2[i], &prev, &x[len - i]);
    }

    fp2_copy(&mut x[0], &t2[len - 1]);

    for i in 1..len {
        let lhs = t1[i - 1];
        let rhs = t2[len - i - 1];
        fp2_mul(&mut x[i], &lhs, &rhs);
    }
}

/// `fp2_pow_vartime(out, x, exp, size)`: square-and-multiply
/// exponentiation. **Not constant time** (the reference comment is
/// "Warning!! Not constant time!"). `exp` is `size` little-endian
/// `digit_t == u64` words; the loop walks each word low to high and
/// each bit within a word low to high, multiplying when the bit is
/// set and squaring on every step. Statement order mirrors
/// `fp2.c:256..275`.
///
/// As for [`fp2_batched_inv`], a differential boundary for this
/// composite over [`fp2_mul`] and [`fp2_sqr`] adds no new algebraic
/// content. Left out of the fp2 mega-batch differential gate
/// deliberately.
pub fn fp2_pow_vartime(out: &mut Fp2, x: &Fp2, exp: &[u64]) {
    let mut acc = *x;
    fp2_set_one(out);

    // The reference loops `for (int i = 0; i < RADIX; i++)` where
    // `RADIX == 64` (the digit_t bit-width set by `tutil.h` under
    // `RADIX_64`, not the field's per-limb radix of 51 which is a
    // different `RADIX` macro in the `fp_p5248_64.h` namespace). The
    // loop walks the 64 bits of each digit low to high.
    const DIGIT_BITS: u32 = 64;

    for &word in exp {
        for i in 0..DIGIT_BITS {
            let bit = (word >> i) & 1;
            if bit == 1 {
                let snap_out = *out;
                fp2_mul(out, &snap_out, &acc);
            }
            let snap_acc = acc;
            fp2_sqr(&mut acc, &snap_acc);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Independent value-level oracle for [`NRES_C`]: applying [`nres`]
    /// to the positional `1` must yield the Montgomery representative of
    /// `1`. This is the open question recorded in
    /// `~/sqisign-port-notes.md` ("Boundaries that hardcode a precomputed
    /// reference constant"): [`fp_set_one`] writes [`MONTGOMERY_ONE`]
    /// directly as a const, deferring the algorithmic confirmation until
    /// `nres` (and the `R^2` constant it relies on) landed. With both
    /// ported as internal helpers alongside [`fp_set_small`], the
    /// independent check is now executable: running the positional
    /// `[1, 0, 0, 0, 0]` through `nres` (which is exactly
    /// `modmul(_, NRES_C, _)`) must produce the same bit pattern the
    /// reference exposes as `extern const ONE`. A mismatch would mean
    /// either [`NRES_C`] is wrong limb-for-limb or [`MONTGOMERY_ONE`] is.
    #[test]
    fn nres_of_positional_one_is_montgomery_one() {
        let positional_one: Fp = [1, 0, 0, 0, 0];
        let mut out: Fp = [0u64; NWORDS_FIELD];
        nres(&positional_one, &mut out);
        assert_eq!(
            out, MONTGOMERY_ONE,
            "nres(positional 1) must equal the Montgomery representative of 1; \
             a mismatch indicates NRES_C or MONTGOMERY_ONE was transcribed wrong"
        );
    }
}
