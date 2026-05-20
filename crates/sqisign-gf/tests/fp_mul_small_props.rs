//! Property tests for `fp_mul_small`.
//!
//! `fp_mul_small(out, a, val)` is the reference's thin wrapper around
//! `modmli(*a, (int)val, *x)`: narrow `val` to `int32`, build the
//! Montgomery representative of the narrowed-and-sign-extended integer
//! via `modint(b, t)` into a five-limb scratch, then run a single
//! `modmul(a, t, x)` cross product. The output is the Montgomery
//! representative of `a * (val as i32 as i64)` in the redundant
//! radix-2^51 representation, with limbs 0..=3 below `2^51` (the
//! column mask) and limb 4 fully unmasked (the underlying `modmul`'s
//! final `c[4] = (spint)t` is a 64-bit truncation, no `& mask`),
//! exactly as the reference leaves it. Raw-limb equality is therefore
//! *not* a sound notion of field equality in general; the properties
//! below were each cross-checked bit-exactly against the full
//! 1132-vector C-derived battery before being committed (the canonical
//! correctness check remains the differential test in
//! `fp_mul_small_vectors.rs`).
//!
//! Four sound raw-limb properties:
//!
//! 1. **`fp_mul_small(out, a, 0)` is the canonical all-zero
//!    representative, for arbitrary `a`.** Verified empirically across
//!    all 12 val == 0 records in the differential battery before
//!    pinning: every column sum inside the inner `modmul(a, [0; 5])`
//!    accumulates only zero partial products (each `a[i] * 0` is
//!    zero), so the per-column 51-bit mask reads off zero into
//!    `v0..v4` and the final per-limb writes `c[0..=3]` are zero; the
//!    final `c[4] = (spint)t` truncates a zero `t` to zero. The
//!    Montgomery image of zero is zero, and so the redundant
//!    representative this boundary produces is the *canonical*
//!    all-zero limb vector, bit-exact, rather than merely a redundant
//!    representative congruent to zero.
//!
//! 2. **`fp_mul_small(out, a, 1) ==_field a` via [`fp_is_equal`].**
//!    Multiplication by `1` in the field is the identity. The
//!    boundary's actual computation is `a * MONTGOMERY_ONE * R^-1
//!    mod p == a mod p`, so the output is a valid (and in general
//!    distinct from `a`) redundant representative of the same residue
//!    class. The sound value-level oracle is therefore [`fp_is_equal`]
//!    (the reference's `modcmp`, which `redc`s both operands to their
//!    canonical reduced form before comparing), not raw-limb equality:
//!    two redundant representatives of the same residue class are
//!    bit-distinct in general but always compare equal under
//!    [`fp_is_equal`].
//!
//! 3. **Cross-oracle equality:** `fp_mul_small(out, a, val) ==
//!    fp_mul(a, &fp_set_small(val as i32 as u64))` bit-exact, for
//!    arbitrary `a` and `val`. The boundary is literally that
//!    `modint + modmul` chain: `modmli(a, b, c)` does `modint(b, t);
//!    modmul(a, t, c);`, and `fp_set_small(val as i32 as u64)`
//!    re-runs the same `modint` (because the public wrapper's `val as
//!    i32 as u64` already canonicalises any 64-bit input to the same
//!    `(int32_t)val` image the inner `modmli`'s `(int)val` cast
//!    produces, sign-extended at the limb-0 write). The pre-canonical
//!    `val as i32 as u64` on the right-hand side is load-bearing: the
//!    cross-oracle would otherwise mismatch on `val` values whose
//!    sign-extension differs from a naive `val as u64` widening
//!    (every val above `2^31 - 1` exercises that branch). Verified
//!    empirically across all 1132 records before pinning.
//!
//! 4. **Structural carry-propagation invariant on limbs 0..=3 for
//!    arbitrary inputs.** Inherited from `fp_mul`: the inner
//!    `modmul`'s per-column writes `c[0..=3] = (t as u64) & MASK51`
//!    apply the per-limb 51-bit mask, so every output has
//!    `out[0..4] < 2^51`. Limb 4 is intentionally *not* asserted: the
//!    reference's final write is the unmasked truncation `c[4] =
//!    (spint)t` and the port faithfully does the same.

use proptest::prelude::*;
use sqisign_gf::{fp_is_equal, fp_mul, fp_mul_small, fp_set_small, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;
const ZERO: Fp = [0u64; NWORDS_FIELD];

fn uniform5() -> impl Strategy<Value = Fp> {
    (
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
    )
        .prop_map(|(a, b, c, d, e)| [a, b, c, d, e])
}

proptest! {
    // (1) val == 0: output is the bit-exact canonical all-zero
    // representative, for arbitrary `a`. The cross product
    // `modmul(a, [0; 5])` accumulates only zero partial products by
    // structural induction on the per-column accumulator, so the
    // output is bit-exactly the canonical zero.
    #[test]
    fn val_zero_is_canonical_zero(a in uniform5()) {
        let mut out: Fp = [0u64; NWORDS_FIELD];
        fp_mul_small(&mut out, &a, 0);
        prop_assert_eq!(out, ZERO);
    }

    // (2) val == 1: output is field-equal to `a` via fp_is_equal.
    // Raw-limb equality is *not* asserted: the boundary's output is
    // the Montgomery-reduced redundant representative of `a`, which
    // in general differs limb-for-limb from `a` itself; only the
    // value-level equality holds.
    #[test]
    fn val_one_is_field_identity(a in uniform5()) {
        let mut out: Fp = [0u64; NWORDS_FIELD];
        fp_mul_small(&mut out, &a, 1);
        prop_assert_eq!(fp_is_equal(&out, &a), 0xffff_ffffu32);
    }

    // (3) Cross-oracle: fp_mul_small(a, val) == fp_mul(a, fp_set_small(val as i32 as u64)).
    // Sign-extend `val` into the u64 fp_set_small accepts so the
    // narrowed-int32 image is reproduced exactly; without that
    // sign-extension the right-hand side would mismatch on any val
    // above 2^31 - 1 (the C wrapper's `(int)val` cast turns those
    // into negative int32s, sign-extending to 0xffffffff_xxxxxxxx at
    // the positional limb-0 write).
    #[test]
    fn cross_oracle_matches_modint_then_modmul(
        a in uniform5(),
        val in any::<u32>(),
    ) {
        let val_sx = val as i32 as u64;
        let mut scratch: Fp = [0u64; NWORDS_FIELD];
        fp_set_small(&mut scratch, val_sx);
        let mut via_chain: Fp = [0u64; NWORDS_FIELD];
        fp_mul(&mut via_chain, &a, &scratch);

        let mut direct: Fp = [0u64; NWORDS_FIELD];
        fp_mul_small(&mut direct, &a, val);
        prop_assert_eq!(direct, via_chain);
    }

    // (4) Structural carry-propagation invariant: intermediate column
    // writes c[0..=3] are masked below 2^51; limb 4 is left unmasked
    // by design (the reference's final c[4] = (spint)t is a full
    // 64-bit truncation, no & mask) and is deliberately not
    // constrained.
    #[test]
    fn limbs_0_3_below_radix(
        a in uniform5(),
        val in any::<u32>(),
    ) {
        let mut out: Fp = [0u64; NWORDS_FIELD];
        fp_mul_small(&mut out, &a, val);
        for (k, &limb) in out.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }
}
