//! Property tests for `fp_set_small`.
//!
//! `fp_set_small(out, val)` is the reference's thin wrapper around
//! `modint((int)val, *x)`: narrow `val` to `int32`, sign-extend to `u64`
//! at limb 0, zero limbs 1..=4, then call `nres` to convert positional
//! to Montgomery. The output is the Montgomery representative of the
//! narrowed-and-sign-extended integer, in the redundant radix-2^51
//! representation.
//!
//! Three sound raw-limb properties:
//!
//! 1. **`fp_set_small(0)` is the canonical all-zero representative.**
//!    Verified empirically before pinning: the Montgomery image of
//!    positional zero is positional zero (`0 * R mod p == 0`), and the
//!    `nres` chain (multiply-then-modfsb) preserves bit-exactness on
//!    the all-zero input. The differential battery's val=0 records
//!    confirm this (`fp_set_small_vectors.rs` pins all of them as
//!    equal to the bit-exact zero), so the property is sound.
//! 2. **`fp_set_small(1)` equals `fp_set_one()` bit-exact.** The
//!    Montgomery representative of positional `1` is the
//!    [`MONTGOMERY_ONE`] constant `fp_set_one` writes directly. Pinned
//!    via `fp_is_equal` for value-level equality, *and* via raw-limb
//!    equality (the redundant form's lone canonical fixed point at
//!    val=1 is the bit-exact Montgomery ONE).
//! 3. **High-bits-ignored narrowing.** For arbitrary `u64` `val`,
//!    `fp_set_small(out, val)` equals `fp_set_small(out, val as i32 as
//!    u64)` bit-for-bit: the C wrapper's `(int)val` cast drops the high
//!    32 bits of `val`, so any two `val`s sharing the same low 32 bits
//!    produce the same output. The `val as i32 as u64` round-trip
//!    canonicalises `val` to its sign-extended-low-32-bits image, the
//!    same image the boundary reduces to before `nres`.
//!
//! Pre-fills are drawn from arbitrary `u64` limb tuples (a setter
//! accepts non-canonical destinations by construction). `val` is drawn
//! from `any::<u64>()` so the high-bits-ignored property is exercised
//! across the full 64-bit width.

use proptest::prelude::*;
use sqisign_gf::{fp_is_equal, fp_set_one, fp_set_small, Fp, NWORDS_FIELD};

/// Montgomery representative of `1`; must match the reference's
/// `extern const ONE` at `fp_p5248_64.c:526..530`.
const MONTGOMERY_ONE: Fp = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
];

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
    // representative, regardless of the destination's prior contents.
    // The Montgomery image of zero is zero, and `nres`'s `modmul`
    // chain on all-zero inputs produces all-zero outputs by structural
    // induction on the per-column accumulator (every partial product is
    // zero, every mask read is zero, every limb write is zero).
    #[test]
    fn val_zero_is_canonical_zero(prefill in uniform5()) {
        let mut out: Fp = prefill;
        fp_set_small(&mut out, 0);
        prop_assert_eq!(out, ZERO);
    }

    // (2) val == 1: output is the Montgomery ONE constant, bit-exact.
    // Pinned both via raw-limb equality (the redundant form's canonical
    // image of positional `1` IS the Montgomery ONE, by the
    // `nres_of_positional_one_is_montgomery_one` unit test) and via the
    // value-level `fp_is_equal` against `fp_set_one`'s output.
    #[test]
    fn val_one_equals_montgomery_one(prefill in uniform5()) {
        let mut from_small: Fp = prefill;
        fp_set_small(&mut from_small, 1);
        prop_assert_eq!(from_small, MONTGOMERY_ONE);

        let mut from_one: Fp = prefill;
        fp_set_one(&mut from_one);
        prop_assert_eq!(from_small, from_one);
        prop_assert_eq!(fp_is_equal(&from_small, &from_one), 0xffff_ffff);
    }

    // (3) High-bits-ignored narrowing: any two `val`s sharing the same
    // low 32 bits (after sign-extension) produce bit-exact equal
    // outputs. The narrowing is `(int)val` in the C wrapper, i.e. a
    // truncation-and-reinterpret to `int32_t` whose subsequent
    // assignment to `spint = uint64_t` sign-extends. The
    // `val as i32 as u64` round-trip canonicalises any `u64` to its
    // sign-extended-low-32-bits image, the same image the boundary
    // reduces to before `nres`.
    #[test]
    fn high_bits_ignored(prefill in uniform5(), val in any::<u64>()) {
        let narrowed = (val as i32) as u64;
        let mut from_full: Fp = prefill;
        fp_set_small(&mut from_full, val);
        let mut from_narrowed: Fp = prefill;
        fp_set_small(&mut from_narrowed, narrowed);
        prop_assert_eq!(from_full, from_narrowed);
    }
}
