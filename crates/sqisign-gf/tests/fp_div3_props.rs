//! Property tests for `fp_div3`.
//!
//! `fp_div3` operates on the same **redundant, non-canonical** radix-2^51
//! form as `fp_add`/`fp_sub`/`fp_neg`/`fp_mul`/`fp_sqr`/`fp_half`: it is
//! the one-liner `modmul(THREE_INV, *a, *out);` (see
//! `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:658..662`), so a
//! residue class has many limb encodings and the underlying `modmul`
//! reduces only to "less than 2p" (in the Montgomery domain), leaving
//! limb 4 *fully unmasked* (the reference's final write is
//! `c[4] = (spint)t`, a 64-bit truncation with **no** `& mask`). Raw-limb
//! equality is therefore *not* a sound notion of field equality on
//! `fp_div3`'s output. Only the properties below are sound; each was
//! cross-checked bit-exactly against the full 1012-vector C-derived
//! battery before being committed (the canonical correctness check
//! remains the differential test in `fp_div3_vectors.rs`):
//!
//! 1. **`fp_div3(a) == fp_mul(&THREE_INV, a)` bit-exact, for arbitrary
//!    limb inputs.** Sound by construction: the port of `fp_div3` *is*
//!    literally the one call `modmul(&THREE_INV, a, out)`, and `fp_mul`
//!    is the public wrapper for the same `modmul`. The two execute the
//!    bit-identical column-by-column accumulator trace, so the recorded
//!    `c[0..=4]` are bit-identical regardless of input. The property is
//!    therefore a *tautology* at the source level (both call sites
//!    invoke the same `modmul` with the same operands); its value is as
//!    a **sanity oracle for the `THREE_INV` constant**: any
//!    mis-transcribed limb in the internal `THREE_INV` would diverge
//!    from the local duplicate below on the very first input. (The
//!    duplicate is here because the proptest crate cannot access the
//!    private internal `THREE_INV` in `lib.rs`; the cross-oracle keeps
//!    the two byte-for-byte in sync.) Verified empirically across all
//!    1012 committed vectors (matches: 1012, divergences: 0).
//!
//! 2. **`3 * fp_div3(a) ==_field a`, for arbitrary limb inputs.** The
//!    value-level identity: dividing by three and then tripling is the
//!    identity on the field. Sound on the *redundant* form because the
//!    equality is expressed via [`fp_is_equal`] (which `redc`s both
//!    operands to their canonical representatives before comparing per
//!    limb), not raw-limb equality of the two redundant representatives.
//!    In symbols: `fp_div3(a)` represents `A / 3 mod p` in the Montgomery
//!    domain, and `fp_add(fp_add(d, d), d)` where `d = fp_div3(a)`
//!    represents `(A / 3) + (A / 3) + (A / 3) == A mod p`, so
//!    [`fp_is_equal`] returns the all-ones mask. This is the strongest
//!    sound value-level law for the division-by-three operation and the
//!    one that exercises the `THREE_INV` constant's defining property
//!    (its Montgomery rep is `3^-1 R`); a `THREE_INV` perturbed by
//!    anything other than a `p`-multiple would break the tripling-back
//!    identity. Verified empirically across all 1012 committed vectors
//!    (matches: 1012, divergences: 0).
//!
//!    The tripling-back oracle subsumes the simpler raw-limb-unsound
//!    "`fp_div3(MONTGOMERY_ONE)` is the Montgomery rep of `1/3`": both
//!    rely on the same `THREE_INV * R` identity, and (2) checks it
//!    value-wise on every input rather than positionally on one.
//!
//! 3. **Structural carry-propagation invariant for limbs 0..=3, for
//!    arbitrary inputs.** The intermediate column writes in `modmul`
//!    apply `(t as u64) & MASK51` per limb, so every output of
//!    `fp_div3` (which is `modmul`) has `out[0..4] < 2^51`. Limb 4 is
//!    intentionally *not* asserted: the reference's final write is the
//!    unmasked truncation `c[4] = (spint)t` and the port faithfully
//!    does the same. Verified: 0 violations across the 1012 committed
//!    vectors.
//!
//! ## What was considered and *omitted* as unsound
//!
//! - **`fp_div3(MONTGOMERY_ONE) == Montgomery rep of 1/3`** as raw-limb
//!   equality. `1/3 mod p` is a single fixed residue class, and
//!   `fp_div3` should produce its Montgomery representative. *Omitted*
//!   in favour of property (2)'s value-level `fp_is_equal` formulation,
//!   which is sound on the redundant form and exercises every input
//!   rather than one. The raw-limb pinning adds no information beyond
//!   what the differential vector battery already provides on the
//!   canonical edge cases.
//! - **Value-level laws on the canonical (non-Montgomery) domain.** The
//!   same reasoning that justifies `fp_mul_props` / `fp_sqr_props` /
//!   `fp_half_props` omitting positional value laws applies here:
//!   `fp_div3(a)` is `(a/3) * R^-1 mod p` relative to the *positional*
//!   a, not relative to the Montgomery-domain `A` it represents. Stating
//!   the law in Montgomery-domain terms ("tripling-back is the identity
//!   mod p") is the soundest formulation, which is what (2) does.

use proptest::prelude::*;
use sqisign_gf::{fp_add, fp_div3, fp_is_equal, fp_mul, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;

/// Montgomery representative of `3^-1 mod p`, transcribed verbatim from
/// `vendor/the-sqisign/src/gf/ref/lvl1/fp_p5248_64.c:538..542`. Kept in
/// sync with the internal `THREE_INV` constant in `sqisign_gf::lib`; the
/// test crate cannot access private items so the constant is duplicated
/// here. Any drift between the two is exactly what property (1) catches
/// on the first input.
const THREE_INV: Fp = [
    0x0005_5555_5555_555d,
    0x0002_aaaa_aaaa_aaaa,
    0x0005_5555_5555_5555,
    0x0002_aaaa_aaaa_aaaa,
    0x0000_4555_5555_5555,
];

fn div3(a: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_div3(&mut c, a);
    c
}

fn mul(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_mul(&mut c, a, b);
    c
}

fn add(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_add(&mut c, a, b);
    c
}

proptest! {
    // (1) fp_div3(a) bit-exact equals fp_mul(&THREE_INV, a). The two
    // call sites invoke the same modmul with the same operands, so the
    // property is a source-level tautology that doubles as a sanity
    // oracle for the THREE_INV constant: any mis-transcribed limb in
    // the internal THREE_INV would diverge from the local duplicate on
    // the first input. Verified empirically against the full 1012-vector
    // battery (1012/1012 matches, 0 divergences) before pinning.
    #[test]
    fn div3_equals_mul_with_three_inv_bit_exact(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(div3(&a), mul(&THREE_INV, &a));
    }

    // (2) 3 * fp_div3(a) ==_field a. The value-level identity: dividing
    // by three and tripling is the identity. Sound on the redundant
    // form because the comparison uses fp_is_equal (which redc's both
    // sides to canonical before per-limb compare), not raw-limb
    // equality. This is the strongest sound value-level law on
    // division-by-three, and it exercises the THREE_INV constant's
    // defining property (its Montgomery rep is 3^-1 R; perturbing it
    // by anything other than a p-multiple breaks tripling-back).
    // Verified empirically against the full 1012-vector battery
    // (1012/1012 matches, 0 divergences) before pinning.
    #[test]
    fn tripling_div3_recovers_value(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        let d = div3(&a);
        let tripled = add(&add(&d, &d), &d);
        prop_assert_eq!(fp_is_equal(&tripled, &a), 0xFFFF_FFFF);
    }

    // (3) Structural carry-propagation invariant: intermediate column
    // writes in modmul mask c[0..=3] below 2^51; limb 4 is left
    // unmasked by design (the reference's final c[4] = (spint)t is a
    // full 64-bit truncation, no & mask) and is deliberately not
    // constrained, exactly as fp_mul / fp_sqr / fp_half leave it.
    #[test]
    fn limbs_0_3_below_radix(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        let c = div3(&a);
        for (k, &limb) in c.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }
}
