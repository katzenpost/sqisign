//! Property tests for `fp_sub`.
//!
//! `fp_sub` operates on the same **redundant, non-canonical** radix-2^51
//! form as `fp_add`: a residue class has many limb encodings and `modsub`
//! reduces only to "less than 2p", leaving limb 4 unmasked. Raw-limb
//! equality is therefore *not* a sound notion of field equality, and the
//! reference's own equality (`modcmp`) is not ported yet. Subtraction is
//! also not commutative, so `fp_add`'s strongest raw-limb law (bit-exact
//! commutativity) does not transfer. Only the properties below are sound;
//! each was cross-checked bit-exactly against the full 1144-vector
//! C-derived battery before being committed (the canonical correctness
//! check remains the differential test in `fp_sub_vectors.rs`):
//!
//! 1. **`fp_sub(a, a)` is the canonical all-zero representative, bit-exact,
//!    for arbitrary limb inputs.** `modsub` forms `n[i] = a[i] - a[i]`,
//!    which is exactly `0` in every limb (unsigned wraparound cancels);
//!    `prop` then sees an all-zero value, returns a zero carry mask, and
//!    no `2p` correction fires, leaving all five limbs zero. Verified
//!    empirically: across all 1144 committed vectors, every `a == b`
//!    input yields all-zero output (12 cases), and conversely an all-zero
//!    output occurs *only* when `a == b` (no other vector outputs zero).
//!    This is the sharpest sound raw-limb law for subtraction and pins
//!    the zero representative.
//!
//! 2. **Structural carry-propagation invariant, for arbitrary inputs.**
//!    The final `prop` masks limbs 0..=3 with `(1<<51)-1`, so every output
//!    has `out[0..4] < 2^51`. Limb 4 is intentionally *not* asserted: the
//!    reference leaves it unmasked and the port faithfully does too.
//!    Verified: 0 violations across the 1144 committed vectors.
//!
//! ## Why no value-level (mod-p) property is asserted
//!
//! `fp_add_props` asserts a canonical-domain value law
//! (`value_mod_p` of the output equals the sum mod `p`), reading the
//! output positionally as `sum limb[i] * 2^(51*i) mod p`. That positional
//! reading is sound for `modadd` *only* because `canonical + canonical`
//! has all-positive limbs and stays in `[0, 2p)`, so `modadd`'s single
//! conditional `-p`-shaped correction leaves a representative whose plain
//! positional value already equals the field element. `modsub` breaks
//! that: `canonical - canonical` can be negative, and `modsub`'s `+2p`
//! correction (the redundant `+2 at limb 0`, `-2*p4 at limb 4` pair)
//! produces a representative whose plain positional value is *not*
//! congruent to the field value: recovering it requires the reference's
//! own `redc`/`modfsb` canonicalization, which is not ported yet.
//! Verified directly: for `a = 0`, `b` a canonical encoding of `2`,
//! `modsub` leaves limb 0 = `0xffff_ffff_ffff_fffe` whose positional
//! value mod `p` is `2^64 - 2`, not the field value `p - 2`. Asserting an
//! fp_add-style value law here would therefore be *unsound*, so it is
//! deliberately omitted rather than weakened or silently corrected; the
//! differential test against the C oracle remains the value-correctness
//! authority until `redc` lands.

use proptest::prelude::*;
use sqisign_gf::{fp_sub, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;

fn sub(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_sub(&mut c, a, b);
    c
}

proptest! {
    // (1) fp_sub(a, a) is the canonical all-zero representative, bit-exact,
    // arbitrary (possibly non-canonical) limb inputs. Sound: modsub's
    // limbwise difference is exactly zero in every limb and the prop tail
    // applies no correction. Verified against the full battery.
    #[test]
    fn self_difference_is_zero(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(sub(&a, &a), [0u64; NWORDS_FIELD]);
    }

    // (2) Structural carry-propagation invariant: limbs 0..=3 are masked
    // below 2^51 by the final prop; limb 4 is left unmasked by design and
    // is deliberately not constrained.
    #[test]
    fn limbs_0_3_below_radix(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        let c = sub(&a, &b);
        for (k, &limb) in c.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }
}
