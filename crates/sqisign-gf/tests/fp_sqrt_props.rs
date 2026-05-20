//! Property tests for `fp_sqrt`. The directive limits the property
//! surface to a single named cross-validation: `fp_sqr(fp_sqrt(a))
//! ==_field a` for any `a` that is a quadratic residue (as determined
//! by `fp_is_square`). This is the value-level law that defines the
//! square root; it cross-validates the just-landed `fp_sqrt` against
//! the equally just-landed `fp_is_square` and `fp_is_equal`, and the
//! previously-landed `fp_sqr`, without relying on a re-implementation
//! of Tonelli-Shanks or the progenitor chain in the test crate.
//!
//! Soundness:
//! - The identity is stated in the value domain (`==_field`), not on
//!   raw limbs, so the redundant Montgomery representation is fine:
//!   `fp_is_equal` `redc`s both sides to canonical before comparing.
//! - `a` is constrained to a quadratic residue by `fp_is_square(&a) ==
//!   0xFFFF_FFFF`; on a non-residue, the reference (and the port)
//!   return meaningless garbage and the squared-back value will not
//!   equal `a`, so the filter is necessary for soundness. Per the
//!   reference convention, `fp_is_square` returns the positive mask on
//!   the field zero too, and the identity holds there trivially
//!   (`sqrt(0) * sqrt(0) == 0 * 0 == 0`).
//! - `a` is additionally constrained to be canonical (a well-formed
//!   Montgomery representative of an integer in `[0, p)`): arbitrary
//!   `[u64; 5]` limb patterns whose limb 4 is outside the redundant
//!   `[0, 2p)` range yield deterministic-but-meaningless `fp_sqrt`
//!   output (the shape the differential vectors pin) and break the
//!   round-trip identity. The test constructs canonical inputs by
//!   routing arbitrary bytes through [`fp_decode_reduce`] (the
//!   public path that produces a well-formed Montgomery
//!   representative of an integer mod `p`).

use proptest::prelude::*;
use sqisign_gf::{fp_decode_reduce, fp_is_equal, fp_is_square, fp_sqr, fp_sqrt, Fp, NWORDS_FIELD};

fn sqr(a: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_sqr(&mut c, a);
    c
}

/// Build a canonical (in-range) Montgomery `Fp` from an arbitrary
/// byte seed by feeding it through `fp_decode_reduce`, the public
/// path that reduces an arbitrary-length byte string mod `p` into
/// Montgomery form. Every output is a well-formed representative of
/// an integer in `[0, p)`, the domain on which the `sqrt`/`sqr`
/// round-trip identity is sound.
fn canonical_from(seed: &[u8]) -> Fp {
    let mut a = [0u64; NWORDS_FIELD];
    fp_decode_reduce(&mut a, seed);
    a
}

proptest! {
    // fp_sqr(fp_sqrt(a)) ==_field a, for any a that is a quadratic
    // residue (or the field zero, per the reference's convention that
    // zero is a square). The square-back identity, the value-level
    // law that defines modsqrt for p == 3 mod 4. Sound on the
    // redundant Montgomery form because the comparison uses
    // fp_is_equal (which redc's both sides to canonical before
    // per-limb compare), not raw-limb equality. The QR filter is
    // necessary: on a non-residue the reference returns garbage and
    // the round-trip identity does not hold; the filter cuts the
    // input space to the domain on which the contract is meaningful.
    #[test]
    fn sqrt_squared_recovers_value(
        seed in proptest::collection::vec(any::<u8>(), 0..=128),
    ) {
        let a = canonical_from(&seed);
        prop_assume!(fp_is_square(&a) == 0xFFFF_FFFF);
        let mut root = a;
        fp_sqrt(&mut root);
        let back = sqr(&root);
        prop_assert_eq!(fp_is_equal(&back, &a), 0xFFFF_FFFF);
    }
}
