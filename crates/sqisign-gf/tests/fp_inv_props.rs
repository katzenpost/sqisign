//! Property tests for `fp_inv`. The directive limits the property
//! surface to a single named cross-validation: the multiplicative
//! identity `fp_mul(fp_inv(a), a) ==_field 1` for canonical-nonzero
//! `a`. This is the value-level law that defines the modular inverse;
//! it cross-validates the just-landed `fp_inv` against the equally
//! just-landed `fp_is_equal` predicate and the previously-landed
//! `fp_mul` and `fp_set_one`, without relying on a re-implementation
//! of Fermat in the test crate.
//!
//! Soundness:
//! - The identity is stated in the value domain (`==_field`), not on
//!   raw limbs, so the redundant Montgomery representation is fine:
//!   `fp_is_equal` `redc`s both sides to canonical before comparing.
//! - `a` is constrained to be canonical-nonzero: arbitrary `[u64; 5]`
//!   limb patterns are *not* a valid Montgomery representative domain
//!   for the inversion identity (limb 4 above the redundant `[0, 2p)`
//!   range yields a deterministic-but-meaningless inverse output, the
//!   shape the differential vector battery pins). The test
//!   constructs canonical inputs by routing arbitrary bytes through
//!   [`fp_decode_reduce`] (the public path that produces a
//!   well-formed Montgomery representative of an integer mod `p`)
//!   and then filters out the field zero (on which Fermat
//!   squares-and-multiplies down to zero, breaking the identity).

use proptest::prelude::*;
use sqisign_gf::{
    fp_decode_reduce, fp_inv, fp_is_equal, fp_is_zero, fp_mul, fp_set_one, Fp, NWORDS_FIELD,
};

fn mul(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_mul(&mut c, a, b);
    c
}

fn one() -> Fp {
    let mut o = [0u64; NWORDS_FIELD];
    fp_set_one(&mut o);
    o
}

/// Build a canonical (in-range) Montgomery `Fp` from an arbitrary byte
/// seed by feeding it through `fp_decode_reduce`, the public path that
/// reduces an arbitrary-length byte string mod `p` into Montgomery
/// form. Every output is a well-formed representative of an integer in
/// `[0, p)`, the domain on which the inversion identity is sound.
fn canonical_from(seed: &[u8]) -> Fp {
    let mut a = [0u64; NWORDS_FIELD];
    fp_decode_reduce(&mut a, seed);
    a
}

proptest! {
    // fp_mul(fp_inv(a), a) ==_field 1, for any canonical-nonzero a.
    // The multiplicative-inverse identity, the value-level law that
    // defines modinv. Sound on the redundant Montgomery form because
    // the comparison uses fp_is_equal (which redc's both sides to
    // canonical before per-limb compare), not raw-limb equality. The
    // zero filter is necessary: on the field zero, Fermat returns
    // zero (the chain squares-and-multiplies zero down to zero), and
    // 0 * 0 = 0 != 1.
    #[test]
    fn inv_times_a_is_one(
        seed in proptest::collection::vec(any::<u8>(), 0..=128),
    ) {
        let a = canonical_from(&seed);
        prop_assume!(fp_is_zero(&a) == 0);
        let mut inv_a = a;
        fp_inv(&mut inv_a);
        let product = mul(&inv_a, &a);
        prop_assert_eq!(fp_is_equal(&product, &one()), 0xFFFF_FFFF);
    }
}
