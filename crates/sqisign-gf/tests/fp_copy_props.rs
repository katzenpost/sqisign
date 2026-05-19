//! Property tests for `fp_copy`.
//!
//! `fp_copy` is the reference's `modcpy`: a plain five-limb assignment,
//! `out[i] = a[i]` for `i` in `0..5`. No `prop`, no `2p` correction, no
//! reduction. Unlike `fp_add`/`fp_sub`/`fp_neg` it is therefore not
//! constrained to the redundant-mod-`p` notion of equality: the output is
//! **bit-exact equal to the input** on every input, canonical or not.
//!
//! Two sound raw-limb properties:
//!
//! 1. **Identity.** For arbitrary five-limb `a`, `fp_copy(out, a)` leaves
//!    `out == a`, bit for bit. Verified against the full 1012-vector
//!    C-derived battery (`fp_copy_vectors.rs` pins this as a
//!    `count == total` assertion).
//! 2. **Chained idempotence.** Copying through any number of intermediate
//!    destinations leaves the value bit-exact unchanged: the second copy
//!    has the first copy's output as its input, which is by (1) the
//!    original `a`, so the second copy's output is also `a`. A chained
//!    copy is therefore observationally indistinguishable from a single
//!    copy. Asserted at chain length 3 with a deliberately distinct
//!    pre-fill in each destination to expose a no-op or aliased write.
//!
//! Pre-fills are non-trivial (a striped pattern) so a port that left some
//! limbs untouched would diverge visibly. `fp_copy` accepts non-canonical
//! limbs by construction, so the inputs are drawn from arbitrary `u64`
//! limb tuples, not the canonical subdomain.

use proptest::prelude::*;
use sqisign_gf::{fp_copy, Fp, NWORDS_FIELD};

fn copy(a: &Fp) -> Fp {
    let mut c: Fp = [0xa5a5_a5a5_a5a5_a5a5u64; NWORDS_FIELD];
    fp_copy(&mut c, a);
    c
}

proptest! {
    // (1) Identity: out == a, bit-exact, for arbitrary five-limb a.
    #[test]
    fn copies_exactly(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(copy(&a), a);
    }

    // (1') Destination prior contents are fully overwritten: a port that
    // wrote only some limbs would leave the fill byte visible. The fill
    // is itself fuzzed so no single sentinel is privileged.
    #[test]
    fn overwrites_destination(
        a in proptest::array::uniform5(any::<u64>()),
        fill in any::<u64>(),
    ) {
        let mut c: Fp = [fill; NWORDS_FIELD];
        fp_copy(&mut c, &a);
        prop_assert_eq!(c, a);
    }

    // (2) Chained copies are idempotent: the value at the end of a chain
    // of three copies, each through a freshly distinct pre-filled buffer,
    // equals the original input bit for bit.
    #[test]
    fn chained_copies_idempotent(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        let mut b: Fp = [0x1111_1111_1111_1111u64; NWORDS_FIELD];
        let mut c: Fp = [0x2222_2222_2222_2222u64; NWORDS_FIELD];
        let mut d: Fp = [0x3333_3333_3333_3333u64; NWORDS_FIELD];
        fp_copy(&mut b, &a);
        fp_copy(&mut c, &b);
        fp_copy(&mut d, &c);
        prop_assert_eq!(b, a);
        prop_assert_eq!(c, a);
        prop_assert_eq!(d, a);
    }
}
