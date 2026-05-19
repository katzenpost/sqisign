//! Property tests for `fp_is_zero`.
//!
//! `fp_is_zero` is the first **predicate** boundary in the gf battery. It
//! takes an `fp_t` in the redundant radix-2^51 representation and returns
//! a `uint32_t` mask: `0xFFFFFFFF` for "represents the field zero", `0`
//! for nonzero. The mask shape is load-bearing: the rest of the codebase
//! ANDs it with field limbs (`fp_select(d, _, _, ctl)` in particular),
//! so a port that returned `0x1` rather than `0xFFFFFFFF` for the
//! positive outcome would silently zero out the field arithmetic that
//! consumes it.
//!
//! Only the properties below are sound on the redundant representation;
//! each was cross-checked against the full 1012-vector C-derived
//! battery before being committed (the canonical correctness check
//! remains the differential test in `fp_is_zero_vectors.rs`):
//!
//! 1. **`fp_is_zero([0, 0, 0, 0, 0])` returns the all-ones mask.** The
//!    canonical zero is the easiest fixed point: `modis0` runs `redc`
//!    (which leaves all-zero limbs untouched modulo `p`), OR-folds to
//!    `d = 0`, and the bit-twiddle `((d - 1) >> 51) & 1` yields `1`
//!    (`0 - 1` wraps to all-ones, whose `>> 51` is `0x1FFF`, whose low
//!    bit is `1`); the `-(uint32_t)` wrapper turns the `1` into
//!    `0xFFFFFFFF`. Pinned as a fixed case, not a proptest, because it
//!    is a single point.
//! 2. **`fp_is_zero(MONTGOMERY_ONE)` returns `0`.** The Montgomery
//!    representative of `1` is nonzero modulo `p` (it is `R mod p`,
//!    `R = 2^255 mod p`), so its canonical reduction via `redc` is the
//!    nonzero limb pattern `[1, 0, 0, 0, 0]`; OR-fold gives `d = 1`,
//!    the bit-twiddle yields `0`, the wrapper yields `0`. This is the
//!    second pinned fixed point and the natural negative-case witness:
//!    if a port accidentally tested raw-limb equality to the canonical
//!    zero, this case would slip through (Montgomery one's limbs are
//!    not all zero), but the predicate over the *redc*'d value still
//!    correctly says "nonzero".
//! 3. **Returned mask is always `0` or `0xFFFFFFFF`, for arbitrary
//!    limb inputs.** This is the C `-(uint32_t)int01` invariant: the
//!    inner `modis0` returns `{0, 1}` and the wrapper's unary minus on
//!    `uint32_t` widens to `{0, 0xFFFFFFFF}` with no intermediate bit
//!    patterns possible. A port that forgot the negation (returning
//!    `{0, 1}`), shifted in the wrong direction, or used a signed cast
//!    would diverge here; verified empirically across the full
//!    1012-vector battery (all outputs are `0` or `0xFFFFFFFF`).
//!
//! ## What was considered and *omitted* as unsound
//!
//! - **`fp_is_zero(fp_neg(a)) == fp_is_zero(a)`** (the negation
//!   preserves zero-ness). Sound as a value-level law (since `-0 == 0`
//!   and `-a` is nonzero iff `a` is nonzero modulo `p`), and indeed
//!   `redc(fp_neg(a))` and `redc(a)` are congruent. But asserting it on
//!   *arbitrary* (non-canonical) inputs depends on the redundant form
//!   not introducing an artificial zero; this is the case here
//!   (`modis0` canonicalises first), but the law is then subsumed by
//!   the differential test (every recorded `(a, result)` pair is
//!   already pinned, and `fp_neg(a)` is a different `a`). *Omitted as
//!   redundant.*
//! - **`fp_is_zero(fp_sub(a, a)) == 0xFFFFFFFF`.** Same caveat: sound
//!   as a value-level law, but redundantly checked via the differential
//!   suite for `fp_sub` (which pins `fp_sub(a, a) == [0, 0, 0, 0, 0]`
//!   for canonical inputs) composed with the canonical-zero case here.
//!   *Omitted.*
//! - **Value-level laws on arbitrary redundant representatives of `0`.**
//!   E.g. "for any representative of `0 mod p`, `fp_is_zero` returns
//!   `0xFFFFFFFF`". Sound in principle, but constructing redundant
//!   zero representatives from a proptest seed requires a model of the
//!   reduction algebra, and the differential battery already pins two
//!   distinct redundant zeros (the canonical `[0, 0, 0, 0, 0]` and
//!   the radix-2^51 encoding of `p`); generating a wider family in the
//!   property suite would duplicate that pin without strengthening it.
//!   *Omitted as redundant.*

use proptest::prelude::*;
use sqisign_gf::{fp_is_zero, fp_set_one, Fp, NWORDS_FIELD};

#[test]
fn canonical_zero_returns_all_ones_mask() {
    // (1) The lone sound fixed-point positive: fp_is_zero of the
    // canonical all-zero limb pattern is the all-ones mask. Pinned as a
    // fixed case, not a proptest, because it is a single point.
    let zero: Fp = [0u64; NWORDS_FIELD];
    assert_eq!(fp_is_zero(&zero), 0xFFFF_FFFF);
}

#[test]
fn montgomery_one_returns_zero_mask() {
    // (2) The natural fixed-point negative: fp_is_zero of the
    // Montgomery representative of 1 is the all-zero mask. We obtain
    // the Montgomery one via the public fp_set_one (which writes the
    // exposed const ONE directly), avoiding any dependency on
    // MONTGOMERY_ONE being public.
    let mut one: Fp = [0u64; NWORDS_FIELD];
    fp_set_one(&mut one);
    assert_eq!(fp_is_zero(&one), 0);
}

proptest! {
    // (3) The returned mask is always 0 or 0xFFFFFFFF for arbitrary
    // (possibly non-canonical) limb inputs. The C `-(uint32_t)int01`
    // invariant: any other bit pattern would indicate the port forgot
    // the negation, used a different shift, or mishandled the cast
    // chain. Verified empirically across the full 1012-vector battery.
    #[test]
    fn mask_is_zero_or_all_ones(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        let m = fp_is_zero(&a);
        prop_assert!(
            m == 0 || m == 0xFFFF_FFFF,
            "fp_is_zero returned non-mask value {m:#x}"
        );
    }
}
