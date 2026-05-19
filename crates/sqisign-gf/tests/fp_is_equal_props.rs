//! Property tests for `fp_is_equal`.
//!
//! `fp_is_equal` is the **binary predicate** counterpart of `fp_is_zero`:
//! it takes two `fp_t` operands in the redundant radix-2^51
//! representation and returns a `uint32_t` mask: `0xFFFFFFFF` for "the
//! two operands represent the same field element", `0` for "they do
//! not". The mask shape is load-bearing in the same way as for
//! `fp_is_zero`: the rest of the codebase ANDs it with field limbs (via
//! `fp_select(d, _, _, ctl)`), so a port that returned `0x1` rather
//! than `0xFFFFFFFF` for the positive outcome would silently zero out
//! the field arithmetic that consumes it.
//!
//! Only the properties below are sound on the redundant representation;
//! each was cross-checked against the full 1144-vector C-derived
//! battery before being committed (the canonical correctness check
//! remains the differential test in `fp_is_equal_vectors.rs`):
//!
//! 1. **Reflexivity, `fp_is_equal(a, a) == 0xFFFFFFFF`, for arbitrary
//!    `a`.** `modcmp` runs `redc` on each operand independently; `redc`
//!    is deterministic with no global state, so `redc(a) == redc(a)`
//!    bit-for-bit, every per-limb XOR is `0`, every per-limb
//!    `(x - 1) >> 51 & 1` is `1`, the AND-fold is `1`, and the
//!    `-(uint32_t)` wrapper yields the all-ones mask. Sound on
//!    *arbitrary* limb inputs (including non-canonical garbage),
//!    because the per-operand `redc` canonicalises before the
//!    comparison.
//!
//! 2. **Montgomery one is not equal to canonical zero,
//!    `fp_is_equal(fp_set_one_out, [0; 5]) == 0`.** The Montgomery
//!    representative of `1` is nonzero modulo `p` (it is `R mod p`,
//!    `R = 2^255 mod p`), so its canonical reduction is `[1, 0, 0, 0,
//!    0]`; the canonical zero reduces to `[0, 0, 0, 0, 0]`; the
//!    per-limb XOR at limb 0 is `1`, the per-limb `(1 - 1) >> 51 & 1`
//!    is `0`, the AND-fold is `0`, the wrapper yields `0`. The
//!    natural negative-case witness, mirroring the analogous fixed
//!    case in `fp_is_zero_props.rs`.
//!
//! 3. **Returned mask is always `0` or `0xFFFFFFFF`, for arbitrary
//!    limb inputs.** This is the C `-(uint32_t)int01` invariant: the
//!    inner `modcmp` returns `{0, 1}` (each per-limb bit it ANDs into
//!    `eq` is `{0, 1}` by the zero-detect trick, and the initial
//!    `eq = 1` is `{0, 1}`, so the AND-fold stays `{0, 1}`) and the
//!    wrapper's unary minus on `uint32_t` widens to `{0, 0xFFFFFFFF}`
//!    with no intermediate bit patterns possible. A port that forgot
//!    the negation (returning `{0, 1}`), shifted in the wrong
//!    direction, or used a signed cast would diverge here; verified
//!    empirically across the full 1144-vector battery (all outputs are
//!    `0` or `0xFFFFFFFF`).
//!
//! 4. **Symmetry, `fp_is_equal(a, b) == fp_is_equal(b, a)`, bit-exact,
//!    for arbitrary inputs.** `modcmp` runs `redc` on each operand
//!    independently, then per-limb XORs the canonical limbs; XOR is
//!    symmetric (`c[i] ^ d[i] == d[i] ^ c[i]`), the zero-detect trick
//!    applies the same shift and mask irrespective of operand order,
//!    and the AND-fold is commutative. So the returned mask is
//!    invariant under operand swap, bit-for-bit, not merely on the
//!    value level.
//!
//! ## What was considered and *omitted* as unsound or redundant
//!
//! - **Transitivity, `fp_is_equal(a, b) == 0xFFFFFFFF &&
//!   fp_is_equal(b, c) == 0xFFFFFFFF => fp_is_equal(a, c) ==
//!   0xFFFFFFFF`.** Sound as a value-level law (field equality is
//!   transitive), and indeed verifiable bit-exact since `redc` is a
//!   function. Omitted because the proptest seed for three independent
//!   fp_t's almost never has any pair equal: the property would
//!   trivially hold by both implications being false. A constructive
//!   proptest (generate two redundant representatives of the same
//!   field element, then a third) needs more reduction-algebra
//!   machinery than the property suite carries today; the differential
//!   battery already pins the canonical-zero / p-encoded-zero
//!   diagonal-and-off-diagonal pair, which exercises the
//!   non-canonical-equality branch.
//! - **`fp_is_equal(a, b) == fp_is_zero(fp_sub(a, b))`.** Sound as a
//!   value-level law (a - b == 0 iff a == b), and bit-exact on the
//!   canonical domain. But raw-limb equivalence requires that
//!   `fp_sub`'s redundant `[0, 2p)` representative reduces to the
//!   canonical zero under `modis0`'s `redc`; while this is true
//!   (modis0 calls redc first), the law is then subsumed by the
//!   differential test for both fp_is_equal and fp_sub. *Omitted as
//!   redundant.*
//! - **Value-level laws on arbitrary redundant representatives of the
//!   same field element.** E.g. "for any two representatives of the
//!   same field element, `fp_is_equal` returns `0xFFFFFFFF`". Sound
//!   in principle, but constructing redundant-representative pairs
//!   from a proptest seed requires a model of the reduction algebra,
//!   and the differential battery already pins one such pair (the
//!   canonical zero and the radix-2^51 encoding of `p`). *Omitted as
//!   redundant.*

use proptest::prelude::*;
use sqisign_gf::{fp_is_equal, fp_set_one, Fp, NWORDS_FIELD};

#[test]
fn montgomery_one_not_equal_to_canonical_zero() {
    // (b) The natural fixed-point negative: fp_is_equal of the
    // Montgomery representative of 1 and the canonical zero is the
    // all-zero mask. We obtain the Montgomery one via the public
    // fp_set_one (which writes the exposed const ONE directly),
    // avoiding any dependency on MONTGOMERY_ONE being public.
    let mut one: Fp = [0u64; NWORDS_FIELD];
    fp_set_one(&mut one);
    let zero: Fp = [0u64; NWORDS_FIELD];
    assert_eq!(fp_is_equal(&one, &zero), 0);
}

proptest! {
    // (a) Reflexivity: fp_is_equal(a, a) is the all-ones mask for
    // arbitrary (possibly non-canonical) limb inputs. Sound: redc is
    // deterministic with no global state, so redc(a) == redc(a)
    // bit-for-bit, every per-limb XOR is 0, the AND-fold is 1, the
    // wrapper yields the all-ones mask.
    #[test]
    fn reflexive_is_all_ones(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(fp_is_equal(&a, &a), 0xFFFF_FFFF);
    }

    // (c) The returned mask is always 0 or 0xFFFFFFFF for arbitrary
    // (possibly non-canonical) limb inputs. The C `-(uint32_t)int01`
    // invariant: any other bit pattern would indicate the port forgot
    // the negation, used a different shift, or mishandled the cast
    // chain. Verified empirically across the full 1144-vector battery.
    #[test]
    fn mask_is_zero_or_all_ones(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        let m = fp_is_equal(&a, &b);
        prop_assert!(
            m == 0 || m == 0xFFFF_FFFF,
            "fp_is_equal returned non-mask value {m:#x}"
        );
    }

    // (d) Symmetry: fp_is_equal(a, b) == fp_is_equal(b, a) bit-exact
    // for arbitrary inputs. Sound: XOR is symmetric (c[i] ^ d[i] ==
    // d[i] ^ c[i]), the zero-detect trick applies the same shift and
    // mask irrespective of operand order, and the AND-fold is
    // commutative.
    #[test]
    fn symmetric_bit_exact(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(fp_is_equal(&a, &b), fp_is_equal(&b, &a));
    }
}
