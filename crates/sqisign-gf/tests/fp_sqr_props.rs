//! Property tests for `fp_sqr`.
//!
//! `fp_sqr` operates on the same **redundant, non-canonical** radix-2^51
//! form as `fp_add`/`fp_sub`/`fp_neg`/`fp_mul`: a residue class has many
//! limb encodings and `modsqr` reduces only to "less than 2p" (in the
//! Montgomery domain), leaving limb 4 *fully unmasked* (the reference's
//! final write is `c[4] = (spint)t`, a 64-bit truncation with **no**
//! `& mask`, identical to `modmul`). Raw-limb equality is therefore *not*
//! a sound notion of field equality, and the reference's own equality
//! (`modcmp`) is not ported yet. Only the properties below are sound on
//! the redundant representation; each was cross-checked bit-exactly
//! against the full 1012-vector C-derived battery before being committed
//! (the canonical correctness check remains the differential test in
//! `fp_sqr_vectors.rs`):
//!
//! 1. **`fp_sqr(a) == fp_mul(a, a)` bit-exact, for arbitrary limb
//!    inputs.** `modsqr` is the squaring specialisation of `modmul`: per
//!    column `k`, `modmul` sums every `a[i] * b[j]` with `i + j == k`
//!    into the 128-bit accumulator `t`, while `modsqr` (with `b = a`)
//!    instead builds a per-column accumulator `tot` from the off-diagonal
//!    `a[i] * a[j]` (i < j) products, doubles it (`tot *= 2;`) to account
//!    for the symmetric `a[j] * a[i]`, adds any diagonal `a[i] * a[i]`
//!    (when present in the column) un-doubled, and folds `t += tot;`.
//!    Algebraically the column sum is the same value as `modmul(a, a)`'s
//!    column sum; the Montgomery reduction folds (`v_{k-4} * p4`) and the
//!    unmasked limb-4 write are identical. Because `wrapping_add` and the
//!    underlying `__uint128_t` addition are associative and commutative
//!    modulo `2^128`, the running `t` is bit-equal at every masking point
//!    under the two orderings, so the recorded `v0..v4` and `c[0..=4]`
//!    are bit-identical. Verified empirically: across all 1012 committed
//!    vectors, `fp_sqr(a)` is raw-limb equal to `fp_mul(a, a)` for every
//!    record (matches: 1012, divergences: 0). The strongest sound
//!    raw-limb law for squaring on the redundant form, and the only one
//!    that exercises the squaring optimisation specifically.
//!
//! 2. **Structural carry-propagation invariant for limbs 0..=3, for
//!    arbitrary inputs.** The intermediate column writes
//!    `c[0..=3] = (t as u64) & MASK51` apply the per-limb mask, so every
//!    output has `out[0..4] < 2^51`. Limb 4 is intentionally *not*
//!    asserted: the reference's final write is the unmasked truncation
//!    `c[4] = (spint)t` and the port faithfully does the same.
//!    Verified: 0 violations across the 1012 committed vectors.
//!
//! ## What was considered and *omitted* as unsound
//!
//! - **`fp_sqr(MONTGOMERY_ONE) == MONTGOMERY_ONE`** (the Montgomery
//!   identity element is a multiplicative fixed point). Sound in
//!   principle but raw-limb-unsound: like `fp_mul(MONTGOMERY_ONE,
//!   MONTGOMERY_ONE)`, the output is `R^-1 mod p` positionally rather
//!   than `1`, and value-level comparison needs `modcmp` or `redc`
//!   (neither ported). *Omitted* until `redc` lands.
//! - **Value-level laws on the canonical Montgomery domain.** The same
//!   reasoning that justifies `fp_mul_props` omitting positional value
//!   laws applies here: `fp_sqr(a)` is `a^2 * R^-1 mod p`, not the
//!   positional `a^2 mod p`, and recovering the latter needs the inverse
//!   Montgomery factor `R`. *Omitted*.
//! - **`fp_sqr(0) == 0` (raw-limb).** Sound as a value, and indeed
//!   `modsqr` accumulates all-zero columns into an all-zero `t`, so the
//!   output is the bit-exact all-zero limb vector. *Omitted as
//!   redundant*: the differential battery's first record is the all-zero
//!   edge pattern (the same first edge `fp_neg`'s vectors use), so
//!   asserting it again at the property level would double-record the
//!   same bit. The `fp_sqr == fp_mul(a, a)` law also subsumes it via the
//!   `fp_mul(0, 0)` case.

use proptest::prelude::*;
use sqisign_gf::{fp_mul, fp_sqr, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;

fn sqr(a: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_sqr(&mut c, a);
    c
}

fn mul(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_mul(&mut c, a, b);
    c
}

proptest! {
    // (1) fp_sqr(a) bit-exact equals fp_mul(a, a) for arbitrary (possibly
    // non-canonical) limb inputs. Sound: the per-column set of partial
    // products is the same multiset in both cases (modsqr just bundles
    // the symmetric pair (i, j)/(j, i) into a single doubled term), and
    // u128 wrapping_add is associative and commutative, so t is bit-equal
    // at every masking point. Verified empirically against the full
    // 1012-vector battery (1012/1012 matches, 0 divergences) before
    // pinning.
    #[test]
    fn sqr_equals_mul_with_self_bit_exact(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(sqr(&a), mul(&a, &a));
    }

    // (2) Structural carry-propagation invariant: intermediate column
    // writes c[0..=3] are masked below 2^51; limb 4 is left unmasked by
    // design (the reference's final c[4] = (spint)t is a full 64-bit
    // truncation, no & mask) and is deliberately not constrained.
    #[test]
    fn limbs_0_3_below_radix(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        let c = sqr(&a);
        for (k, &limb) in c.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }
}
