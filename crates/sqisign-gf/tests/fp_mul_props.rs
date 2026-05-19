//! Property tests for `fp_mul`.
//!
//! `fp_mul` operates on the same **redundant, non-canonical** radix-2^51
//! form as `fp_add`/`fp_sub`/`fp_neg`: a residue class has many limb
//! encodings and `modmul` reduces only to "less than 2p" (in the
//! Montgomery domain), leaving limb 4 *fully unmasked* (the reference's
//! final write is `c[4] = (spint)t`, a 64-bit truncation with **no**
//! `& mask`). Raw-limb equality is therefore *not* a sound notion of
//! field equality, and the reference's own equality (`modcmp`) is not
//! ported yet. Only the properties below are sound on the redundant
//! representation; each was cross-checked bit-exactly against the full
//! 1144-vector C-derived battery before being committed (the canonical
//! correctness check remains the differential test in
//! `fp_mul_vectors.rs`):
//!
//! 1. **Commutativity, bit-exact, for arbitrary limb inputs.** Although
//!    `modmul`'s source is *not* visibly symmetric in `a` and `b` (it
//!    organises by `a[i] * b[j]` diagonals in a specific order, and the
//!    `v_k * p4` Montgomery folds depend on the running accumulator at
//!    each masking point), the per-column *set* of partial products
//!    `{ a[i] * b[j] : i + j == k }` is operand-symmetric (swap is the
//!    same set with `(i, j)` -> `(j, i)`), and the 128-bit accumulator
//!    additions are associative and commutative; the running `t` value
//!    at every masking point is therefore the same under operand swap,
//!    so `v0..v4` and `c[0..=4]` come out bit-identical. Verified
//!    empirically: across all 1144 committed vectors, `fp_mul(a, b)` is
//!    raw-limb equal to `fp_mul(b, a)` for every record. The strongest
//!    sound raw-limb law for multiplication.
//!
//! 2. **Structural carry-propagation invariant for limbs 0..=3, for
//!    arbitrary inputs.** The intermediate column writes
//!    `c[0..=3] = (t as u64) & MASK51` apply the per-limb mask, so
//!    every output has `out[0..4] < 2^51`. Limb 4 is intentionally *not*
//!    asserted: the reference's final write is the unmasked truncation
//!    `c[4] = (spint)t` and the port faithfully does the same.
//!    Verified: 0 violations across the 1144 committed vectors.
//!
//! ## What was considered and *omitted* as unsound
//!
//! - **`fp_mul(a, MONTGOMERY_ONE) == redc(a)`** (and its mirror
//!   `fp_mul(MONTGOMERY_ONE, a) == redc(a)`). Sound in principle: this
//!   is the defining identity of the Montgomery domain. *Omitted* here
//!   because `redc` is not yet ported, so there is no canonical
//!   oracle to compare against from within the Rust port. It will be
//!   added when `redc` lands; for now the differential vectors are the
//!   value-correctness authority.
//! - **Value-level associativity `fp_mul(fp_mul(a, b), c) ==
//!   fp_mul(a, fp_mul(b, c))`.** Sound at the field-value level on the
//!   canonical Montgomery domain, but raw-limb equality is *not* the
//!   correct comparator: distinct redundant representatives of the same
//!   residue class are bit-distinct, and any value-level comparison
//!   needs `modcmp` or `redc`/`modfsb`, neither of which is ported.
//!   *Omitted* rather than weakened.
//! - **Value-level distributivity over `fp_add`/`fp_sub`.** Same reason:
//!   needs `modcmp` or canonicalisation to compare at the value level.
//!   *Omitted*.
//! - **`fp_mul(a, ZERO) == ZERO` (raw-limb).** Sound as a value, but
//!   raw-limb-unsound: `modmul` accumulates all-zero columns into an
//!   all-zero `t`, so the result *is* the bit-exact all-zero limb
//!   vector; the verified empirical behaviour matches. *Omitted as
//!   redundant*: the differential battery already pins this for every
//!   `(a, 0)` and `(0, b)` record in the edge battery (column-zero
//!   patterns), and asserting it again at the property level would
//!   double-record the same bit.
//!
//! ## Why no fp_add-style canonical-domain value law is asserted
//!
//! `fp_add_props` asserts `value_mod_p(fp_add(a, b)) == (a + b) mod p`
//! on the canonical sub-domain, reading the output positionally as
//! `sum limb[i] * 2^(51*i) mod p`. That positional reading is sound for
//! `modadd` because addition stays in the additive group regardless of
//! Montgomery factor. `modmul` breaks the analogue: even on the
//! canonical Montgomery sub-domain, the output is `a * b * R^-1 mod p`,
//! not `a * b mod p`; recovering the latter from the former needs the
//! inverse Montgomery factor `R`, which in turn needs `redc`/`nres`,
//! neither ported. Asserting a positional value law without that
//! correction would be unsound (the witness is any non-trivial product:
//! e.g. `fp_mul(MONTGOMERY_ONE, MONTGOMERY_ONE)` is *not* the positional
//! `1`, it is the Montgomery `R^-1 mod p`). The differential test
//! against the C oracle therefore remains the value-correctness
//! authority until `redc` lands.

use proptest::prelude::*;
use sqisign_gf::{fp_mul, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;

fn mul(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_mul(&mut c, a, b);
    c
}

proptest! {
    // (1) Bit-exact commutativity, arbitrary (possibly non-canonical)
    // limb inputs. Sound: the per-column set of partial products is
    // operand-symmetric and u128 accumulator addition is associative and
    // commutative, so t is bit-equal at every masking point under operand
    // swap. Verified against the full 1144-vector battery.
    #[test]
    fn commutative_bit_exact(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(mul(&a, &b), mul(&b, &a));
    }

    // (2) Structural carry-propagation invariant: intermediate column
    // writes c[0..=3] are masked below 2^51; limb 4 is left unmasked by
    // design (the reference's final c[4] = (spint)t is a full 64-bit
    // truncation, no & mask) and is deliberately not constrained.
    #[test]
    fn limbs_0_3_below_radix(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        let c = mul(&a, &b);
        for (k, &limb) in c.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }
}
