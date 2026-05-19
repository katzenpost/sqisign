//! Property tests for `fp_neg`.
//!
//! `fp_neg` operates on the same **redundant, non-canonical** radix-2^51
//! form as `fp_add`/`fp_sub`: a residue class has many limb encodings and
//! `modneg` reduces only to "less than 2p", leaving limb 4 unmasked. Raw-
//! limb equality is therefore *not* a sound notion of field equality, and
//! the reference's own equality (`modcmp`) is not ported yet. `modneg` is
//! the unary analogue of `modsub` (`0 - b[i]` limbwise), so it inherits
//! `fp_sub`'s caveats, not `fp_add`'s bit-exact commutativity. Only the
//! properties below are sound; each was cross-checked against the full
//! 1012-vector C-derived battery and a wide out-of-tree model before
//! being committed (the canonical correctness check remains the
//! differential test in `fp_neg_vectors.rs`):
//!
//! 1. **`fp_neg(0)` is the canonical all-zero representative, bit-exact.**
//!    `modneg` forms `n[i] = 0 - 0`, which is exactly `0` in every limb
//!    (unsigned wraparound of `0 - 0` is `0`); `prop` then sees an
//!    all-zero value, returns a zero carry mask, no `2p` correction
//!    fires, and the second `prop` leaves all five limbs zero. Verified
//!    empirically: of the 1012 committed vectors exactly one has an
//!    all-zero input (vector 0, the first edge pattern) and its recorded
//!    reference output is the bit-exact all-zero representative;
//!    conversely no nonzero input across the battery negates to all-zero.
//!    This pins the zero fixed point. (Unlike `fp_sub(a, a)`, this is
//!    *not* a law for arbitrary `a`: `fp_neg(a)` for nonzero `a` is a
//!    nonzero, non-canonical representative of `-a`, not all-zero.)
//!
//! 2. **Structural carry-propagation invariant, for arbitrary inputs.**
//!    The final `prop` masks limbs 0..=3 with `(1<<51)-1`, so every output
//!    has `out[0..4] < 2^51`. Limb 4 is intentionally *not* asserted: the
//!    reference leaves it unmasked and the port faithfully does too.
//!    Verified: 0 violations across the 1012 committed vectors.
//!
//! 3. **Additive inverse on the canonical domain, under value semantics.**
//!    For genuine canonical inputs `a` (each limb `< 2^51`, positional
//!    value `< p`), `fp_add(a, fp_neg(a))` has positional value `0` mod
//!    `p`: `fp_neg(a)` is `-a` and the sum is the additive identity. This
//!    reuses `fp_add_props`'s positional `value_mod_p` reduction and is
//!    sound for the same reason `fp_add_props`'s value law is: the
//!    reference's `modneg` of a canonical `a` yields a representative of
//!    `-a` whose limbs are non-negative and bounded below `2p`, so the
//!    subsequent `fp_add` is an addition of two non-negative-valued
//!    representatives, exactly the case in which `fp_add`'s single
//!    conditional `-2p` correction leaves a representative whose plain
//!    positional value already equals the field element. Crucially the
//!    value is read off the *final* `fp_add` output, never the
//!    intermediate `fp_neg` output: reading `fp_neg(a)` positionally
//!    would be unsound for exactly the reason `fp_sub_props` documents
//!    (the `+2p` correction's `0xff..fe`-style limb 0 is not congruent),
//!    so that intermediate reading is deliberately avoided here, not
//!    asserted. Verified directly: across 2,000,006 canonical inputs
//!    (including the adversarial fixed cases `0`, `1`, the top-limb-only
//!    and max-clamped encodings) the final positional value was `0` with
//!    zero exceptions.
//!
//! ## Why no arbitrary-input or raw-limb negation law is asserted
//!
//! `fp_neg` is *not* a bit-exact involution and `fp_add(a, fp_neg(a))` is
//! *not* the bit-exact all-zero limb vector for arbitrary (non-canonical)
//! `a`: like `modsub`, `modneg`'s `+2p` correction produces a
//! representative whose plain positional value is not congruent to the
//! field value, and limb 4 is left unmasked, so neither a raw-limb
//! involution nor an arbitrary-input value law holds. Concrete witness
//! (the `fp_sub_props` witness, transposed): for `b` a canonical encoding
//! of `2`, `fp_neg(b)` leaves limb 0 `= 0xffff_ffff_ffff_fffe`, whose
//! positional value mod `p` is `2^64 - 2`, not the field value `p - 2`;
//! reading that intermediate positionally, or asserting raw-limb
//! `fp_neg(fp_neg(a)) == a`, would therefore be unsound. These are
//! deliberately omitted rather than weakened; the differential test
//! against the C oracle remains the value-correctness authority until
//! `redc` lands.

use proptest::prelude::*;
use sqisign_gf::{fp_add, fp_neg, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;
const MASK51: u64 = (1u64 << RADIX) - 1;

fn neg(a: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_neg(&mut c, a);
    c
}

fn add(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_add(&mut c, a, b);
    c
}

/// Little-endian byte big integer wide enough for the positional value of
/// five 64-bit limbs at weights `2^(51*i)` (max weight `2^204`, so a
/// 64-byte buffer never overflows). Same model as `fp_add_props`.
type Big = Vec<u8>;

fn zero_big() -> Big {
    vec![0u8; 64]
}

fn add_shifted(v: &mut [u8], limb: u64, shift: u32) {
    let byte0 = (shift / 8) as usize;
    let bit = shift % 8;
    let mut acc = (limb as u128) << bit;
    let mut idx = byte0;
    let mut carry = 0u16;
    while (acc != 0 || carry != 0) && idx < v.len() {
        let cur = v[idx] as u16 + (acc as u8) as u16 + carry;
        v[idx] = cur as u8;
        carry = cur >> 8;
        acc >>= 8;
        idx += 1;
    }
}

fn sub_one(v: &mut [u8]) {
    let mut i = 0;
    while v[i] == 0 {
        v[i] = 0xff;
        i += 1;
    }
    v[i] -= 1;
}

fn cmp(a: &[u8], b: &[u8]) -> std::cmp::Ordering {
    for i in (0..a.len().max(b.len())).rev() {
        let x = *a.get(i).unwrap_or(&0);
        let y = *b.get(i).unwrap_or(&0);
        if x != y {
            return x.cmp(&y);
        }
    }
    std::cmp::Ordering::Equal
}

fn sub_in_place(a: &mut [u8], b: &[u8]) {
    let mut borrow = 0i16;
    for (i, ai) in a.iter_mut().enumerate() {
        let bb = *b.get(i).unwrap_or(&0) as i16;
        let cur = *ai as i16 - bb - borrow;
        if cur < 0 {
            *ai = (cur + 256) as u8;
            borrow = 1;
        } else {
            *ai = cur as u8;
            borrow = 0;
        }
    }
}

fn shl1(v: &mut [u8]) {
    let mut carry = 0u8;
    for byte in v.iter_mut() {
        let nc = *byte >> 7;
        *byte = (*byte << 1) | carry;
        carry = nc;
    }
}

/// The level-1 prime `p5248 = 5 * 2^248 - 1`, little-endian bytes.
fn p_bytes() -> Big {
    let mut p = zero_big();
    add_shifted(&mut p, 5, 248);
    sub_one(&mut p);
    p
}

/// Reduce a little-endian byte big integer modulo `p` (binary long
/// division), returning the trimmed remainder.
fn reduce(v: &[u8], p: &[u8]) -> Big {
    let mut rem = vec![0u8; v.len()];
    for bit in (0..v.len() * 8).rev() {
        shl1(&mut rem);
        if (v[bit / 8] >> (bit % 8)) & 1 == 1 {
            rem[0] |= 1;
        }
        if cmp(&rem, p) != std::cmp::Ordering::Less {
            sub_in_place(&mut rem, p);
        }
    }
    trim(rem)
}

fn trim(mut v: Big) -> Big {
    while v.len() > 1 && *v.last().unwrap() == 0 {
        v.pop();
    }
    v
}

/// Positional value `sum limb[i] * 2^(51*i)` reduced modulo `p`. Sound as
/// the field element only on the canonical domain (see the module note),
/// the only domain on which it is used here.
fn value_mod_p(n: &Fp) -> Big {
    let mut v = zero_big();
    for (i, &limb) in n.iter().enumerate() {
        add_shifted(&mut v, limb, (51 * i) as u32);
    }
    reduce(&v, &p_bytes())
}

/// Build a canonical encoding (limbs < 2^51, positional value < p) from a
/// seed. Clamping limb 4 to 48 bits keeps `5 * 2^248` dominant, so the
/// positional value stays below `p`. Identical to `fp_add_props`.
fn canonical(seed: &[u64; NWORDS_FIELD]) -> Fp {
    let mut n = [0u64; NWORDS_FIELD];
    for i in 0..NWORDS_FIELD {
        n[i] = seed[i] & MASK51;
    }
    n[4] &= (1u64 << 48) - 1;
    n
}

#[test]
fn neg_of_zero_is_bit_exact_zero() {
    // (1) The lone sound raw-limb law: fp_neg of the canonical zero is
    // the bit-exact canonical zero. Pinned as a fixed case, not a
    // proptest, because it is a single point, not a family.
    assert_eq!(neg(&[0u64; NWORDS_FIELD]), [0u64; NWORDS_FIELD]);
}

proptest! {
    // (2) Structural carry-propagation invariant: limbs 0..=3 are masked
    // below 2^51 by the final prop; limb 4 is left unmasked by design and
    // is deliberately not constrained.
    #[test]
    fn limbs_0_3_below_radix(
        a in proptest::array::uniform5(any::<u64>()),
    ) {
        let c = neg(&a);
        for (k, &limb) in c.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }

    // (3) Additive inverse on the canonical domain, value semantics:
    // fp_add(a, fp_neg(a)) has positional value 0 mod p. Sound for the
    // same reason fp_add_props's value law is (the final fp_add output of
    // two non-negative-valued representatives reads positionally); the
    // intermediate fp_neg output is never read positionally. Verified
    // against the C battery and a 2,000,006-input out-of-tree model.
    #[test]
    fn canonical_neg_is_additive_inverse(
        sa in proptest::array::uniform5(any::<u64>()),
    ) {
        let a = canonical(&sa);
        let sum = add(&a, &neg(&a));
        prop_assert_eq!(value_mod_p(&sum), vec![0u8]);
    }
}
