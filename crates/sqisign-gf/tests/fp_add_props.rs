//! Property tests for `fp_add`.
//!
//! `fp_add` operates on a **redundant, non-canonical** radix-2^51 form: a
//! residue class has many limb encodings and `modadd` reduces only to
//! "less than 2p", leaving limb 4 unmasked. Raw-limb equality is therefore
//! *not* a sound notion of field equality, and the reference's own
//! equality (`modcmp`) is not ported yet. So the properties asserted here
//! are exactly the ones that are sound without it, established by study of
//! `modadd`/`prop` and cross-checked against the C-derived vectors (the
//! canonical correctness check remains the differential test):
//!
//! 1. **Commutativity, bit-exact, for arbitrary limb inputs.** `modadd`
//!    forms `n[i] = a[i] + b[i]` (symmetric in `a`,`b`) then runs a
//!    deterministic `prop`/correction sequence independent of operand
//!    order. So `fp_add(a,b)` and `fp_add(b,a)` are identical *at the
//!    raw-limb level*, not merely congruent. The strongest sound raw-limb
//!    law.
//!
//! 2. **Structural carry-propagation invariant, for arbitrary inputs.**
//!    The final `prop` masks limbs 0..=3 with `(1<<51)-1`, so every output
//!    has `out[0..4] < 2^51`. Limb 4 is intentionally *not* asserted: the
//!    reference leaves it unmasked and the port faithfully does too.
//!
//! 3. **Field laws on the canonical domain.** When inputs are genuine
//!    canonical encodings (each limb `< 2^51`, positional value `< p`),
//!    the positional value `V = sum limb[i] * 2^(51*i) mod p` is the field
//!    element and `modadd` realizes `(a + b) mod p`; on that domain the
//!    value semantics are commutative and associative. These do *not*
//!    hold for arbitrary non-canonical limb garbage (which the
//!    differential vectors also exercise), which is exactly why field
//!    equality needs `modcmp`, not this positional reading.

use proptest::prelude::*;
use sqisign_gf::{fp_add, Fp, NWORDS_FIELD};

const RADIX: u32 = 51;
const MASK51: u64 = (1u64 << RADIX) - 1;

fn add(a: &Fp, b: &Fp) -> Fp {
    let mut c = [0u64; NWORDS_FIELD];
    fp_add(&mut c, a, b);
    c
}

/// Little-endian byte big integer wide enough for the positional value of
/// five 64-bit limbs at weights `2^(51*i)` (max weight `2^204`, so a
/// 64-byte buffer never overflows).
type Big = Vec<u8>;

fn zero_big() -> Big {
    vec![0u8; 64]
}

/// Add `limb << shift` into a little-endian byte big integer in place.
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

/// Positional value `sum limb[i] * 2^(51*i)` reduced modulo `p`. Defined
/// for any limbs; equals the field element only on the canonical domain
/// (see the module note), which is the only domain on which it is used to
/// assert field laws.
fn value_mod_p(n: &Fp) -> Big {
    let mut v = zero_big();
    for (i, &limb) in n.iter().enumerate() {
        add_shifted(&mut v, limb, (51 * i) as u32);
    }
    reduce(&v, &p_bytes())
}

/// Sum of two reduced residues, reduced again modulo `p`.
fn add_mod_p(a: &[u8], b: &[u8]) -> Big {
    let mut s = zero_big();
    for (i, &x) in a.iter().enumerate() {
        add_shifted(&mut s, x as u64, (8 * i) as u32);
    }
    for (i, &x) in b.iter().enumerate() {
        add_shifted(&mut s, x as u64, (8 * i) as u32);
    }
    reduce(&s, &p_bytes())
}

/// Build a canonical encoding (limbs < 2^51, positional value < p) from a
/// seed. Clamping limb 4 to 48 bits keeps `5 * 2^248` dominant, so the
/// positional value stays below `p`.
fn canonical(seed: &[u64; NWORDS_FIELD]) -> Fp {
    let mut n = [0u64; NWORDS_FIELD];
    for i in 0..NWORDS_FIELD {
        n[i] = seed[i] & MASK51;
    }
    n[4] &= (1u64 << 48) - 1;
    n
}

proptest! {
    // (1) Commutativity, bit-exact, arbitrary (possibly non-canonical)
    // limb inputs. Sound: modadd's limbwise add is symmetric and the
    // prop/correction tail is operand-order-independent.
    #[test]
    fn commutative_bit_exact(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        prop_assert_eq!(add(&a, &b), add(&b, &a));
    }

    // (2) Structural carry-propagation invariant: limbs 0..=3 are masked
    // below 2^51 by the final prop; limb 4 is left unmasked by design and
    // is deliberately not constrained.
    #[test]
    fn limbs_0_3_below_radix(
        a in proptest::array::uniform5(any::<u64>()),
        b in proptest::array::uniform5(any::<u64>()),
    ) {
        let c = add(&a, &b);
        for (k, &limb) in c.iter().take(4).enumerate() {
            prop_assert!(limb < (1u64 << RADIX), "limb {k} = {limb:#x} >= 2^51");
        }
    }

    // (3a) Field value identity on the canonical domain: modadd realizes
    // (a + b) mod p when inputs are genuine canonical encodings.
    #[test]
    fn canonical_value_is_sum_mod_p(
        sa in proptest::array::uniform5(any::<u64>()),
        sb in proptest::array::uniform5(any::<u64>()),
    ) {
        let a = canonical(&sa);
        let b = canonical(&sb);
        let c = add(&a, &b);
        prop_assert_eq!(value_mod_p(&c), add_mod_p(&value_mod_p(&a), &value_mod_p(&b)));
    }

    // (3b) Associativity on the canonical domain, under value semantics.
    #[test]
    fn canonical_value_associative(
        sa in proptest::array::uniform5(any::<u64>()),
        sb in proptest::array::uniform5(any::<u64>()),
        sc in proptest::array::uniform5(any::<u64>()),
    ) {
        let a = canonical(&sa);
        let b = canonical(&sb);
        let c = canonical(&sc);
        prop_assert_eq!(
            value_mod_p(&add(&add(&a, &b), &c)),
            value_mod_p(&add(&a, &add(&b, &c)))
        );
    }
}
