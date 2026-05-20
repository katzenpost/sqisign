//! SQIsign `ec`: elliptic curve point operations on the Kummer line and
//! Jacobian coordinates, plus degree-2 and degree-4 isogeny construction
//! and evaluation.
//!
//! Mirrors `the-sqisign/src/ec`. Phase 1 unit 4 ports the core
//! curve and point arithmetic from `lvlx/ec.c`, the Jacobian arithmetic
//! from `lvlx/ec_jac.c`, the short isogeny primitives from `lvlx/xeval.c`
//! and `lvlx/xisog.c`, the long 2-isogeny chains and curve isomorphisms
//! from `lvlx/isog_chains.c`, the torsion basis generation and y-recovery
//! from `lvlx/basis.c`, and the Weil and reduced Tate pairings (and the
//! E[2^e] discrete logs derived from them) from `lvlx/biextension.c`.
//! The precomputed E0 torsion basis constants `BASIS_E0_PX` and
//! `BASIS_E0_QX` from `precomp/ref/lvl1/e0_basis.c` are also included.
//!
//! The differential boundary is the raw fp2 limb representation: every
//! ported function replays the recorded C output bit-for-bit on the
//! committed vector battery, identical in spirit to the gf port.
//!
//! ## Defects observed
//!
//! - `xeval_4_singular` and `xisog_4_singular` are declared in
//!   `the-sqisign/src/ec/ref/include/isog.h` but never defined
//!   anywhere in the upstream tree (no callers either): they are dead
//!   declarations. The port skips them; if upstream ever adds bodies,
//!   we add the boundaries then.
//! - `ec_is_equal` (in `lvlx/ec.c`) returns the value of
//!   `(l_zero & r_zero) | (~l_zero & ~r_zero * lr_equal)`, where C
//!   operator precedence binds `*` more tightly than `&`. The intended
//!   expression is almost certainly `~l_zero & ~r_zero & lr_equal`; the
//!   `*` is a typo. For all-ones masks the typo collapses
//!   `0xFFFFFFFF * 0xFFFFFFFF` to `1`, so a positive equality result
//!   comes back as `1` rather than the all-ones mask. The reference's
//!   `assert(ec_is_equal(...))` callers in `biextension.c` are
//!   unaffected because C `assert` only checks for non-zero; our Rust
//!   `debug_assert!` for the same calls likewise checks `!= 0` rather
//!   than `== 0xFFFFFFFF` so we faithfully mirror the C semantics. A
//!   call site that relies on the all-ones contract would be broken by
//!   the same typo on the C side. Not patched here; the literal port is
//!   the contract.
//! - `ec_dlog_2_tate_to_full` is declared in
//!   `the-sqisign/src/ec/ref/include/biextension.h` but never
//!   defined anywhere in the upstream tree. Skipped on the same
//!   grounds as `xeval_4_singular` above.

#![forbid(unsafe_code)]
#![allow(non_snake_case)]

use sqisign_gf::{
    fp2_add, fp2_add_one, fp2_batched_inv, fp2_copy, fp2_cswap, fp2_inv, fp2_is_equal, fp2_is_one,
    fp2_is_square, fp2_is_zero, fp2_mul, fp2_mul_small, fp2_neg, fp2_select, fp2_set_one,
    fp2_set_small, fp2_set_zero, fp2_sqr, fp2_sqrt, fp2_sqrt_verify, fp2_sub, fp_add, fp_copy,
    fp_div3, fp_is_square, fp_neg, fp_set_one, fp_set_small, fp_sub, Fp, Fp2, NWORDS_FIELD,
};
use sqisign_mp::{mp_add, mp_shiftr, mp_sub, multiple_mp_shiftl, select_ct, swap_ct};

// ---------------------------------------------------------------------------
// ec_params.c (lvl1)
// ---------------------------------------------------------------------------

/// The power of two dividing `p + 1`. Mirrors the C macro
/// `TORSION_EVEN_POWER` defined in
/// `the-sqisign/src/precomp/ref/lvl1/include/ec_params.h`.
pub const TORSION_EVEN_POWER: usize = 248;

/// The odd cofactor `(p + 1) / 2^TORSION_EVEN_POWER`, kept as a
/// single-element array to mirror the C declaration
/// `const digit_t p_cofactor_for_2f[1] = {5};` in
/// `the-sqisign/src/precomp/ref/lvl1/ec_params.c`.
pub const P_COFACTOR_FOR_2F: [u64; 1] = [5];

/// Bitlength of `P_COFACTOR_FOR_2F`. Mirrors
/// `#define P_COFACTOR_FOR_2F_BITLENGTH 3`.
pub const P_COFACTOR_FOR_2F_BITLENGTH: u32 = 3;

/// Number of limbs in an order-sized scalar at level 1
/// (`fp_constants.h`).
pub const NWORDS_ORDER: usize = 4;

/// Bit width of the field at level 1 (`fp_constants.h`).
pub const BITS: usize = 256;

/// Radix of the multiprecision representation (matches `tutil.h`
/// `RADIX = 64`).
pub const RADIX: usize = 64;

/// Log2 of `RADIX`, mirroring `LOG2RADIX = 6`.
pub const LOG2RADIX: usize = 6;

// ---------------------------------------------------------------------------
// Data structures (ec.h, isog.h)
// ---------------------------------------------------------------------------

/// Projective point on the Kummer line `E / pm 1` in Montgomery
/// coordinates `(X : Z)`. Mirrors
/// `struct ec_point_t { fp2_t x; fp2_t z; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcPoint {
    pub x: Fp2,
    pub z: Fp2,
}

impl EcPoint {
    /// All-zero placeholder useful for stack allocation; not a valid
    /// projective point (use [`ec_point_init`] for the identity).
    pub const fn zero() -> Self {
        EcPoint {
            x: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            z: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
        }
    }
}

/// Projective point in Jacobian Montgomery coordinates `(X : Y : Z)`.
/// Mirrors `struct jac_point_t { fp2_t x; fp2_t y; fp2_t z; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JacPoint {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
}

impl JacPoint {
    pub const fn zero() -> Self {
        JacPoint {
            x: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            y: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            z: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
        }
    }
}

/// Addition components: three values `(u, v, w)` such that
/// `x(P+Q) = (u-v : w)` and `x(P-Q) = (u+v : w)`. Mirrors
/// `struct add_components_t { fp2_t u; fp2_t v; fp2_t w; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AddComponents {
    pub u: Fp2,
    pub v: Fp2,
    pub w: Fp2,
}

impl AddComponents {
    pub const fn zero() -> Self {
        AddComponents {
            u: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            v: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            w: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
        }
    }
}

/// A torsion-subgroup basis: a pair of points and their difference.
/// Mirrors
/// `struct ec_basis_t { ec_point_t P; ec_point_t Q; ec_point_t PmQ; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcBasis {
    pub P: EcPoint,
    pub Q: EcPoint,
    pub PmQ: EcPoint,
}

impl EcBasis {
    pub const fn zero() -> Self {
        EcBasis {
            P: EcPoint::zero(),
            Q: EcPoint::zero(),
            PmQ: EcPoint::zero(),
        }
    }
}

/// An elliptic curve in projective Montgomery form. Mirrors
/// `struct ec_curve_t { fp2_t A; fp2_t C; ec_point_t A24;
///                     bool is_A24_computed_and_normalized; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcCurve {
    pub A: Fp2,
    pub C: Fp2,
    pub A24: EcPoint,
    pub is_A24_computed_and_normalized: bool,
}

impl EcCurve {
    pub const fn zero() -> Self {
        EcCurve {
            A: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            C: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            A24: EcPoint::zero(),
            is_A24_computed_and_normalized: false,
        }
    }
}

/// KPS structure for a degree-2 isogeny. Mirrors
/// `typedef struct { ec_point_t K; } ec_kps2_t;`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcKps2 {
    pub K: EcPoint,
}

impl EcKps2 {
    pub const fn zero() -> Self {
        EcKps2 { K: EcPoint::zero() }
    }
}

/// KPS structure for a degree-4 isogeny. Mirrors
/// `typedef struct { ec_point_t K[3]; } ec_kps4_t;`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcKps4 {
    pub K: [EcPoint; 3],
}

impl EcKps4 {
    pub const fn zero() -> Self {
        EcKps4 {
            K: [EcPoint::zero(); 3],
        }
    }
}

// ---------------------------------------------------------------------------
// Small helpers (ec.h inline definitions)
// ---------------------------------------------------------------------------

/// Mirrors the inline `copy_point(P, Q)` in `ec.h`.
pub fn copy_point(dst: &mut EcPoint, src: &EcPoint) {
    fp2_copy(&mut dst.x, &src.x);
    fp2_copy(&mut dst.z, &src.z);
}

/// Mirrors the inline `copy_basis(B1, B0)` in `ec.h`.
pub fn copy_basis(dst: &mut EcBasis, src: &EcBasis) {
    copy_point(&mut dst.P, &src.P);
    copy_point(&mut dst.Q, &src.Q);
    copy_point(&mut dst.PmQ, &src.PmQ);
}

/// Mirrors the inline `copy_curve(E1, E2)` in `ec.h`.
pub fn copy_curve(dst: &mut EcCurve, src: &EcCurve) {
    fp2_copy(&mut dst.A, &src.A);
    fp2_copy(&mut dst.C, &src.C);
    dst.is_A24_computed_and_normalized = src.is_A24_computed_and_normalized;
    copy_point(&mut dst.A24, &src.A24);
}

/// Mirrors the inline `AC_to_A24(A24, E)` in `ec.h`.
pub fn ac_to_a24(a24: &mut EcPoint, e: &EcCurve) {
    if e.is_A24_computed_and_normalized {
        copy_point(a24, &e.A24);
        return;
    }
    // A24 = (A+2C : 4C)
    fp2_add(&mut a24.z, &e.C, &e.C);
    let az = a24.z;
    fp2_add(&mut a24.x, &e.A, &az);
    let az = a24.z;
    fp2_add(&mut a24.z, &az, &az);
}

/// Mirrors the inline `A24_to_AC(E, A24)` in `ec.h`.
pub fn a24_to_ac(e: &mut EcCurve, a24: &EcPoint) {
    // (A : C) = ((A+2C)*2 - 4C : 4C)
    fp2_add(&mut e.A, &a24.x, &a24.x);
    let a = e.A;
    fp2_sub(&mut e.A, &a, &a24.z);
    let a = e.A;
    fp2_add(&mut e.A, &a, &a);
    fp2_copy(&mut e.C, &a24.z);
}

// ---------------------------------------------------------------------------
// ec.c: init, normalize, predicates, j-invariant
// ---------------------------------------------------------------------------

/// Initialize point as the identity element `(1 : 0)`.
pub fn ec_point_init(p: &mut EcPoint) {
    fp2_set_one(&mut p.x);
    fp2_set_zero(&mut p.z);
}

/// Initialize the curve struct: `A = 0`, `C = 1`, `A24 = (1 : 0)`,
/// `is_A24_computed_and_normalized = false`.
pub fn ec_curve_init(e: &mut EcCurve) {
    fp2_set_zero(&mut e.A);
    fp2_set_one(&mut e.C);
    ec_point_init(&mut e.A24);
    e.is_A24_computed_and_normalized = false;
}

/// Branchless constant-time select between two points. Mirrors
/// `select_point`. `option = 0` selects `P1`; `option = 0xFF..FF`
/// selects `P2`.
pub fn select_point(q: &mut EcPoint, p1: &EcPoint, p2: &EcPoint, option: u64) {
    fp2_select(&mut q.x, &p1.x, &p2.x, option as u32);
    fp2_select(&mut q.z, &p1.z, &p2.z, option as u32);
}

/// Branchless constant-time conditional swap. Mirrors `cswap_points`.
pub fn cswap_points(p: &mut EcPoint, q: &mut EcPoint, option: u64) {
    fp2_cswap(&mut p.x, &mut q.x, option as u32);
    fp2_cswap(&mut p.z, &mut q.z, option as u32);
}

/// Mirrors `ec_normalize_point`: in-place `(X : Z) -> (X/Z : 1)`.
pub fn ec_normalize_point(p: &mut EcPoint) {
    fp2_inv(&mut p.z);
    let px = p.x;
    let pz = p.z;
    fp2_mul(&mut p.x, &px, &pz);
    fp2_set_one(&mut p.z);
}

/// Mirrors `ec_normalize_curve`: in-place `(A : C) -> (A/C : 1)`.
pub fn ec_normalize_curve(e: &mut EcCurve) {
    fp2_inv(&mut e.C);
    let ea = e.A;
    let ec = e.C;
    fp2_mul(&mut e.A, &ea, &ec);
    fp2_set_one(&mut e.C);
}

/// Mirrors `ec_curve_normalize_A24`: ensures `(A+2)/4 : 1` is cached.
pub fn ec_curve_normalize_a24(e: &mut EcCurve) {
    if !e.is_A24_computed_and_normalized {
        let mut tmp = EcPoint::zero();
        ac_to_a24(&mut tmp, e);
        e.A24 = tmp;
        ec_normalize_point(&mut e.A24);
        e.is_A24_computed_and_normalized = true;
    }
    debug_assert!(fp2_is_one(&e.A24.z) == 0xFFFFFFFF);
}

/// Mirrors `ec_normalize_curve_and_A24`: brings both `(A : C)` and
/// `(A+2 : 4C)` to normalized form.
pub fn ec_normalize_curve_and_a24(e: &mut EcCurve) {
    if fp2_is_one(&e.C) != 0xFFFFFFFF {
        ec_normalize_curve(e);
    }
    if !e.is_A24_computed_and_normalized {
        // re(A24.x) = re(A) + 1
        let ea = e.A;
        fp2_add_one(&mut e.A24.x, &ea);
        // re(A24.x) = re(A) + 2
        let xx = e.A24.x;
        fp2_add_one(&mut e.A24.x, &xx);
        // im(A24.x) = im(A)
        fp_copy(&mut e.A24.x.im, &e.A.im);
        // (A + 2) / 2
        let xx = e.A24.x;
        sqisign_gf::fp2_half(&mut e.A24.x, &xx);
        // (A + 2) / 4
        let xx = e.A24.x;
        sqisign_gf::fp2_half(&mut e.A24.x, &xx);
        fp2_set_one(&mut e.A24.z);
        e.is_A24_computed_and_normalized = true;
    }
}

/// Mirrors `ec_is_zero`.
pub fn ec_is_zero(p: &EcPoint) -> u32 {
    fp2_is_zero(&p.z)
}

/// Mirrors `ec_has_zero_coordinate`.
pub fn ec_has_zero_coordinate(p: &EcPoint) -> u32 {
    fp2_is_zero(&p.x) | fp2_is_zero(&p.z)
}

/// Mirrors `ec_is_equal`. The upstream expression
/// `(l_zero & r_zero) | (~l_zero & ~r_zero * lr_equal)` is reproduced
/// literally (`*` binds tighter than `&`).
pub fn ec_is_equal(p: &EcPoint, q: &EcPoint) -> u32 {
    let l_zero = ec_is_zero(p);
    let r_zero = ec_is_zero(q);

    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    fp2_mul(&mut t0, &p.x, &q.z);
    fp2_mul(&mut t1, &p.z, &q.x);
    let lr_equal = fp2_is_equal(&t0, &t1);

    // Wrapping mul reproduces unsigned C multiplication.
    (l_zero & r_zero) | ((!l_zero) & (!r_zero).wrapping_mul(lr_equal))
}

/// Mirrors `ec_is_two_torsion`.
pub fn ec_is_two_torsion(p: &EcPoint, e: &EcCurve) -> u32 {
    if ec_is_zero(p) == 0xFFFFFFFF {
        return 0;
    }

    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    fp2_add(&mut t0, &p.x, &p.z);
    let tt = t0;
    fp2_sqr(&mut t0, &tt);
    fp2_sub(&mut t1, &p.x, &p.z);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    fp2_sub(&mut t2, &t0, &t1);
    let ta = t0;
    let tb = t1;
    fp2_add(&mut t1, &ta, &tb);
    let tt = t2;
    fp2_mul(&mut t2, &tt, &e.A);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &e.C);
    let tt = t1;
    fp2_add(&mut t1, &tt, &tt);
    fp2_add(&mut t0, &t1, &t2);

    let x_is_zero = fp2_is_zero(&p.x);
    let tmp_is_zero = fp2_is_zero(&t0);
    x_is_zero | tmp_is_zero
}

/// Mirrors `ec_is_four_torsion`.
pub fn ec_is_four_torsion(p: &EcPoint, e: &EcCurve) -> u32 {
    let mut test = EcPoint::zero();
    x_dbl_a24(&mut test, p, &e.A24, e.is_A24_computed_and_normalized);
    ec_is_two_torsion(&test, e)
}

/// Mirrors `ec_is_basis_four_torsion`.
pub fn ec_is_basis_four_torsion(b: &EcBasis, e: &EcCurve) -> u32 {
    let mut p2 = EcPoint::zero();
    let mut q2 = EcPoint::zero();
    x_dbl_a24(&mut p2, &b.P, &e.A24, e.is_A24_computed_and_normalized);
    x_dbl_a24(&mut q2, &b.Q, &e.A24, e.is_A24_computed_and_normalized);
    ec_is_two_torsion(&p2, e) & ec_is_two_torsion(&q2, e) & !ec_is_equal(&p2, &q2)
}

/// Mirrors `ec_curve_verify_A`: `1` if the Montgomery coefficient `A`
/// is valid (`A^2 - 4 != 0`), `0` otherwise.
pub fn ec_curve_verify_a(a: &Fp2) -> i32 {
    let mut t = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_set_one(&mut t);
    let tr = t.re;
    fp_add(&mut t.re, &tr, &tr); // t = 2
    if fp2_is_equal(a, &t) == 0xFFFFFFFF {
        return 0;
    }
    let tr = t.re;
    fp_neg(&mut t.re, &tr); // t = -2
    if fp2_is_equal(a, &t) == 0xFFFFFFFF {
        return 0;
    }
    1
}

/// Mirrors `ec_curve_init_from_A`.
pub fn ec_curve_init_from_a(e: &mut EcCurve, a: &Fp2) -> i32 {
    ec_curve_init(e);
    fp2_copy(&mut e.A, a);
    ec_curve_verify_a(a)
}

/// Mirrors `ec_j_inv`: writes the j-invariant of the curve into
/// `j_inv`.
pub fn ec_j_inv(j_inv: &mut Fp2, curve: &EcCurve) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;

    fp2_sqr(&mut t1, &curve.C);
    fp2_sqr(j_inv, &curve.A);
    fp2_add(&mut t0, &t1, &t1);
    let jj = *j_inv;
    let tt = t0;
    fp2_sub(&mut t0, &jj, &tt);
    let tt = t0;
    fp2_sub(&mut t0, &tt, &t1);
    fp2_sub(j_inv, &t0, &t1);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    let jj = *j_inv;
    fp2_mul(j_inv, &jj, &t1);
    let tt = t0;
    fp2_add(&mut t0, &tt, &tt);
    let tt = t0;
    fp2_add(&mut t0, &tt, &tt);
    fp2_sqr(&mut t1, &t0);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t1);
    let tt = t0;
    fp2_add(&mut t0, &tt, &tt);
    let tt = t0;
    fp2_add(&mut t0, &tt, &tt);
    fp2_inv(j_inv);
    let jj = *j_inv;
    fp2_mul(j_inv, &t0, &jj);
}

// ---------------------------------------------------------------------------
// ec.c: doubling, addition, scalar ladder
// ---------------------------------------------------------------------------

/// Doubling on the curve `E0` with `(A : C) = (0 : 1)`. Mirrors
/// `xDBL_E0`.
pub fn x_dbl_e0(q: &mut EcPoint, p: &EcPoint) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    let tt = t0;
    fp2_sqr(&mut t0, &tt);
    fp2_sub(&mut t1, &p.x, &p.z);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    fp2_sub(&mut t2, &t0, &t1);
    let tt = t1;
    fp2_add(&mut t1, &tt, &tt);
    fp2_mul(&mut q.x, &t0, &t1);
    fp2_add(&mut q.z, &t1, &t2);
    let qz = q.z;
    fp2_mul(&mut q.z, &qz, &t2);
}

/// Doubling on a general Montgomery curve with on-the-fly recovery of
/// `(A+2C, 4C)`. Mirrors `xDBL`. The reference takes the curve as
/// `const ec_point_t *AC` aliased onto the `(A, C)` prefix of an
/// `ec_curve_t`; the port takes the full `EcCurve` and reads `A` and
/// `C` directly.
pub fn x_dbl(q: &mut EcPoint, p: &EcPoint, ac: &EcCurve) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    let tt = t0;
    fp2_sqr(&mut t0, &tt);
    fp2_sub(&mut t1, &p.x, &p.z);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    fp2_sub(&mut t2, &t0, &t1);
    fp2_add(&mut t3, &ac.C, &ac.C);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &t3);
    let tt = t1;
    fp2_add(&mut t1, &tt, &tt);
    fp2_mul(&mut q.x, &t0, &t1);
    fp2_add(&mut t0, &t3, &ac.A);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t2);
    let tt = t0;
    fp2_add(&mut t0, &tt, &t1);
    fp2_mul(&mut q.z, &t0, &t2);
}

/// Doubling taking `A24 = (A+2C : 4C)` (or `(A+2C/4C : 1)` when
/// normalized). Mirrors `xDBL_A24`.
pub fn x_dbl_a24(q: &mut EcPoint, p: &EcPoint, a24: &EcPoint, a24_normalized: bool) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    let tt = t0;
    fp2_sqr(&mut t0, &tt);
    fp2_sub(&mut t1, &p.x, &p.z);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    fp2_sub(&mut t2, &t0, &t1);
    if !a24_normalized {
        let tt = t1;
        fp2_mul(&mut t1, &tt, &a24.z);
    }
    fp2_mul(&mut q.x, &t0, &t1);
    fp2_mul(&mut t0, &t2, &a24.x);
    let tt = t0;
    fp2_add(&mut t0, &tt, &t1);
    fp2_mul(&mut q.z, &t0, &t2);
}

/// Differential addition on the Kummer line. Mirrors `xADD`.
pub fn x_add(r: &mut EcPoint, p: &EcPoint, q: &EcPoint, pq: &EcPoint) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_add(&mut t2, &q.x, &q.z);
    fp2_sub(&mut t3, &q.x, &q.z);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t3);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &t2);
    fp2_add(&mut t2, &t0, &t1);
    fp2_sub(&mut t3, &t0, &t1);
    let tt = t2;
    fp2_sqr(&mut t2, &tt);
    let tt = t3;
    fp2_sqr(&mut t3, &tt);
    let tt = t2;
    fp2_mul(&mut t2, &pq.z, &tt);
    fp2_mul(&mut r.z, &pq.x, &t3);
    fp2_copy(&mut r.x, &t2);
}

/// Simultaneous double-and-differential-add. Mirrors `xDBLADD`.
pub fn x_dbl_add(
    r: &mut EcPoint,
    s: &mut EcPoint,
    p: &EcPoint,
    q: &EcPoint,
    pq: &EcPoint,
    a24: &EcPoint,
    a24_normalized: bool,
) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_sqr(&mut r.x, &t0);
    fp2_sub(&mut t2, &q.x, &q.z);
    fp2_add(&mut s.x, &q.x, &q.z);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t2);
    fp2_sqr(&mut r.z, &t1);
    let tt = t1;
    let sx = s.x;
    fp2_mul(&mut t1, &tt, &sx);
    let rx = r.x;
    let rz = r.z;
    fp2_sub(&mut t2, &rx, &rz);
    if !a24_normalized {
        let rz = r.z;
        fp2_mul(&mut r.z, &rz, &a24.z);
    }
    let rz = r.z;
    let rx = r.x;
    fp2_mul(&mut r.x, &rx, &rz);
    fp2_mul(&mut s.x, &a24.x, &t2);
    fp2_sub(&mut s.z, &t0, &t1);
    let rz = r.z;
    let sx = s.x;
    fp2_add(&mut r.z, &rz, &sx);
    fp2_add(&mut s.x, &t0, &t1);
    let rz = r.z;
    fp2_mul(&mut r.z, &rz, &t2);
    let sz = s.z;
    fp2_sqr(&mut s.z, &sz);
    let sx = s.x;
    fp2_sqr(&mut s.x, &sx);
    let sz = s.z;
    fp2_mul(&mut s.z, &sz, &pq.x);
    let sx = s.x;
    fp2_mul(&mut s.x, &sx, &pq.z);
}

/// The Montgomery ladder `Q = k * P`. Mirrors `xMUL`. Takes a scalar
/// slice of length `>= NWORDS_ORDER`; only the first `NWORDS_ORDER`
/// limbs are read.
pub fn x_mul(q: &mut EcPoint, p: &EcPoint, k: &[u64], kbits: i32, curve: &EcCurve) {
    let mut r0 = EcPoint::zero();
    let mut r1 = EcPoint::zero();
    let mut a24 = EcPoint::zero();
    let mut prevbit: u32 = 0;

    if !curve.is_A24_computed_and_normalized {
        // A24 = (A+2C : 4C)
        fp2_add(&mut a24.x, &curve.C, &curve.C);
        fp2_add(&mut a24.z, &a24.x, &a24.x);
        let ax = a24.x;
        fp2_add(&mut a24.x, &ax, &curve.A);
    } else {
        fp2_copy(&mut a24.x, &curve.A24.x);
        fp2_copy(&mut a24.z, &curve.A24.z);
        debug_assert!(fp2_is_one(&a24.z) == 0xFFFFFFFF);
    }

    // R0 <- (1 : 0), R1 <- P
    ec_point_init(&mut r0);
    fp2_copy(&mut r1.x, &p.x);
    fp2_copy(&mut r1.z, &p.z);

    for i in (0..kbits).rev() {
        let idx = (i as usize) >> LOG2RADIX;
        let shift = (i as usize) & (RADIX - 1);
        let bit = ((k[idx] >> shift) & 1) as u32;
        let swap = bit ^ prevbit;
        prevbit = bit;
        let mask = 0u64.wrapping_sub(swap as u64);

        cswap_points(&mut r0, &mut r1, mask);
        let r0_in = r0;
        let r1_in = r1;
        x_dbl_add(&mut r0, &mut r1, &r0_in, &r1_in, p, &a24, true);
    }
    let swap = prevbit;
    let mask = 0u64.wrapping_sub(swap as u64);
    cswap_points(&mut r0, &mut r1, mask);

    fp2_copy(&mut q.x, &r0.x);
    fp2_copy(&mut q.z, &r0.z);
}

/// The Montgomery biladder `S = k*P + l*Q`. Mirrors `xDBLMUL`. Returns
/// `0` if the input is invalid (any of `P`, `Q`, `PQ` has a zero
/// coordinate, or `P+Q` does), else `1`.
#[allow(clippy::too_many_arguments)]
pub fn x_dbl_mul(
    s_out: &mut EcPoint,
    p: &EcPoint,
    k: &[u64],
    q: &EcPoint,
    l: &[u64],
    pq: &EcPoint,
    kbits: i32,
    curve: &EcCurve,
) -> i32 {
    if ec_has_zero_coordinate(p) != 0
        || ec_has_zero_coordinate(q) != 0
        || ec_has_zero_coordinate(pq) != 0
    {
        return 0;
    }

    let bitk0 = k[0] & 1;
    let bitl0 = l[0] & 1;
    let mut maskk = 0u64.wrapping_sub(bitk0);
    let maskl = 0u64.wrapping_sub(bitl0);
    let mut sigma = [bitk0 ^ 1, bitl0 ^ 1];
    let evens = sigma[0] + sigma[1];
    let mevens = 0u64.wrapping_sub(evens & 1);

    sigma[0] &= mevens;
    sigma[1] = (sigma[1] & mevens) | (1 & !mevens);

    // Convert even scalars to odd
    let mut k_t = [0u64; NWORDS_ORDER];
    let mut l_t = [0u64; NWORDS_ORDER];
    let mut one = [0u64; NWORDS_ORDER];
    one[0] = 1;
    {
        let mut k_minus_one = [0u64; NWORDS_ORDER];
        let mut l_minus_one = [0u64; NWORDS_ORDER];
        let k_slice = &k[..NWORDS_ORDER];
        let l_slice = &l[..NWORDS_ORDER];
        mp_sub(&mut k_minus_one, k_slice, &one);
        mp_sub(&mut l_minus_one, l_slice, &one);
        let k_min = k_minus_one;
        let l_min = l_minus_one;
        select_ct(&mut k_t, &k_min, k_slice, maskk);
        select_ct(&mut l_t, &l_min, l_slice, maskl);
    }

    let mut r: Vec<u64> = vec![0; 2 * BITS];
    let mut pre_sigma: u64 = 0;
    for i in 0..(kbits as usize) {
        maskk = 0u64.wrapping_sub(sigma[0] ^ pre_sigma);
        swap_ct(&mut k_t, &mut l_t, maskk);

        let (bs1_ip1, bs2_ip1) = if i as i32 == kbits - 1 {
            (0u64, 0u64)
        } else {
            let a = mp_shiftr(&mut k_t, 1);
            let b = mp_shiftr(&mut l_t, 1);
            (a, b)
        };
        let bs1_i = k_t[0] & 1;
        let bs2_i = l_t[0] & 1;

        r[2 * i] = bs1_i ^ bs1_ip1;
        r[2 * i + 1] = bs2_i ^ bs2_ip1;

        pre_sigma = sigma[0];
        maskk = 0u64.wrapping_sub(r[2 * i + 1]);
        let mut temp = [0u64; 1];
        let s1_arr = [sigma[1]];
        let s0_arr = [sigma[0]];
        select_ct(&mut temp, &s0_arr, &s1_arr, maskk);
        let mut new_s1 = [0u64; 1];
        select_ct(&mut new_s1, &s1_arr, &s0_arr, maskk);
        sigma[1] = new_s1[0];
        sigma[0] = temp[0];
    }

    // Point initialization
    let mut r_pts = [EcPoint::zero(); 3];
    ec_point_init(&mut r_pts[0]);
    maskk = 0u64.wrapping_sub(sigma[0]);
    let mut r1_tmp = EcPoint::zero();
    select_point(&mut r1_tmp, p, q, maskk);
    let mut r2_tmp = EcPoint::zero();
    select_point(&mut r2_tmp, q, p, maskk);
    r_pts[1] = r1_tmp;
    r_pts[2] = r2_tmp;

    let mut diff1a = EcPoint::zero();
    let mut diff1b = EcPoint::zero();
    let mut diff2a = EcPoint::zero();
    let mut diff2b = EcPoint::zero();
    fp2_copy(&mut diff1a.x, &r_pts[1].x);
    fp2_copy(&mut diff1a.z, &r_pts[1].z);
    fp2_copy(&mut diff1b.x, &r_pts[2].x);
    fp2_copy(&mut diff1b.z, &r_pts[2].z);

    // DIFF2a <- P+Q, DIFF2b <- P-Q
    let r1_in = r_pts[1];
    let r2_in = r_pts[2];
    x_add(&mut r_pts[2], &r1_in, &r2_in, pq);
    if ec_has_zero_coordinate(&r_pts[2]) != 0 {
        return 0;
    }
    fp2_copy(&mut diff2a.x, &r_pts[2].x);
    fp2_copy(&mut diff2a.z, &r_pts[2].z);
    fp2_copy(&mut diff2b.x, &pq.x);
    fp2_copy(&mut diff2b.z, &pq.z);

    let a_is_zero = fp2_is_zero(&curve.A);

    for i in (0..kbits as usize).rev() {
        let h = r[2 * i] + r[2 * i + 1];
        let mut t = [EcPoint::zero(); 3];
        maskk = 0u64.wrapping_sub(h & 1);
        select_point(&mut t[0], &r_pts[0], &r_pts[1], maskk);
        maskk = 0u64.wrapping_sub(h >> 1);
        let t0_in = t[0];
        select_point(&mut t[0], &t0_in, &r_pts[2], maskk);
        if a_is_zero == 0xFFFFFFFF {
            let t0_in = t[0];
            x_dbl_e0(&mut t[0], &t0_in);
        } else {
            debug_assert!(fp2_is_one(&curve.A24.z) == 0xFFFFFFFF);
            let t0_in = t[0];
            x_dbl_a24(&mut t[0], &t0_in, &curve.A24, true);
        }

        maskk = 0u64.wrapping_sub(r[2 * i + 1]);
        select_point(&mut t[1], &r_pts[0], &r_pts[1], maskk);
        select_point(&mut t[2], &r_pts[1], &r_pts[2], maskk);

        cswap_points(&mut diff1a, &mut diff1b, maskk);
        let t1_in = t[1];
        let t2_in = t[2];
        x_add(&mut t[1], &t1_in, &t2_in, &diff1a);
        x_add(&mut t[2], &r_pts[0], &r_pts[2], &diff2a);

        maskk = 0u64.wrapping_sub(h & 1);
        cswap_points(&mut diff2a, &mut diff2b, maskk);

        copy_point(&mut r_pts[0], &t[0]);
        copy_point(&mut r_pts[1], &t[1]);
        copy_point(&mut r_pts[2], &t[2]);
    }

    select_point(s_out, &r_pts[0], &r_pts[1], mevens);
    maskk = 0u64.wrapping_sub(bitk0 & bitl0);
    let s_in = *s_out;
    select_point(s_out, &s_in, &r_pts[2], maskk);
    1
}

/// The 3-point Montgomery ladder `R = P + m * Q`. Mirrors
/// `ec_ladder3pt`.
pub fn ec_ladder3pt(
    r: &mut EcPoint,
    m: &[u64],
    p: &EcPoint,
    q: &EcPoint,
    pq: &EcPoint,
    e: &EcCurve,
) -> i32 {
    debug_assert!(e.is_A24_computed_and_normalized);
    if fp2_is_one(&e.A24.z) != 0xFFFFFFFF {
        return 0;
    }
    if ec_has_zero_coordinate(pq) != 0 {
        return 0;
    }

    let mut x0 = EcPoint::zero();
    let mut x1 = EcPoint::zero();
    let mut x2 = EcPoint::zero();
    copy_point(&mut x0, q);
    copy_point(&mut x1, p);
    copy_point(&mut x2, pq);

    for mi in m.iter().take(NWORDS_ORDER) {
        let mut t: u64 = 1;
        for _ in 0..RADIX {
            let cond = ((t & *mi) == 0) as u64;
            let mask = 0u64.wrapping_sub(cond);
            cswap_points(&mut x1, &mut x2, mask);
            let x0_in = x0;
            let x1_in = x1;
            x_dbl_add(&mut x0, &mut x1, &x0_in, &x1_in, &x2, &e.A24, true);
            cswap_points(&mut x1, &mut x2, mask);
            t = t.wrapping_shl(1);
        }
    }
    copy_point(r, &x1);
    1
}

/// Mirrors `ec_dbl`.
pub fn ec_dbl(res: &mut EcPoint, p: &EcPoint, curve: &EcCurve) {
    if curve.is_A24_computed_and_normalized {
        debug_assert!(fp2_is_one(&curve.A24.z) == 0xFFFFFFFF);
        x_dbl_a24(res, p, &curve.A24, true);
    } else {
        x_dbl(res, p, curve);
    }
}

/// Mirrors `ec_dbl_iter`. May mutate `curve` (normalizes `A24` when
/// `n > 50`).
pub fn ec_dbl_iter(res: &mut EcPoint, n: i32, p: &EcPoint, curve: &mut EcCurve) {
    if n == 0 {
        copy_point(res, p);
        return;
    }
    if n > 50 {
        ec_curve_normalize_a24(curve);
    }
    if curve.is_A24_computed_and_normalized {
        debug_assert!(fp2_is_one(&curve.A24.z) == 0xFFFFFFFF);
        x_dbl_a24(res, p, &curve.A24, true);
        for _ in 0..(n - 1) {
            debug_assert!(fp2_is_one(&curve.A24.z) == 0xFFFFFFFF);
            let r_in = *res;
            x_dbl_a24(res, &r_in, &curve.A24, true);
        }
    } else {
        x_dbl(res, p, curve);
        for _ in 0..(n - 1) {
            let r_in = *res;
            x_dbl(res, &r_in, curve);
        }
    }
}

/// Mirrors `ec_dbl_iter_basis`.
pub fn ec_dbl_iter_basis(res: &mut EcBasis, n: i32, b: &EcBasis, curve: &mut EcCurve) {
    ec_dbl_iter(&mut res.P, n, &b.P, curve);
    ec_dbl_iter(&mut res.Q, n, &b.Q, curve);
    ec_dbl_iter(&mut res.PmQ, n, &b.PmQ, curve);
}

/// Mirrors `ec_mul`. May mutate `curve` (normalizes `A24` when
/// `kbits > 50`).
pub fn ec_mul(res: &mut EcPoint, scalar: &[u64], kbits: i32, p: &EcPoint, curve: &mut EcCurve) {
    if kbits > 50 {
        ec_curve_normalize_a24(curve);
    }
    x_mul(res, p, scalar, kbits, curve);
}

/// Mirrors `ec_biscalar_mul`.
pub fn ec_biscalar_mul(
    res: &mut EcPoint,
    scalar_p: &[u64],
    scalar_q: &[u64],
    kbits: i32,
    pq: &EcBasis,
    curve: &EcCurve,
) -> i32 {
    if fp2_is_zero(&pq.PmQ.z) == 0xFFFFFFFF {
        return 0;
    }

    if kbits == 1 {
        if ec_is_two_torsion(&pq.P, curve) != 0xFFFFFFFF
            || ec_is_two_torsion(&pq.Q, curve) != 0xFFFFFFFF
            || ec_is_two_torsion(&pq.PmQ, curve) != 0xFFFFFFFF
        {
            return 0;
        }
        let b_p = scalar_p[0] & 1;
        let b_q = scalar_q[0] & 1;
        if b_p == 0 && b_q == 0 {
            ec_point_init(res);
        } else if b_p == 1 && b_q == 0 {
            copy_point(res, &pq.P);
        } else if b_p == 0 && b_q == 1 {
            copy_point(res, &pq.Q);
        } else {
            copy_point(res, &pq.PmQ);
        }
        return 1;
    }

    let mut e = EcCurve::zero();
    copy_curve(&mut e, curve);
    if fp2_is_zero(&curve.A) != 0xFFFFFFFF {
        ec_curve_normalize_a24(&mut e);
    }
    x_dbl_mul(res, &pq.P, scalar_p, &pq.Q, scalar_q, &pq.PmQ, kbits, &e)
}

// ---------------------------------------------------------------------------
// ec_jac.c: Jacobian point arithmetic
// ---------------------------------------------------------------------------

/// Mirrors `jac_init`: identity element `(0 : 1 : 0)`.
pub fn jac_init(p: &mut JacPoint) {
    fp2_set_zero(&mut p.x);
    fp2_set_one(&mut p.y);
    fp2_set_zero(&mut p.z);
}

/// Mirrors `copy_jac_point`.
pub fn copy_jac_point(dst: &mut JacPoint, src: &JacPoint) {
    fp2_copy(&mut dst.x, &src.x);
    fp2_copy(&mut dst.y, &src.y);
    fp2_copy(&mut dst.z, &src.z);
}

/// Mirrors `jac_is_equal`.
pub fn jac_is_equal(p: &JacPoint, q: &JacPoint) -> u32 {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;

    fp2_sqr(&mut t0, &q.z);
    fp2_mul(&mut t2, &p.x, &t0);
    fp2_sqr(&mut t1, &p.z);
    fp2_mul(&mut t3, &q.x, &t1);
    let tt = t2;
    fp2_sub(&mut t2, &tt, &t3);

    let tt = t0;
    fp2_mul(&mut t0, &tt, &q.z);
    let tt = t0;
    fp2_mul(&mut t0, &p.y, &tt);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &p.z);
    let tt = t1;
    fp2_mul(&mut t1, &q.y, &tt);
    let tt = t0;
    fp2_sub(&mut t0, &tt, &t1);

    fp2_is_zero(&t0) & fp2_is_zero(&t2)
}

/// Mirrors `jac_to_xz`.
pub fn jac_to_xz(p: &mut EcPoint, xy_p: &JacPoint) {
    fp2_copy(&mut p.x, &xy_p.x);
    fp2_copy(&mut p.z, &xy_p.z);
    let pz = p.z;
    fp2_sqr(&mut p.z, &pz);

    let mut one = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_set_one(&mut one);

    let c1 = fp2_is_zero(&p.x);
    let c2 = fp2_is_zero(&p.z);
    let px = p.x;
    fp2_select(&mut p.x, &px, &one, c1 & c2);
}

/// Mirrors `jac_to_ws`.
pub fn jac_to_ws(q: &mut JacPoint, t: &mut Fp2, ao3: &mut Fp2, p: &JacPoint, curve: &EcCurve) {
    let mut one_fp: Fp = [0u64; NWORDS_FIELD];
    fp_set_one(&mut one_fp);
    let mut a = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    if fp2_is_zero(&curve.A) != 0xFFFFFFFF {
        fp_div3(&mut ao3.re, &curve.A.re);
        fp_div3(&mut ao3.im, &curve.A.im);
        fp2_sqr(t, &p.z);
        fp2_mul(&mut q.x, ao3, t);
        let qx = q.x;
        fp2_add(&mut q.x, &qx, &p.x);
        let tt = *t;
        fp2_sqr(t, &tt);
        fp2_mul(&mut a, ao3, &curve.A);
        let ar = a.re;
        fp_sub(&mut a.re, &one_fp, &ar);
        let ai = a.im;
        fp_neg(&mut a.im, &ai);
        let tt = *t;
        fp2_mul(t, &tt, &a);
    } else {
        fp2_copy(&mut q.x, &p.x);
        fp2_sqr(t, &p.z);
        let tt = *t;
        fp2_sqr(t, &tt);
    }
    fp2_copy(&mut q.y, &p.y);
    fp2_copy(&mut q.z, &p.z);
}

/// Mirrors `jac_from_ws`.
pub fn jac_from_ws(q: &mut JacPoint, p: &JacPoint, ao3: &Fp2, curve: &EcCurve) {
    let mut t = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    if fp2_is_zero(&curve.A) != 0xFFFFFFFF {
        fp2_sqr(&mut t, &p.z);
        let tt = t;
        fp2_mul(&mut t, &tt, ao3);
        fp2_sub(&mut q.x, &p.x, &t);
    }
    fp2_copy(&mut q.y, &p.y);
    fp2_copy(&mut q.z, &p.z);
}

/// Mirrors `jac_neg`.
pub fn jac_neg(q: &mut JacPoint, p: &JacPoint) {
    fp2_copy(&mut q.x, &p.x);
    fp2_neg(&mut q.y, &p.y);
    fp2_copy(&mut q.z, &p.z);
}

/// Mirrors `DBL`: doubling on a Montgomery curve in Jacobian
/// coordinates.
pub fn dbl(q: &mut JacPoint, p: &JacPoint, ac: &EcCurve) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;

    let flag = fp2_is_zero(&p.x) & fp2_is_zero(&p.z);

    fp2_sqr(&mut t0, &p.x);
    fp2_add(&mut t1, &t0, &t0);
    let tt = t0;
    fp2_add(&mut t0, &tt, &t1);
    fp2_sqr(&mut t1, &p.z);
    fp2_mul(&mut t2, &p.x, &ac.A);
    let tt = t2;
    fp2_add(&mut t2, &tt, &tt);
    let tt = t2;
    fp2_add(&mut t2, &t1, &tt);
    let tt = t2;
    fp2_mul(&mut t2, &t1, &tt);
    let tt = t2;
    fp2_add(&mut t2, &t0, &tt);
    fp2_mul(&mut q.z, &p.y, &p.z);
    let qz = q.z;
    fp2_add(&mut q.z, &qz, &qz);
    fp2_sqr(&mut t0, &q.z);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &ac.A);
    fp2_sqr(&mut t1, &p.y);
    let tt = t1;
    fp2_add(&mut t1, &tt, &tt);
    fp2_add(&mut t3, &p.x, &p.x);
    let tt = t3;
    fp2_mul(&mut t3, &t1, &tt);
    fp2_sqr(&mut q.x, &t2);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &t0);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &t3);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &t3);
    fp2_sub(&mut q.y, &t3, &q.x);
    let qy = q.y;
    fp2_mul(&mut q.y, &qy, &t2);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    let qy = q.y;
    fp2_sub(&mut q.y, &qy, &t1);
    let qy = q.y;
    fp2_sub(&mut q.y, &qy, &t1);

    let qx = q.x;
    let qz = q.z;
    fp2_select(&mut q.x, &qx, &p.x, flag.wrapping_neg());
    fp2_select(&mut q.z, &qz, &p.z, flag.wrapping_neg());
}

/// Mirrors `DBLW`: doubling on a Weierstrass curve in modified
/// Jacobian coordinates `(X : Y : Z : T = a*Z^4)`.
pub fn dblw(q: &mut JacPoint, u: &mut Fp2, p: &JacPoint, t: &Fp2) {
    let flag = fp2_is_zero(&p.x) & fp2_is_zero(&p.z);

    let mut xx = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut c = xx;
    let mut cc = xx;
    let mut r = xx;
    let mut s = xx;
    let mut m = xx;

    fp2_sqr(&mut xx, &p.x);
    fp2_sqr(&mut c, &p.y);
    let cv = c;
    fp2_add(&mut c, &cv, &cv);
    fp2_sqr(&mut cc, &c);
    fp2_add(&mut r, &cc, &cc);
    fp2_add(&mut s, &p.x, &c);
    let sv = s;
    fp2_sqr(&mut s, &sv);
    let sv = s;
    fp2_sub(&mut s, &sv, &xx);
    let sv = s;
    fp2_sub(&mut s, &sv, &cc);
    fp2_add(&mut m, &xx, &xx);
    let mv = m;
    fp2_add(&mut m, &mv, &xx);
    let mv = m;
    fp2_add(&mut m, &mv, t);
    fp2_sqr(&mut q.x, &m);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &s);
    let qx = q.x;
    fp2_sub(&mut q.x, &qx, &s);
    fp2_mul(&mut q.z, &p.y, &p.z);
    let qz = q.z;
    fp2_add(&mut q.z, &qz, &qz);
    fp2_sub(&mut q.y, &s, &q.x);
    let qy = q.y;
    fp2_mul(&mut q.y, &qy, &m);
    let qy = q.y;
    fp2_sub(&mut q.y, &qy, &r);
    fp2_mul(u, t, &r);
    let uv = *u;
    fp2_add(u, &uv, &uv);

    let qx = q.x;
    let qz = q.z;
    fp2_select(&mut q.x, &qx, &p.x, flag.wrapping_neg());
    fp2_select(&mut q.z, &qz, &p.z, flag.wrapping_neg());
}

/// Constant-time conditional select for Jacobian points. Mirrors
/// `select_jac_point`.
pub fn select_jac_point(q: &mut JacPoint, p1: &JacPoint, p2: &JacPoint, option: u64) {
    fp2_select(&mut q.x, &p1.x, &p2.x, option as u32);
    fp2_select(&mut q.y, &p1.y, &p2.y, option as u32);
    fp2_select(&mut q.z, &p1.z, &p2.z, option as u32);
}

/// Mirrors `ADD`: addition on a Montgomery curve in Jacobian
/// coordinates.
pub fn add(r: &mut JacPoint, p: &JacPoint, q: &JacPoint, ac: &EcCurve) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;
    let mut u1 = t0;
    let mut u2 = t0;
    let mut v1 = t0;
    let mut dx = t0;
    let mut dy = t0;

    let ctl1 = fp2_is_zero(&p.z);
    let ctl2 = fp2_is_zero(&q.z);

    fp2_sqr(&mut t0, &p.z);
    fp2_sqr(&mut t1, &q.z);

    fp2_mul(&mut v1, &t1, &q.z);
    fp2_mul(&mut t2, &t0, &p.z);
    let vv = v1;
    fp2_mul(&mut v1, &vv, &p.y);
    let tt = t2;
    fp2_mul(&mut t2, &tt, &q.y);
    fp2_sub(&mut dy, &t2, &v1);
    fp2_mul(&mut u2, &t0, &q.x);
    fp2_mul(&mut u1, &t1, &p.x);
    fp2_sub(&mut dx, &u2, &u1);

    fp2_add(&mut t1, &p.y, &p.y);
    fp2_add(&mut t2, &ac.A, &ac.A);
    let tt = t2;
    fp2_mul(&mut t2, &tt, &p.x);
    let tt = t2;
    fp2_add(&mut t2, &tt, &t0);
    let tt = t2;
    fp2_mul(&mut t2, &tt, &t0);
    fp2_sqr(&mut t0, &p.x);
    let tt = t2;
    fp2_add(&mut t2, &tt, &t0);
    let tt = t2;
    fp2_add(&mut t2, &tt, &t0);
    let tt = t2;
    fp2_add(&mut t2, &tt, &t0);
    let tt = t2;
    fp2_mul(&mut t2, &tt, &q.z);

    let ctl = fp2_is_zero(&dx) & fp2_is_zero(&dy);
    let dxv = dx;
    fp2_select(&mut dx, &dxv, &t1, ctl);
    let dyv = dy;
    fp2_select(&mut dy, &dyv, &t2, ctl);

    fp2_mul(&mut t0, &p.z, &q.z);
    fp2_sqr(&mut t1, &t0);
    fp2_sqr(&mut t2, &dx);
    fp2_sqr(&mut t3, &dy);

    fp2_mul(&mut r.x, &ac.A, &t1);
    let rx = r.x;
    fp2_add(&mut r.x, &rx, &u1);
    let rx = r.x;
    fp2_add(&mut r.x, &rx, &u2);
    let rx = r.x;
    fp2_mul(&mut r.x, &rx, &t2);
    let rx = r.x;
    fp2_sub(&mut r.x, &t3, &rx);

    fp2_mul(&mut r.y, &u1, &t2);
    let ry = r.y;
    fp2_sub(&mut r.y, &ry, &r.x);
    let ry = r.y;
    fp2_mul(&mut r.y, &ry, &dy);
    let tt = t2;
    fp2_mul(&mut t3, &tt, &dx);
    let tt = t3;
    fp2_mul(&mut t3, &tt, &v1);
    let ry = r.y;
    fp2_sub(&mut r.y, &ry, &t3);

    fp2_mul(&mut r.z, &dx, &t0);

    let r_in = *r;
    select_jac_point(r, &r_in, q, ctl1 as u64);
    let r_in = *r;
    select_jac_point(r, &r_in, p, ctl2 as u64);
}

/// Mirrors `jac_to_xz_add_components`.
pub fn jac_to_xz_add_components(
    add_comp: &mut AddComponents,
    p: &JacPoint,
    q: &JacPoint,
    ac: &EcCurve,
) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;
    let mut t4 = t0;
    let mut t5 = t0;
    let mut t6 = t0;

    fp2_sqr(&mut t0, &p.z);
    fp2_sqr(&mut t1, &q.z);
    fp2_mul(&mut t2, &p.x, &t1);
    fp2_mul(&mut t3, &t0, &q.x);
    fp2_mul(&mut t4, &p.y, &q.z);
    let tt = t4;
    fp2_mul(&mut t4, &tt, &t1);
    fp2_mul(&mut t5, &p.z, &q.y);
    let tt = t5;
    fp2_mul(&mut t5, &tt, &t0);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t1);
    fp2_mul(&mut t6, &t4, &t5);
    fp2_add(&mut add_comp.v, &t6, &t6);
    let tt = t4;
    fp2_sqr(&mut t4, &tt);
    let tt = t5;
    fp2_sqr(&mut t5, &tt);
    let tt = t4;
    fp2_add(&mut t4, &tt, &t5);
    fp2_add(&mut t5, &t2, &t3);
    fp2_add(&mut t6, &t3, &t3);
    let tt = t6;
    fp2_sub(&mut t6, &t5, &tt);
    let tt = t6;
    fp2_sqr(&mut t6, &tt);
    fp2_mul(&mut t1, &ac.A, &t0);
    let tt = t1;
    fp2_add(&mut t1, &t5, &tt);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &t6);
    fp2_sub(&mut add_comp.u, &t4, &t1);
    fp2_mul(&mut add_comp.w, &t6, &t0);
}

// ---------------------------------------------------------------------------
// xeval.c: degree-2 and degree-4 isogeny evaluation
// ---------------------------------------------------------------------------

/// Degree-2 isogeny evaluation. Mirrors `xeval_2`.
pub fn xeval_2(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps2) {
    debug_assert!(r.len() >= q.len());
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    for j in 0..q.len() {
        fp2_add(&mut t0, &q[j].x, &q[j].z);
        fp2_sub(&mut t1, &q[j].x, &q[j].z);
        fp2_mul(&mut t2, &kps.K.x, &t1);
        fp2_mul(&mut t1, &kps.K.z, &t0);
        fp2_add(&mut t0, &t2, &t1);
        let tt = t1;
        fp2_sub(&mut t1, &t2, &tt);
        fp2_mul(&mut r[j].x, &q[j].x, &t0);
        fp2_mul(&mut r[j].z, &q[j].z, &t1);
    }
}

/// Degree-2 isogeny evaluation with kernel `(0 : 0)`. Mirrors
/// `xeval_2_singular`.
pub fn xeval_2_singular(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps2) {
    debug_assert!(r.len() >= q.len());
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    for i in 0..q.len() {
        fp2_mul(&mut t0, &q[i].x, &q[i].z);
        fp2_mul(&mut t1, &kps.K.x, &q[i].z);
        let tt = t1;
        fp2_add(&mut t1, &tt, &q[i].x);
        let tt = t1;
        fp2_mul(&mut t1, &tt, &q[i].x);
        fp2_sqr(&mut r[i].x, &q[i].z);
        let rx = r[i].x;
        fp2_add(&mut r[i].x, &rx, &t1);
        fp2_mul(&mut r[i].z, &t0, &kps.K.z);
    }
}

/// Degree-4 isogeny evaluation. Mirrors `xeval_4`.
pub fn xeval_4(r: &mut [EcPoint], q: &[EcPoint], kps: &EcKps4) {
    debug_assert!(r.len() >= q.len());
    let k = &kps.K;
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;

    for i in 0..q.len() {
        fp2_add(&mut t0, &q[i].x, &q[i].z);
        fp2_sub(&mut t1, &q[i].x, &q[i].z);
        fp2_mul(&mut r[i].x, &t0, &k[1].x);
        fp2_mul(&mut r[i].z, &t1, &k[2].x);
        let tt = t0;
        fp2_mul(&mut t0, &tt, &t1);
        let tt = t0;
        fp2_mul(&mut t0, &tt, &k[0].x);
        fp2_add(&mut t1, &r[i].x, &r[i].z);
        let rx = r[i].x;
        let rz = r[i].z;
        fp2_sub(&mut r[i].z, &rx, &rz);
        let tt = t1;
        fp2_sqr(&mut t1, &tt);
        let rz = r[i].z;
        fp2_sqr(&mut r[i].z, &rz);
        fp2_add(&mut r[i].x, &t0, &t1);
        let tt = t0;
        fp2_sub(&mut t0, &tt, &r[i].z);
        let rx = r[i].x;
        fp2_mul(&mut r[i].x, &rx, &t1);
        let rz = r[i].z;
        fp2_mul(&mut r[i].z, &rz, &t0);
    }
}

// ---------------------------------------------------------------------------
// xisog.c: degree-2 and degree-4 isogeny construction
// ---------------------------------------------------------------------------

/// Degree-2 isogeny construction. Mirrors `xisog_2`.
pub fn xisog_2(kps: &mut EcKps2, b: &mut EcPoint, p: EcPoint) {
    fp2_sqr(&mut b.x, &p.x);
    fp2_sqr(&mut b.z, &p.z);
    let bz = b.z;
    let bx = b.x;
    fp2_sub(&mut b.x, &bz, &bx);
    fp2_add(&mut kps.K.x, &p.x, &p.z);
    fp2_sub(&mut kps.K.z, &p.x, &p.z);
}

/// Degree-2 isogeny construction with kernel above `(0, 0)`. Mirrors
/// `xisog_2_singular`. The kernel point is identified by the input
/// coefficient `A24` (passed by value, matching the C signature
/// `ec_point_t A24`).
pub fn xisog_2_singular(kps: &mut EcKps2, b24: &mut EcPoint, mut a24: EcPoint) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut four = t0;
    fp2_set_small(&mut four, 4);
    fp2_add(&mut t0, &a24.x, &a24.x);
    let tt = t0;
    fp2_sub(&mut t0, &tt, &a24.z);
    let tt = t0;
    fp2_add(&mut t0, &tt, &tt);
    fp2_inv(&mut a24.z);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &a24.z);
    fp2_copy(&mut kps.K.x, &t0);
    fp2_add(&mut b24.x, &t0, &t0);
    let tt = t0;
    fp2_sqr(&mut t0, &tt);
    let tt = t0;
    fp2_sub(&mut t0, &tt, &four);
    fp2_sqrt(&mut t0);
    fp2_neg(&mut kps.K.z, &t0);
    fp2_add(&mut b24.z, &t0, &t0);
    let bx = b24.x;
    fp2_add(&mut b24.x, &bx, &b24.z);
    let bz = b24.z;
    fp2_add(&mut b24.z, &bz, &bz);
}

/// Degree-4 isogeny construction. Mirrors `xisog_4`.
pub fn xisog_4(kps: &mut EcKps4, b: &mut EcPoint, p: EcPoint) {
    fp2_sqr(&mut kps.K[0].x, &p.x);
    fp2_sqr(&mut kps.K[0].z, &p.z);
    let k0_x = kps.K[0].x;
    let k0_z = kps.K[0].z;
    fp2_add(&mut kps.K[1].x, &k0_z, &k0_x);
    fp2_sub(&mut kps.K[1].z, &k0_z, &k0_x);
    let k1_x = kps.K[1].x;
    let k1_z = kps.K[1].z;
    fp2_mul(&mut b.x, &k1_x, &k1_z);
    fp2_sqr(&mut b.z, &k0_z);

    // Constants for xeval_4
    fp2_add(&mut kps.K[2].x, &p.x, &p.z);
    fp2_sub(&mut kps.K[1].x, &p.x, &p.z);
    fp2_add(&mut kps.K[0].x, &k0_z, &k0_z);
    let kx = kps.K[0].x;
    fp2_add(&mut kps.K[0].x, &kx, &kx);
}

// ---------------------------------------------------------------------------
// e0_basis.c (precomp lvl1)
// ---------------------------------------------------------------------------

/// Precomputed x-coordinate of the basis point P on the curve E0 (A = 0).
/// Mirrors the `RADIX == 64`, non-Broadwell branch of
/// `const fp2_t BASIS_E0_PX` in
/// `the-sqisign/src/precomp/ref/lvl1/e0_basis.c`.
pub const BASIS_E0_PX: Fp2 = Fp2 {
    re: [
        0x5bcab12000c08,
        0x452654b56d052,
        0x26f81b5190a0a,
        0x36cfd66a361eb,
        0x12726610d11b,
    ],
    im: [
        0x6b96065c83efc,
        0x29da1d4a82cd9,
        0x190797ab98bdf,
        0x6841aa6eeee05,
        0x1377c5431166,
    ],
};

/// Precomputed x-coordinate of the basis point Q on the curve E0 (A = 0).
/// Mirrors the `RADIX == 64`, non-Broadwell branch of
/// `const fp2_t BASIS_E0_QX` in
/// `the-sqisign/src/precomp/ref/lvl1/e0_basis.c`.
pub const BASIS_E0_QX: Fp2 = Fp2 {
    re: [
        0x21dd55b97832f,
        0x210f2d30b26ad,
        0x680bcfcf6396,
        0x27b318ec126a7,
        0x4ffba5956012,
    ],
    im: [
        0x74590149117e3,
        0x4982edefcc606,
        0x2ae3db0cc6884,
        0x7d0384872f5ec,
        0x4fbb0fcb5a52,
    ],
};

// ---------------------------------------------------------------------------
// isog_chains.c: long isogeny chains and curve isomorphisms
// ---------------------------------------------------------------------------

/// Isomorphism descriptor between two Montgomery curves: the map
/// `(X : Z) |-> ((Nx X + Nz Z) : (D Z))`. Mirrors
/// `struct ec_isom_t { fp2_t Nx; fp2_t Nz; fp2_t D; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcIsom {
    pub Nx: Fp2,
    pub Nz: Fp2,
    pub D: Fp2,
}

impl EcIsom {
    pub const fn zero() -> Self {
        EcIsom {
            Nx: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            Nz: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            D: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
        }
    }
}

/// Even-degree isogeny: a domain curve, a kernel generator, and the
/// 2-isogeny walk length. Mirrors
/// `struct ec_isog_even_t { ec_curve_t curve; ec_point_t kernel; unsigned length; }`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct EcIsogEven {
    pub curve: EcCurve,
    pub kernel: EcPoint,
    pub length: u32,
}

impl EcIsogEven {
    pub const fn zero() -> Self {
        EcIsogEven {
            curve: EcCurve::zero(),
            kernel: EcPoint::zero(),
            length: 0,
        }
    }
}

/// Mirrors `ec_isomorphism`. Returns `0xFFFFFFFF` on error (`Nx` or `D`
/// zero after construction), `0` otherwise.
pub fn ec_isomorphism(isom: &mut EcIsom, from: &EcCurve, to: &EcCurve) -> u32 {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;
    let mut t4 = t0;

    fp2_mul(&mut t0, &from.A, &from.C);
    fp2_mul(&mut t1, &to.A, &to.C);

    fp2_mul(&mut t2, &t1, &to.C); // toA*toC^2
    fp2_add(&mut t3, &t2, &t2);
    let tt = t3;
    fp2_add(&mut t3, &tt, &tt);
    let tt = t3;
    fp2_add(&mut t3, &tt, &tt);
    let tt = t2;
    fp2_add(&mut t2, &tt, &t3); // 9*toA*toC^2
    fp2_sqr(&mut t3, &to.A);
    let tt = t3;
    fp2_mul(&mut t3, &tt, &to.A); // toA^3
    let tt = t3;
    fp2_add(&mut t3, &tt, &tt);
    fp2_sub(&mut isom.Nx, &t3, &t2); // 2*toA^3 - 9*toA*toC^2
    fp2_mul(&mut t2, &t0, &from.A); // fromA^2 * fromC
    fp2_sqr(&mut t3, &from.C);
    let tt = t3;
    fp2_mul(&mut t3, &tt, &from.C); // fromC^3
    fp2_add(&mut t4, &t3, &t3);
    let tt = t3;
    fp2_add(&mut t3, &t4, &tt); // 3*fromC^3
    let tt = t3;
    fp2_sub(&mut t3, &tt, &t2); // 3*fromC^3 - fromA^2*fromC
    let nx = isom.Nx;
    fp2_mul(&mut isom.Nx, &nx, &t3);

    fp2_mul(&mut t2, &t0, &from.C); // fromA*fromC^2
    fp2_add(&mut t3, &t2, &t2);
    let tt = t3;
    fp2_add(&mut t3, &tt, &tt);
    let tt = t3;
    fp2_add(&mut t3, &tt, &tt);
    let tt = t2;
    fp2_add(&mut t2, &tt, &t3); // 9*fromA*fromC^2
    fp2_sqr(&mut t3, &from.A);
    let tt = t3;
    fp2_mul(&mut t3, &tt, &from.A); // fromA^3
    let tt = t3;
    fp2_add(&mut t3, &tt, &tt);
    fp2_sub(&mut isom.D, &t3, &t2); // 2*fromA^3 - 9*fromA*fromC^2
    fp2_mul(&mut t2, &t1, &to.A); // toA^2 * toC
    fp2_sqr(&mut t3, &to.C);
    let tt = t3;
    fp2_mul(&mut t3, &tt, &to.C); // toC^3
    fp2_add(&mut t4, &t3, &t3);
    let tt = t3;
    fp2_add(&mut t3, &t4, &tt); // 3*toC^3
    let tt = t3;
    fp2_sub(&mut t3, &tt, &t2); // 3*toC^3 - toA^2*toC
    let d = isom.D;
    fp2_mul(&mut isom.D, &d, &t3);

    // Mont -> SW -> SW -> Mont
    fp2_mul(&mut t0, &to.C, &from.A);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &isom.Nx); // lambda_x * toC * fromA
    fp2_mul(&mut t1, &from.C, &to.A);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &isom.D); // lambda_z * fromC * toA
    fp2_sub(&mut isom.Nz, &t0, &t1); // lambda_x*toC*fromA - lambda_z*fromC*toA
    fp2_mul(&mut t0, &from.C, &to.C);
    fp2_add(&mut t1, &t0, &t0);
    let tt = t0;
    fp2_add(&mut t0, &tt, &t1); // 3 * fromC * toC
    let d = isom.D;
    fp2_mul(&mut isom.D, &d, &t0); // 3 * lambda_z * fromC * toC
    let nx = isom.Nx;
    fp2_mul(&mut isom.Nx, &nx, &t0); // 3 * lambda_x * fromC * toC

    fp2_is_zero(&isom.Nx) | fp2_is_zero(&isom.D)
}

/// Mirrors `ec_iso_eval`: in-place isomorphism evaluation on a point.
pub fn ec_iso_eval(p: &mut EcPoint, isom: &EcIsom) {
    let px = p.x;
    fp2_mul(&mut p.x, &px, &isom.Nx);
    let mut tmp = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_mul(&mut tmp, &p.z, &isom.Nz);
    let px = p.x;
    fp2_add(&mut p.x, &px, &tmp);
    let pz = p.z;
    fp2_mul(&mut p.z, &pz, &isom.D);
}

/// Static helper mirroring `ec_eval_even_strategy`.
///
/// Walks a 2^len isogeny chain via degree-4 steps with a degree-2 tail
/// when `len` is odd. Returns `0` on success, `0xFFFFFFFF` on error
/// (kernel of wrong order or special isogeny encountered).
fn ec_eval_even_strategy(
    curve: &mut EcCurve,
    points: &mut [EcPoint],
    kernel: &EcPoint,
    isog_len: i32,
) -> u32 {
    ec_curve_normalize_a24(curve);
    let mut a24 = EcPoint::zero();
    copy_point(&mut a24, &curve.A24);

    // space = 1 + ceil(log2(isog_len)) when isog_len > 1, else 1
    let mut space: usize = 1;
    let mut i = 1i32;
    while i < isog_len {
        space += 1;
        i *= 2;
    }

    let mut splits: Vec<EcPoint> = vec![EcPoint::zero(); space];
    let mut todo: Vec<u16> = vec![0u16; space];
    splits[0] = *kernel;
    todo[0] = isog_len as u16;

    let mut current: i32 = 0;

    // Chain of 4-isogenies
    for j in 0..(isog_len / 2) {
        debug_assert!(current >= 0);
        debug_assert!(todo[current as usize] >= 1);
        while todo[current as usize] != 2 {
            debug_assert!(todo[current as usize] >= 3);
            current += 1;
            debug_assert!((current as usize) < space);
            let src = splits[(current - 1) as usize];
            splits[current as usize] = src;

            let prev_todo = todo[(current - 1) as usize];
            let mut num_dbls: u32 = ((prev_todo as u32) / 4) * 2 + ((prev_todo as u32) % 2);
            todo[current as usize] = prev_todo - num_dbls as u16;
            while num_dbls > 0 {
                let s_in = splits[current as usize];
                x_dbl_a24(&mut splits[current as usize], &s_in, &a24, false);
                num_dbls -= 1;
            }
        }

        if j == 0 {
            debug_assert!(fp2_is_one(&a24.z) == 0xFFFFFFFF);
            // ec_is_four_torsion takes a curve, but the C code calls it with
            // the same `curve` whose A24 we copied to `a24`. Reproduce that.
            if ec_is_four_torsion(&splits[current as usize], curve) != 0xFFFFFFFF {
                return 0xFFFFFFFFu32;
            }
            let mut t = EcPoint::zero();
            x_dbl_a24(&mut t, &splits[current as usize], &a24, false);
            if fp2_is_zero(&t.x) == 0xFFFFFFFF {
                return 0xFFFFFFFFu32;
            }
        }

        // Evaluate 4-isogeny
        let mut kps4 = EcKps4::zero();
        let k = splits[current as usize];
        xisog_4(&mut kps4, &mut a24, k);
        // Evaluate over the first `current` split points
        if current > 0 {
            let mut tmp: Vec<EcPoint> = vec![EcPoint::zero(); current as usize];
            xeval_4(&mut tmp, &splits[..current as usize], &kps4);
            splits[..current as usize].copy_from_slice(&tmp);
        }
        for t in todo.iter_mut().take(current as usize) {
            *t = t.wrapping_sub(2);
        }
        // Evaluate over the user's points
        let mut tmp: Vec<EcPoint> = vec![EcPoint::zero(); points.len()];
        xeval_4(&mut tmp, points, &kps4);
        points.copy_from_slice(&tmp);

        current -= 1;
    }
    debug_assert!(if isog_len % 2 != 0 {
        current == 0
    } else {
        current == -1
    });

    // Final 2-isogeny
    if isog_len % 2 != 0 {
        if isog_len == 1 && ec_is_two_torsion(&splits[0], curve) != 0xFFFFFFFF {
            return 0xFFFFFFFFu32;
        }
        if fp2_is_zero(&splits[0].x) == 0xFFFFFFFF {
            return 0xFFFFFFFFu32;
        }

        let mut kps2 = EcKps2::zero();
        let k = splits[0];
        xisog_2(&mut kps2, &mut a24, k);
        let mut tmp: Vec<EcPoint> = vec![EcPoint::zero(); points.len()];
        xeval_2(&mut tmp, points, &kps2);
        points.copy_from_slice(&tmp);
    }

    // Output curve in the form (A : C)
    a24_to_ac(curve, &a24);
    curve.is_A24_computed_and_normalized = false;
    0
}

/// Mirrors `ec_eval_even`. Returns `0` on success or `0xFFFFFFFF` on
/// error.
pub fn ec_eval_even(image: &mut EcCurve, phi: &EcIsogEven, points: &mut [EcPoint]) -> u32 {
    copy_curve(image, &phi.curve);
    ec_eval_even_strategy(image, points, &phi.kernel, phi.length as i32)
}

/// Mirrors `ec_eval_small_chain`. Returns `0` on success or `0xFFFFFFFF`
/// on error.
pub fn ec_eval_small_chain(
    curve: &mut EcCurve,
    kernel: &EcPoint,
    len: i32,
    points: &mut [EcPoint],
    special: bool,
) -> u32 {
    let mut a24 = EcPoint::zero();
    ac_to_a24(&mut a24, curve);

    let mut kps = EcKps2::zero();
    let mut small_k = EcPoint::zero();
    let mut big_k = EcPoint::zero();
    copy_point(&mut big_k, kernel);

    for i in 0..len {
        copy_point(&mut small_k, &big_k);
        for _ in 0..(len - i - 1) {
            let s_in = small_k;
            x_dbl_a24(&mut small_k, &s_in, &a24, false);
        }
        if i == 0 && ec_is_two_torsion(&small_k, curve) != 0xFFFFFFFF {
            return 0xFFFFFFFFu32;
        }
        if fp2_is_zero(&small_k.x) == 0xFFFFFFFF {
            if special {
                let mut b24 = EcPoint::zero();
                xisog_2_singular(&mut kps, &mut b24, a24);
                let big_in = [big_k];
                let mut big_out = [EcPoint::zero(); 1];
                xeval_2_singular(&mut big_out, &big_in, &kps);
                big_k = big_out[0];
                let mut tmp: Vec<EcPoint> = vec![EcPoint::zero(); points.len()];
                xeval_2_singular(&mut tmp, points, &kps);
                points.copy_from_slice(&tmp);
                copy_point(&mut a24, &b24);
            } else {
                return 0xFFFFFFFFu32;
            }
        } else {
            xisog_2(&mut kps, &mut a24, small_k);
            let big_in = [big_k];
            let mut big_out = [EcPoint::zero(); 1];
            xeval_2(&mut big_out, &big_in, &kps);
            big_k = big_out[0];
            let mut tmp: Vec<EcPoint> = vec![EcPoint::zero(); points.len()];
            xeval_2(&mut tmp, points, &kps);
            points.copy_from_slice(&tmp);
        }
    }
    a24_to_ac(curve, &a24);
    curve.is_A24_computed_and_normalized = false;
    0
}

// ---------------------------------------------------------------------------
// basis.c: torsion basis generation, x-only point lifting, y-recovery
// ---------------------------------------------------------------------------

/// Mirrors `ec_recover_y`. Returns `0xFFFFFFFF` if the point lies on the
/// curve, `0` otherwise.
pub fn ec_recover_y(y: &mut Fp2, px: &Fp2, curve: &EcCurve) -> u32 {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_sqr(&mut t0, px);
    fp2_mul(y, &t0, &curve.A); // Ax^2
    let yv = *y;
    fp2_add(y, &yv, px); // Ax^2 + x
    let tt = t0;
    fp2_mul(&mut t0, &tt, px);
    let yv = *y;
    fp2_add(y, &yv, &t0); // x^3 + Ax^2 + x
    fp2_sqrt_verify(y)
}

/// Static helper mirroring `difference_point` in basis.c.
fn difference_point(pq: &mut EcPoint, p: &EcPoint, q: &EcPoint, curve: &EcCurve) {
    let mut bxx = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut bxz = bxx;
    let mut bzz = bxx;
    let mut t0 = bxx;
    let mut t1 = bxx;

    fp2_mul(&mut t0, &p.x, &q.x);
    fp2_mul(&mut t1, &p.z, &q.z);
    fp2_sub(&mut bxx, &t0, &t1);
    let bv = bxx;
    fp2_sqr(&mut bxx, &bv);
    let bv = bxx;
    fp2_mul(&mut bxx, &bv, &curve.C); // C*(P.x*Q.x - P.z*Q.z)^2
    fp2_add(&mut bxz, &t0, &t1);
    fp2_mul(&mut t0, &p.x, &q.z);
    fp2_mul(&mut t1, &p.z, &q.x);
    fp2_add(&mut bzz, &t0, &t1);
    let bv = bxz;
    fp2_mul(&mut bxz, &bv, &bzz);
    fp2_sub(&mut bzz, &t0, &t1);
    let bv = bzz;
    fp2_sqr(&mut bzz, &bv);
    let bv = bzz;
    fp2_mul(&mut bzz, &bv, &curve.C); // C*(P.x*Q.z - P.z*Q.x)^2
    let bv = bxz;
    fp2_mul(&mut bxz, &bv, &curve.C); // C*(P.x*Q.x + P.z*Q.z)(P.x*Q.z + P.z*Q.x)
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t1);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &curve.A);
    let tt = t0;
    fp2_add(&mut t0, &tt, &tt);
    let bv = bxz;
    fp2_add(&mut bxz, &bv, &t0);

    // Normalize by C*C_bar^2*(P.z)_bar^2*(Q.z)_bar^2
    fp_copy(&mut t0.re, &curve.C.re);
    fp_neg(&mut t0.im, &curve.C.im);
    let tt = t0;
    fp2_sqr(&mut t0, &tt);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &curve.C);
    fp_copy(&mut t1.re, &p.z.re);
    fp_neg(&mut t1.im, &p.z.im);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t1);
    fp_copy(&mut t1.re, &q.z.re);
    fp_neg(&mut t1.im, &q.z.im);
    let tt = t1;
    fp2_sqr(&mut t1, &tt);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t1);
    let bv = bxx;
    fp2_mul(&mut bxx, &bv, &t0);
    let bv = bxz;
    fp2_mul(&mut bxz, &bv, &t0);
    let bv = bzz;
    fp2_mul(&mut bzz, &bv, &t0);

    // Solve the quadratic
    fp2_sqr(&mut t0, &bxz);
    fp2_mul(&mut t1, &bxx, &bzz);
    let tt = t0;
    fp2_sub(&mut t0, &tt, &t1);
    fp2_sqrt(&mut t0);
    fp2_add(&mut pq.x, &bxz, &t0);
    fp2_copy(&mut pq.z, &bzz);
}

/// Mirrors `lift_basis_normalized`. Assumes `B->P.z == 1` and
/// `E->C == 1`.
pub fn lift_basis_normalized(p: &mut JacPoint, q: &mut JacPoint, b: &EcBasis, e: &EcCurve) -> u32 {
    debug_assert!(fp2_is_one(&b.P.z) == 0xFFFFFFFF);
    debug_assert!(fp2_is_one(&e.C) == 0xFFFFFFFF);

    fp2_copy(&mut p.x, &b.P.x);
    fp2_copy(&mut q.x, &b.Q.x);
    fp2_copy(&mut q.z, &b.Q.z);
    fp2_set_one(&mut p.z);
    let ret = ec_recover_y(&mut p.y, &p.x, e);

    // Okeya-Sakurai recovery for Q.
    let mut v1 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut v2 = v1;
    let mut v3 = v1;
    let mut v4 = v1;
    fp2_mul(&mut v1, &p.x, &q.z);
    fp2_add(&mut v2, &q.x, &v1);
    fp2_sub(&mut v3, &q.x, &v1);
    let vv = v3;
    fp2_sqr(&mut v3, &vv);
    let vv = v3;
    fp2_mul(&mut v3, &vv, &b.PmQ.x);
    fp2_add(&mut v1, &e.A, &e.A);
    let vv = v1;
    fp2_mul(&mut v1, &vv, &q.z);
    let vv = v2;
    fp2_add(&mut v2, &vv, &v1);
    fp2_mul(&mut v4, &p.x, &q.x);
    let vv = v4;
    fp2_add(&mut v4, &vv, &q.z);
    let vv = v2;
    fp2_mul(&mut v2, &vv, &v4);
    let vv = v1;
    fp2_mul(&mut v1, &vv, &q.z);
    let vv = v2;
    fp2_sub(&mut v2, &vv, &v1);
    let vv = v2;
    fp2_mul(&mut v2, &vv, &b.PmQ.z);
    fp2_sub(&mut q.y, &v3, &v2);
    fp2_add(&mut v1, &p.y, &p.y);
    let vv = v1;
    fp2_mul(&mut v1, &vv, &q.z);
    let vv = v1;
    fp2_mul(&mut v1, &vv, &b.PmQ.z);
    let qx = q.x;
    fp2_mul(&mut q.x, &qx, &v1);
    let qz = q.z;
    fp2_mul(&mut q.z, &qz, &v1);

    // Transform to Jacobian coordinates.
    fp2_sqr(&mut v1, &q.z);
    let qy = q.y;
    fp2_mul(&mut q.y, &qy, &v1);
    let qx = q.x;
    fp2_mul(&mut q.x, &qx, &q.z);
    ret
}

/// Mirrors `lift_basis`. Normalizes the curve and `B->P.z` in place,
/// then calls [`lift_basis_normalized`].
pub fn lift_basis(p: &mut JacPoint, q: &mut JacPoint, b: &mut EcBasis, e: &mut EcCurve) -> u32 {
    let mut inverses: [Fp2; 2] = [
        Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        },
        Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        },
    ];
    fp2_copy(&mut inverses[0], &b.P.z);
    fp2_copy(&mut inverses[1], &e.C);

    fp2_batched_inv(&mut inverses);
    fp2_set_one(&mut b.P.z);
    fp2_set_one(&mut e.C);

    let bx = b.P.x;
    fp2_mul(&mut b.P.x, &bx, &inverses[0]);
    let ea = e.A;
    fp2_mul(&mut e.A, &ea, &inverses[1]);

    lift_basis_normalized(p, q, b, e)
}

/// Static helper mirroring `is_on_curve` in basis.c. Assumes `C == 1`.
fn is_on_curve(x: &Fp2, curve: &EcCurve) -> u32 {
    debug_assert!(fp2_is_one(&curve.C) == 0xFFFFFFFF);
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_add(&mut t0, x, &curve.A);
    let tt = t0;
    fp2_mul(&mut t0, &tt, x);
    let tt = t0;
    fp2_add_one(&mut t0, &tt);
    let tt = t0;
    fp2_mul(&mut t0, &tt, x);
    fp2_is_square(&t0)
}

/// Static helper mirroring `clear_cofactor_for_maximal_even_order` in
/// basis.c.
fn clear_cofactor_for_maximal_even_order(p: &mut EcPoint, curve: &mut EcCurve, f: i32) {
    let p_in = *p;
    ec_mul(
        p,
        &P_COFACTOR_FOR_2F,
        P_COFACTOR_FOR_2F_BITLENGTH as i32,
        &p_in,
        curve,
    );

    let lim = TORSION_EVEN_POWER as i32 - f;
    for _ in 0..lim {
        let p_in = *p;
        x_dbl_a24(p, &p_in, &curve.A24, curve.is_A24_computed_and_normalized);
    }
}

/// Static helper mirroring `find_nqr_factor` in basis.c.
fn find_nqr_factor(x: &mut Fp2, curve: &EcCurve, start: u8) -> u8 {
    let mut n: u16 = start as u16;
    let mut qr_b: bool = true;
    let mut found: u32 = 0;
    let mut b: Fp = [0u64; NWORDS_FIELD];
    let mut z = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t0 = z;
    let mut t1 = z;

    while found == 0 {
        while qr_b {
            let mut tmp: Fp = [0u64; NWORDS_FIELD];
            let val: u32 = (n as u32) * (n as u32) + 1;
            fp_set_small(&mut tmp, val as u64);
            qr_b = fp_is_square(&tmp) != 0;
            n = n.wrapping_add(1);
        }

        let val: u32 = (n as u32) - 1;
        fp_set_small(&mut b, val as u64);
        fp2_set_zero(&mut t0);
        fp2_set_one(&mut z);
        fp_copy(&mut z.im, &b);
        fp_copy(&mut t0.im, &b);

        fp2_sqr(&mut t1, &curve.A);
        let tt = t0;
        fp2_mul(&mut t0, &tt, &t1); // A^2 * (z - 1)
        fp2_sqr(&mut t1, &z);
        let tt = t0;
        fp2_sub(&mut t0, &tt, &t1); // A^2 * (z - 1) - z^2
        found = if fp2_is_square(&t0) == 0 { 1 } else { 0 };

        qr_b = true;
    }

    // x = -A / (1 + i*b)
    fp2_copy(x, &z);
    fp2_inv(x);
    let xv = *x;
    fp2_mul(x, &xv, &curve.A);
    let xv = *x;
    fp2_neg(x, &xv);

    if n <= 128 {
        (n - 1) as u8
    } else {
        0
    }
}

/// Static helper mirroring `find_nA_x_coord` in basis.c.
fn find_nA_x_coord(x: &mut Fp2, curve: &EcCurve, start: u8) -> u8 {
    debug_assert!(fp2_is_square(&curve.A) == 0);

    let mut n: u8 = start;
    if n == 1 {
        fp2_copy(x, &curve.A);
    } else {
        fp2_mul_small(x, &curve.A, n as u32);
    }
    while is_on_curve(x, curve) == 0 {
        let xv = *x;
        fp2_add(x, &xv, &curve.A);
        n = n.wrapping_add(1);
    }
    if n < 128 {
        n
    } else {
        0
    }
}

/// Static helper mirroring `ec_basis_E0_2f` in basis.c.
fn ec_basis_e0_2f(pq2: &mut EcBasis, curve: &EcCurve, f: i32) {
    debug_assert!(fp2_is_zero(&curve.A) == 0xFFFFFFFF);
    let mut p = EcPoint::zero();
    let mut q = EcPoint::zero();

    fp2_copy(&mut p.x, &BASIS_E0_PX);
    fp2_copy(&mut q.x, &BASIS_E0_QX);
    fp2_set_one(&mut p.z);
    fp2_set_one(&mut q.z);

    let lim = TORSION_EVEN_POWER as i32 - f;
    for _ in 0..lim {
        let p_in = p;
        x_dbl_e0(&mut p, &p_in);
        let q_in = q;
        x_dbl_e0(&mut q, &q_in);
    }

    copy_point(&mut pq2.P, &p);
    copy_point(&mut pq2.Q, &q);
    difference_point(&mut pq2.PmQ, &p, &q, curve);
}

/// Mirrors `ec_curve_to_basis_2f_to_hint`. Returns the recomputation
/// hint.
pub fn ec_curve_to_basis_2f_to_hint(pq2: &mut EcBasis, curve: &mut EcCurve, f: i32) -> u8 {
    ec_normalize_curve_and_a24(curve);

    if fp2_is_zero(&curve.A) == 0xFFFFFFFF {
        ec_basis_e0_2f(pq2, curve, f);
        return 0;
    }

    let hint_a: bool = fp2_is_square(&curve.A) != 0;
    let mut p = EcPoint::zero();
    let mut q = EcPoint::zero();

    let hint = if !hint_a {
        find_nA_x_coord(&mut p.x, curve, 1)
    } else {
        find_nqr_factor(&mut p.x, curve, 1)
    };

    fp2_set_one(&mut p.z);
    fp2_add(&mut q.x, &curve.A, &p.x);
    let qx = q.x;
    fp2_neg(&mut q.x, &qx);
    fp2_set_one(&mut q.z);

    clear_cofactor_for_maximal_even_order(&mut p, curve, f);
    clear_cofactor_for_maximal_even_order(&mut q, curve, f);

    difference_point(&mut pq2.Q, &p, &q, curve);
    copy_point(&mut pq2.P, &p);
    copy_point(&mut pq2.PmQ, &q);

    debug_assert!(hint < 128);
    (hint << 1) | (hint_a as u8)
}

/// Mirrors `ec_curve_to_basis_2f_from_hint`. Returns `1` if the basis is
/// valid, `0` otherwise.
pub fn ec_curve_to_basis_2f_from_hint(
    pq2: &mut EcBasis,
    curve: &mut EcCurve,
    f: i32,
    hint: u8,
) -> i32 {
    ec_normalize_curve_and_a24(curve);

    if fp2_is_zero(&curve.A) == 0xFFFFFFFF {
        ec_basis_e0_2f(pq2, curve, f);
        return 1;
    }

    let hint_a: bool = (hint & 1) != 0;
    let hint_p: u8 = hint >> 1;

    let mut p = EcPoint::zero();
    let mut q = EcPoint::zero();

    if hint_p == 0 {
        if !hint_a {
            find_nA_x_coord(&mut p.x, curve, 128);
        } else {
            find_nqr_factor(&mut p.x, curve, 128);
        }
    } else if !hint_a {
        fp2_mul_small(&mut p.x, &curve.A, hint_p as u32);
    } else {
        fp_set_one(&mut p.x.re);
        fp_set_small(&mut p.x.im, hint_p as u64);
        fp2_inv(&mut p.x);
        let px = p.x;
        fp2_mul(&mut p.x, &px, &curve.A);
        let px = p.x;
        fp2_neg(&mut p.x, &px);
    }
    fp2_set_one(&mut p.z);

    fp2_add(&mut q.x, &curve.A, &p.x);
    let qx = q.x;
    fp2_neg(&mut q.x, &qx);
    fp2_set_one(&mut q.z);

    clear_cofactor_for_maximal_even_order(&mut p, curve, f);
    clear_cofactor_for_maximal_even_order(&mut q, curve, f);

    difference_point(&mut pq2.Q, &p, &q, curve);
    copy_point(&mut pq2.P, &p);
    copy_point(&mut pq2.PmQ, &q);

    1
}

// ---------------------------------------------------------------------------
// biextension.c: Weil and Tate pairings via the cubical-torsor biextension
// ladder, plus E[2^e] dlogs derived from them.
// ---------------------------------------------------------------------------

/// Pairing parameters for the Weil / reduced Tate routines. Mirrors
/// `struct pairing_params { uint32_t e; ec_point_t P; ec_point_t Q;
/// ec_point_t PQ; fp2_t ixP; fp2_t ixQ; ec_point_t A24; }`. The struct
/// is internal to the upstream `biextension.c`; the port keeps it
/// internal too (no public callers).
#[derive(Clone, Copy, Debug)]
struct PairingParams {
    e: u32,
    P: EcPoint,
    Q: EcPoint,
    PQ: EcPoint,
    ixP: Fp2,
    ixQ: Fp2,
    A24: EcPoint,
}

impl PairingParams {
    const fn zero() -> Self {
        PairingParams {
            e: 0,
            P: EcPoint::zero(),
            Q: EcPoint::zero(),
            PQ: EcPoint::zero(),
            ixP: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            ixQ: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            A24: EcPoint::zero(),
        }
    }
}

/// Storage for the four x-only difference points needed by the
/// dlog pairing routines. Mirrors `struct pairing_dlog_diff_points`.
#[derive(Clone, Copy, Debug)]
struct PairingDlogDiffPoints {
    PmR: EcPoint,
    PmS: EcPoint,
    RmQ: EcPoint,
    SmQ: EcPoint,
}

impl PairingDlogDiffPoints {
    const fn zero() -> Self {
        PairingDlogDiffPoints {
            PmR: EcPoint::zero(),
            PmS: EcPoint::zero(),
            RmQ: EcPoint::zero(),
            SmQ: EcPoint::zero(),
        }
    }
}

/// Parameters for the dlog pairing routines. Mirrors
/// `struct pairing_dlog_params`.
#[derive(Clone, Copy, Debug)]
struct PairingDlogParams {
    e: u32,
    PQ: EcBasis,
    RS: EcBasis,
    diff: PairingDlogDiffPoints,
    ixP: Fp2,
    ixQ: Fp2,
    ixR: Fp2,
    ixS: Fp2,
    A24: EcPoint,
}

impl PairingDlogParams {
    const fn zero() -> Self {
        PairingDlogParams {
            e: 0,
            PQ: EcBasis::zero(),
            RS: EcBasis::zero(),
            diff: PairingDlogDiffPoints::zero(),
            ixP: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            ixQ: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            ixR: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            ixS: Fp2 {
                re: [0u64; NWORDS_FIELD],
                im: [0u64; NWORDS_FIELD],
            },
            A24: EcPoint::zero(),
        }
    }
}

/// Static helper mirroring `cubicalADD` in biextension.c.
fn cubical_add(r: &mut EcPoint, p: &EcPoint, q: &EcPoint, ix_pq: &Fp2) {
    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_add(&mut t2, &q.x, &q.z);
    fp2_sub(&mut t3, &q.x, &q.z);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t3);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &t2);
    fp2_add(&mut t2, &t0, &t1);
    fp2_sub(&mut t3, &t0, &t1);
    fp2_sqr(&mut r.z, &t3);
    let tt = t2;
    fp2_sqr(&mut t2, &tt);
    fp2_mul(&mut r.x, ix_pq, &t2);
}

/// Static helper mirroring `cubicalDBLADD` in biextension.c.
fn cubical_dbl_add(
    p_pq: &mut EcPoint,
    qq: &mut EcPoint,
    p: &EcPoint,
    q: &EcPoint,
    ix_pq: &Fp2,
    a24: &EcPoint,
) {
    debug_assert!(fp2_is_one(&a24.z) == 0xFFFFFFFF);

    let mut t0 = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut t1 = t0;
    let mut t2 = t0;
    let mut t3 = t0;

    fp2_add(&mut t0, &p.x, &p.z);
    fp2_sub(&mut t1, &p.x, &p.z);
    fp2_add(&mut p_pq.x, &q.x, &q.z);
    fp2_sub(&mut t3, &q.x, &q.z);
    fp2_sqr(&mut t2, &p_pq.x);
    fp2_sqr(&mut qq.z, &t3);
    let tt = t0;
    fp2_mul(&mut t0, &tt, &t3);
    let tt = t1;
    fp2_mul(&mut t1, &tt, &p_pq.x);
    fp2_add(&mut p_pq.x, &t0, &t1);
    fp2_sub(&mut t3, &t0, &t1);
    fp2_sqr(&mut p_pq.z, &t3);
    let pp = p_pq.x;
    fp2_sqr(&mut p_pq.x, &pp);
    let pp = p_pq.x;
    fp2_mul(&mut p_pq.x, ix_pq, &pp);
    fp2_sub(&mut t3, &t2, &qq.z);
    fp2_mul(&mut qq.x, &t2, &qq.z);
    fp2_mul(&mut t0, &t3, &a24.x);
    let tt = t0;
    fp2_add(&mut t0, &tt, &qq.z);
    fp2_mul(&mut qq.z, &t0, &t3);
}

/// Static helper mirroring `biext_ladder_2e` in biextension.c.
fn biext_ladder_2e(
    e: u32,
    p_nq: &mut EcPoint,
    nq: &mut EcPoint,
    pq: &EcPoint,
    q: &EcPoint,
    ix_p: &Fp2,
    a24: &EcPoint,
) {
    copy_point(p_nq, pq);
    copy_point(nq, q);
    for _ in 0..e {
        let pn_in = *p_nq;
        let nq_in = *nq;
        cubical_dbl_add(p_nq, nq, &pn_in, &nq_in, ix_p, a24);
    }
}

/// Static helper mirroring `point_ratio` in biextension.c. The C
/// reference uses `assert(ec_is_equal(...))`, which is satisfied by any
/// non-zero return value; `ec_is_equal` in turn has a known typo where
/// `*` should have been `&` in its return expression (the level-0
/// `ec.c` port records this literally), so a positive equality answer
/// can come back as `1` rather than the all-ones mask. The assertion
/// here therefore checks for non-zero rather than for `0xFFFFFFFF`.
fn point_ratio(r: &mut EcPoint, p_nq: &EcPoint, nq: &EcPoint, p: &EcPoint) {
    debug_assert!(ec_is_zero(nq) != 0);
    debug_assert!(ec_is_equal(p_nq, p) != 0);
    fp2_mul(&mut r.x, &nq.x, &p.x);
    fp2_copy(&mut r.z, &p_nq.x);
}

/// Static helper mirroring `translate` in biextension.c.
fn translate(p: &mut EcPoint, t: &EcPoint) {
    let mut px_new = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut pz_new = px_new;
    {
        let mut t0 = px_new;
        let mut t1 = px_new;
        fp2_mul(&mut t0, &t.x, &p.x);
        fp2_mul(&mut t1, &t.z, &p.z);
        fp2_sub(&mut px_new, &t0, &t1);
        fp2_mul(&mut t0, &t.z, &p.x);
        fp2_mul(&mut t1, &t.x, &p.z);
        fp2_sub(&mut pz_new, &t0, &t1);
    }
    let ta_zero = fp2_is_zero(&t.x);
    let pv = px_new;
    fp2_select(&mut px_new, &pv, &p.z, ta_zero);
    let pv = pz_new;
    fp2_select(&mut pz_new, &pv, &p.x, ta_zero);
    let tb_zero = fp2_is_zero(&t.z);
    let pv = px_new;
    fp2_select(&mut px_new, &pv, &p.x, tb_zero);
    let pv = pz_new;
    fp2_select(&mut pz_new, &pv, &p.z, tb_zero);
    fp2_copy(&mut p.x, &px_new);
    fp2_copy(&mut p.z, &pz_new);
}

/// Static helper mirroring `monodromy_i` in biextension.c.
fn monodromy_i(r: &mut EcPoint, pairing_data: &PairingParams, swap_pq: bool) {
    let mut ix_p = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut p = EcPoint::zero();
    let mut q = EcPoint::zero();
    let mut p_nq = EcPoint::zero();
    let mut nq = EcPoint::zero();

    if !swap_pq {
        copy_point(&mut p, &pairing_data.P);
        copy_point(&mut q, &pairing_data.Q);
        fp2_copy(&mut ix_p, &pairing_data.ixP);
    } else {
        copy_point(&mut p, &pairing_data.Q);
        copy_point(&mut q, &pairing_data.P);
        fp2_copy(&mut ix_p, &pairing_data.ixQ);
    }

    biext_ladder_2e(
        pairing_data.e - 1,
        &mut p_nq,
        &mut nq,
        &pairing_data.PQ,
        &q,
        &ix_p,
        &pairing_data.A24,
    );
    let nq_const = nq;
    translate(&mut p_nq, &nq_const);
    let nq_const = nq;
    translate(&mut nq, &nq_const);
    point_ratio(r, &p_nq, &nq, &p);
}

/// Static helper mirroring `cubical_normalization` in biextension.c.
fn cubical_normalization(pairing_data: &mut PairingParams, p: &EcPoint, q: &EcPoint) {
    let mut t: [Fp2; 4] = [Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    }; 4];
    fp2_copy(&mut t[0], &p.x);
    fp2_copy(&mut t[1], &p.z);
    fp2_copy(&mut t[2], &q.x);
    fp2_copy(&mut t[3], &q.z);
    fp2_batched_inv(&mut t);

    fp2_mul(&mut pairing_data.ixP, &p.z, &t[0]);
    fp2_mul(&mut pairing_data.ixQ, &q.z, &t[2]);

    fp2_mul(&mut pairing_data.P.x, &p.x, &t[1]);
    fp2_mul(&mut pairing_data.Q.x, &q.x, &t[3]);
    fp2_set_one(&mut pairing_data.P.z);
    fp2_set_one(&mut pairing_data.Q.z);
}

/// Static helper mirroring `weil_n` in biextension.c.
fn weil_n(r: &mut Fp2, pairing_data: &PairingParams) {
    let mut r0 = EcPoint::zero();
    let mut r1 = EcPoint::zero();
    monodromy_i(&mut r0, pairing_data, true);
    monodromy_i(&mut r1, pairing_data, false);

    fp2_mul(r, &r0.x, &r1.z);
    fp2_inv(r);
    let rv = *r;
    fp2_mul(r, &rv, &r0.z);
    let rv = *r;
    fp2_mul(r, &rv, &r1.x);
}

/// Mirrors `weil`. Computes `e_{2^e}(P, Q)` using the cubical biextension
/// ladder.
pub fn weil(r: &mut Fp2, e: u32, p: &EcPoint, q: &EcPoint, pq: &EcPoint, ec: &mut EcCurve) {
    let mut pairing_data = PairingParams::zero();
    pairing_data.e = e;
    cubical_normalization(&mut pairing_data, p, q);
    copy_point(&mut pairing_data.PQ, pq);

    ec_curve_normalize_a24(ec);
    copy_point(&mut pairing_data.A24, &ec.A24);

    weil_n(r, &pairing_data);
}

/// Mirrors `clear_cofac` in biextension.c. Raises `a` to the power
/// `(p_cofactor_for_2f >> 1)`.
fn clear_cofac(r: &mut Fp2, a: &Fp2) {
    let mut exp: u64 = P_COFACTOR_FOR_2F[0];
    exp >>= 1;

    let mut x = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_copy(&mut x, a);
    fp2_copy(r, a);
    while exp > 0 {
        let rv = *r;
        fp2_sqr(r, &rv);
        if exp & 1 != 0 {
            let rv = *r;
            fp2_mul(r, &rv, &x);
        }
        exp >>= 1;
    }
}

/// Mirrors `fp2_frob` in biextension.c. Applies `a + i*b -> a - i*b`.
fn fp2_frob(out: &mut Fp2, inp: &Fp2) {
    fp_copy(&mut out.re, &inp.re);
    fp_neg(&mut out.im, &inp.im);
}

/// Mirrors `reduced_tate`. Computes the reduced Tate pairing
/// `t_{2^e}(P, Q)`.
pub fn reduced_tate(r: &mut Fp2, e: u32, p: &EcPoint, q: &EcPoint, pq: &EcPoint, ec: &mut EcCurve) {
    let e_full: u32 = TORSION_EVEN_POWER as u32;
    let e_diff: u32 = e_full - e;

    let mut pairing_data = PairingParams::zero();
    pairing_data.e = e;
    cubical_normalization(&mut pairing_data, p, q);
    copy_point(&mut pairing_data.PQ, pq);

    ec_curve_normalize_a24(ec);
    copy_point(&mut pairing_data.A24, &ec.A24);

    let mut rpt = EcPoint::zero();
    monodromy_i(&mut rpt, &pairing_data, true);

    let mut frob = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    let mut tmp = frob;
    fp2_copy(&mut tmp, &rpt.x);
    fp2_frob(&mut frob, &rpt.x);
    fp2_mul(&mut rpt.x, &rpt.z, &frob);
    fp2_frob(&mut frob, &rpt.z);
    fp2_mul(&mut rpt.z, &tmp, &frob);
    fp2_inv(&mut rpt.x);
    fp2_mul(r, &rpt.x, &rpt.z);

    let rv = *r;
    clear_cofac(r, &rv);
    for _ in 0..e_diff {
        let rv = *r;
        fp2_sqr(r, &rv);
    }
}

/// Static helper mirroring `fp2_dlog_2e_rec`.
#[allow(clippy::needless_range_loop)]
fn fp2_dlog_2e_rec(
    a: &mut [u64],
    len: i64,
    pows_f: &mut [Fp2],
    pows_g: &mut [Fp2],
    stacklen: usize,
) -> bool {
    if len == 0 {
        for i in 0..NWORDS_ORDER {
            a[i] = 0;
        }
        return true;
    }
    if len == 1 {
        if fp2_is_one(&pows_f[stacklen - 1]) == 0xFFFFFFFF {
            for i in 0..NWORDS_ORDER {
                a[i] = 0;
            }
            for i in 0..stacklen - 1 {
                let gv = pows_g[i];
                fp2_sqr(&mut pows_g[i], &gv);
            }
            return true;
        } else if fp2_is_equal(&pows_f[stacklen - 1], &pows_g[stacklen - 1]) == 0xFFFFFFFF {
            a[0] = 1;
            for i in 1..NWORDS_ORDER {
                a[i] = 0;
            }
            for i in 0..stacklen - 1 {
                let fv = pows_f[i];
                fp2_mul(&mut pows_f[i], &fv, &pows_g[i]);
                let gv = pows_g[i];
                fp2_sqr(&mut pows_g[i], &gv);
            }
            return true;
        } else {
            return false;
        }
    }
    // The reference uses `(double)len * 0.5`. For integer len this is
    // equivalent to integer division by 2 (truncation toward zero); the
    // result is the same value whether computed in floating point or
    // integer arithmetic.
    let right: i64 = len / 2;
    let left: i64 = len - right;
    pows_f[stacklen] = pows_f[stacklen - 1];
    pows_g[stacklen] = pows_g[stacklen - 1];
    for _ in 0..left {
        let fv = pows_f[stacklen];
        fp2_sqr(&mut pows_f[stacklen], &fv);
        let gv = pows_g[stacklen];
        fp2_sqr(&mut pows_g[stacklen], &gv);
    }
    let mut dlp1 = [0u64; NWORDS_ORDER];
    let mut dlp2 = [0u64; NWORDS_ORDER];
    if !fp2_dlog_2e_rec(&mut dlp1, right, pows_f, pows_g, stacklen + 1) {
        return false;
    }
    if !fp2_dlog_2e_rec(&mut dlp2, left, pows_f, pows_g, stacklen) {
        return false;
    }
    // a = dlp1 + 2^right * dlp2
    multiple_mp_shiftl(&mut dlp2, right as u32);
    mp_add(&mut a[..NWORDS_ORDER], &dlp2, &dlp1);
    true
}

/// Static helper mirroring `fp2_dlog_2e`.
fn fp2_dlog_2e(scal: &mut [u64], f: &Fp2, g_inverse: &Fp2, e: i32) -> bool {
    let mut log: i64 = 0;
    let mut len: i64 = e as i64;
    while len > 1 {
        log += 1;
        len >>= 1;
    }
    log += 1;

    let mut pows_f: Vec<Fp2> = vec![
        Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD]
        };
        log as usize
    ];
    let mut pows_g: Vec<Fp2> = vec![
        Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD]
        };
        log as usize
    ];
    pows_f[0] = *f;
    pows_g[0] = *g_inverse;

    for s in scal.iter_mut().take(NWORDS_ORDER) {
        *s = 0;
    }

    let ok = fp2_dlog_2e_rec(scal, e as i64, &mut pows_f, &mut pows_g, 1);
    debug_assert!(ok);
    ok
}

/// Static helper mirroring `cubical_normalization_dlog` in
/// biextension.c.
fn cubical_normalization_dlog(pairing_dlog_data: &mut PairingDlogParams, curve: &mut EcCurve) {
    let mut t: [Fp2; 11] = [Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    }; 11];
    fp2_copy(&mut t[0], &pairing_dlog_data.PQ.P.x);
    fp2_copy(&mut t[1], &pairing_dlog_data.PQ.P.z);
    fp2_copy(&mut t[2], &pairing_dlog_data.PQ.Q.x);
    fp2_copy(&mut t[3], &pairing_dlog_data.PQ.Q.z);
    fp2_copy(&mut t[4], &pairing_dlog_data.PQ.PmQ.x);
    fp2_copy(&mut t[5], &pairing_dlog_data.PQ.PmQ.z);
    fp2_copy(&mut t[6], &pairing_dlog_data.RS.P.x);
    fp2_copy(&mut t[7], &pairing_dlog_data.RS.P.z);
    fp2_copy(&mut t[8], &pairing_dlog_data.RS.Q.x);
    fp2_copy(&mut t[9], &pairing_dlog_data.RS.Q.z);
    fp2_copy(&mut t[10], &curve.C);

    fp2_batched_inv(&mut t);

    fp2_mul(&mut pairing_dlog_data.ixP, &pairing_dlog_data.PQ.P.z, &t[0]);
    let px = pairing_dlog_data.PQ.P.x;
    fp2_mul(&mut pairing_dlog_data.PQ.P.x, &px, &t[1]);
    fp2_set_one(&mut pairing_dlog_data.PQ.P.z);

    fp2_mul(&mut pairing_dlog_data.ixQ, &pairing_dlog_data.PQ.Q.z, &t[2]);
    let qx = pairing_dlog_data.PQ.Q.x;
    fp2_mul(&mut pairing_dlog_data.PQ.Q.x, &qx, &t[3]);
    fp2_set_one(&mut pairing_dlog_data.PQ.Q.z);

    let pmqx = pairing_dlog_data.PQ.PmQ.x;
    fp2_mul(&mut pairing_dlog_data.PQ.PmQ.x, &pmqx, &t[5]);
    fp2_set_one(&mut pairing_dlog_data.PQ.PmQ.z);

    fp2_mul(&mut pairing_dlog_data.ixR, &pairing_dlog_data.RS.P.z, &t[6]);
    let rx = pairing_dlog_data.RS.P.x;
    fp2_mul(&mut pairing_dlog_data.RS.P.x, &rx, &t[7]);
    fp2_set_one(&mut pairing_dlog_data.RS.P.z);

    fp2_mul(&mut pairing_dlog_data.ixS, &pairing_dlog_data.RS.Q.z, &t[8]);
    let sx = pairing_dlog_data.RS.Q.x;
    fp2_mul(&mut pairing_dlog_data.RS.Q.x, &sx, &t[9]);
    fp2_set_one(&mut pairing_dlog_data.RS.Q.z);

    let ea = curve.A;
    fp2_mul(&mut curve.A, &ea, &t[10]);
    fp2_set_one(&mut curve.C);
}

/// Static helper mirroring `compute_difference_points` in biextension.c.
fn compute_difference_points(pairing_dlog_data: &mut PairingDlogParams, curve: &mut EcCurve) {
    let mut xy_p = JacPoint::zero();
    let mut xy_q = JacPoint::zero();
    let mut xy_r = JacPoint::zero();
    let mut xy_s = JacPoint::zero();
    let mut temp = JacPoint::zero();

    lift_basis_normalized(&mut xy_p, &mut xy_q, &pairing_dlog_data.PQ, curve);
    lift_basis_normalized(&mut xy_r, &mut xy_s, &pairing_dlog_data.RS, curve);

    // x(P - R)
    jac_neg(&mut temp, &xy_r);
    let temp_in = temp;
    add(&mut temp, &temp_in, &xy_p, curve);
    jac_to_xz(&mut pairing_dlog_data.diff.PmR, &temp);

    // x(P - S)
    jac_neg(&mut temp, &xy_s);
    let temp_in = temp;
    add(&mut temp, &temp_in, &xy_p, curve);
    jac_to_xz(&mut pairing_dlog_data.diff.PmS, &temp);

    // x(R - Q)
    jac_neg(&mut temp, &xy_q);
    let temp_in = temp;
    add(&mut temp, &temp_in, &xy_r, curve);
    jac_to_xz(&mut pairing_dlog_data.diff.RmQ, &temp);

    // x(S - Q)
    jac_neg(&mut temp, &xy_q);
    let temp_in = temp;
    add(&mut temp, &temp_in, &xy_s, curve);
    jac_to_xz(&mut pairing_dlog_data.diff.SmQ, &temp);
}

/// Static helper mirroring `weil_dlog` in biextension.c.
#[allow(clippy::needless_range_loop)]
fn weil_dlog(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    pairing_dlog_data: &PairingDlogParams,
) {
    let mut n_p = EcPoint::zero();
    let mut n_q = EcPoint::zero();
    let mut n_r = EcPoint::zero();
    let mut n_s = EcPoint::zero();
    let mut n_pq = EcPoint::zero();
    let mut p_nq = EcPoint::zero();
    let mut n_pr = EcPoint::zero();
    let mut p_nr = EcPoint::zero();
    let mut n_ps = EcPoint::zero();
    let mut p_ns = EcPoint::zero();
    let mut n_rq = EcPoint::zero();
    let mut r_nq = EcPoint::zero();
    let mut n_sq = EcPoint::zero();
    let mut s_nq = EcPoint::zero();

    copy_point(&mut n_p, &pairing_dlog_data.PQ.P);
    copy_point(&mut n_q, &pairing_dlog_data.PQ.Q);
    copy_point(&mut n_r, &pairing_dlog_data.RS.P);
    copy_point(&mut n_s, &pairing_dlog_data.RS.Q);
    copy_point(&mut n_pq, &pairing_dlog_data.PQ.PmQ);
    copy_point(&mut p_nq, &pairing_dlog_data.PQ.PmQ);
    copy_point(&mut n_pr, &pairing_dlog_data.diff.PmR);
    copy_point(&mut n_ps, &pairing_dlog_data.diff.PmS);
    copy_point(&mut p_nr, &pairing_dlog_data.diff.PmR);
    copy_point(&mut p_ns, &pairing_dlog_data.diff.PmS);
    copy_point(&mut n_rq, &pairing_dlog_data.diff.RmQ);
    copy_point(&mut n_sq, &pairing_dlog_data.diff.SmQ);
    copy_point(&mut r_nq, &pairing_dlog_data.diff.RmQ);
    copy_point(&mut s_nq, &pairing_dlog_data.diff.SmQ);

    for _ in 0..(pairing_dlog_data.e - 1) {
        let in1 = n_pq;
        cubical_add(&mut n_pq, &in1, &n_p, &pairing_dlog_data.ixQ);
        let in1 = n_pr;
        cubical_add(&mut n_pr, &in1, &n_p, &pairing_dlog_data.ixR);
        let in1 = n_ps;
        let in2 = n_p;
        cubical_dbl_add(
            &mut n_ps,
            &mut n_p,
            &in1,
            &in2,
            &pairing_dlog_data.ixS,
            &pairing_dlog_data.A24,
        );

        let in1 = p_nq;
        cubical_add(&mut p_nq, &in1, &n_q, &pairing_dlog_data.ixP);
        let in1 = r_nq;
        cubical_add(&mut r_nq, &in1, &n_q, &pairing_dlog_data.ixR);
        let in1 = s_nq;
        let in2 = n_q;
        cubical_dbl_add(
            &mut s_nq,
            &mut n_q,
            &in1,
            &in2,
            &pairing_dlog_data.ixS,
            &pairing_dlog_data.A24,
        );

        let in1 = p_nr;
        cubical_add(&mut p_nr, &in1, &n_r, &pairing_dlog_data.ixP);
        let in1 = n_rq;
        let in2 = n_r;
        cubical_dbl_add(
            &mut n_rq,
            &mut n_r,
            &in1,
            &in2,
            &pairing_dlog_data.ixQ,
            &pairing_dlog_data.A24,
        );

        let in1 = p_ns;
        cubical_add(&mut p_ns, &in1, &n_s, &pairing_dlog_data.ixP);
        let in1 = n_sq;
        let in2 = n_s;
        cubical_dbl_add(
            &mut n_sq,
            &mut n_s,
            &in1,
            &in2,
            &pairing_dlog_data.ixQ,
            &pairing_dlog_data.A24,
        );
    }

    let np_const = n_p;
    translate(&mut n_pq, &np_const);
    translate(&mut n_pr, &np_const);
    translate(&mut n_ps, &np_const);
    let nq_const = n_q;
    translate(&mut p_nq, &nq_const);
    translate(&mut r_nq, &nq_const);
    translate(&mut s_nq, &nq_const);
    let nr_const = n_r;
    translate(&mut p_nr, &nr_const);
    translate(&mut n_rq, &nr_const);
    let ns_const = n_s;
    translate(&mut p_ns, &ns_const);
    translate(&mut n_sq, &ns_const);

    let np_const = n_p;
    translate(&mut n_p, &np_const);
    let nq_const = n_q;
    translate(&mut n_q, &nq_const);
    let nr_const = n_r;
    translate(&mut n_r, &nr_const);
    let ns_const = n_s;
    translate(&mut n_s, &ns_const);

    // Reference Weil pairings.
    let mut t0 = EcPoint::zero();
    let mut t1 = EcPoint::zero();
    let mut w1: [Fp2; 5] = [Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    }; 5];
    let mut w2: [Fp2; 5] = w1;

    // e(P, Q) = w0; w1 and w2 are swapped relative to their later use
    // so that fp2_dlog_2e can consume w1 as the inverse.
    point_ratio(&mut t0, &n_pq, &n_p, &pairing_dlog_data.PQ.Q);
    point_ratio(&mut t1, &p_nq, &n_q, &pairing_dlog_data.PQ.P);
    fp2_mul(&mut w2[0], &t0.x, &t1.z);
    fp2_mul(&mut w1[0], &t1.x, &t0.z);

    // e(P, R) = w0^r2
    point_ratio(&mut t0, &n_pr, &n_p, &pairing_dlog_data.RS.P);
    point_ratio(&mut t1, &p_nr, &n_r, &pairing_dlog_data.PQ.P);
    fp2_mul(&mut w1[1], &t0.x, &t1.z);
    fp2_mul(&mut w2[1], &t1.x, &t0.z);

    // e(R, Q) = w0^r1
    point_ratio(&mut t0, &n_rq, &n_r, &pairing_dlog_data.PQ.Q);
    point_ratio(&mut t1, &r_nq, &n_q, &pairing_dlog_data.RS.P);
    fp2_mul(&mut w1[2], &t0.x, &t1.z);
    fp2_mul(&mut w2[2], &t1.x, &t0.z);

    // e(P, S) = w0^s2
    point_ratio(&mut t0, &n_ps, &n_p, &pairing_dlog_data.RS.Q);
    point_ratio(&mut t1, &p_ns, &n_s, &pairing_dlog_data.PQ.P);
    fp2_mul(&mut w1[3], &t0.x, &t1.z);
    fp2_mul(&mut w2[3], &t1.x, &t0.z);

    // e(S, Q) = w0^s1
    point_ratio(&mut t0, &n_sq, &n_s, &pairing_dlog_data.PQ.Q);
    point_ratio(&mut t1, &s_nq, &n_q, &pairing_dlog_data.RS.Q);
    fp2_mul(&mut w1[4], &t0.x, &t1.z);
    fp2_mul(&mut w2[4], &t1.x, &t0.z);

    fp2_batched_inv(&mut w1);
    for i in 0..5 {
        let wv = w1[i];
        fp2_mul(&mut w1[i], &wv, &w2[i]);
    }

    let w1_0 = w1[0];
    fp2_dlog_2e(r2, &w1[1], &w1_0, pairing_dlog_data.e as i32);
    fp2_dlog_2e(r1, &w1[2], &w1_0, pairing_dlog_data.e as i32);
    fp2_dlog_2e(s2, &w1[3], &w1_0, pairing_dlog_data.e as i32);
    fp2_dlog_2e(s1, &w1[4], &w1_0, pairing_dlog_data.e as i32);
}

/// Mirrors `ec_dlog_2_weil`. Given two bases `<P, Q>` and `<R, S>` of
/// E[2^e], compute scalars such that `R = r1*P + r2*Q` and
/// `S = s1*P + s2*Q`.
#[allow(clippy::too_many_arguments)]
pub fn ec_dlog_2_weil(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    pq: &mut EcBasis,
    rs: &EcBasis,
    curve: &mut EcCurve,
    e: i32,
) {
    ec_curve_normalize_a24(curve);

    let mut pairing_dlog_data = PairingDlogParams::zero();
    pairing_dlog_data.e = e as u32;
    pairing_dlog_data.PQ = *pq;
    pairing_dlog_data.RS = *rs;
    pairing_dlog_data.A24 = curve.A24;

    cubical_normalization_dlog(&mut pairing_dlog_data, curve);
    compute_difference_points(&mut pairing_dlog_data, curve);

    weil_dlog(r1, r2, s1, s2, &pairing_dlog_data);
}

/// Static helper mirroring `tate_dlog_partial` in biextension.c.
#[allow(clippy::needless_range_loop)]
fn tate_dlog_partial(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    pairing_dlog_data: &PairingDlogParams,
) {
    let e_full: u32 = TORSION_EVEN_POWER as u32;
    let e_diff: u32 = e_full - pairing_dlog_data.e;

    let mut n_p = EcPoint::zero();
    let mut n_q = EcPoint::zero();
    let mut n_r = EcPoint::zero();
    let mut n_s = EcPoint::zero();
    let mut n_pq = EcPoint::zero();
    let mut p_nr = EcPoint::zero();
    let mut p_ns = EcPoint::zero();
    let mut n_rq = EcPoint::zero();
    let mut n_sq = EcPoint::zero();

    copy_point(&mut n_p, &pairing_dlog_data.PQ.P);
    copy_point(&mut n_q, &pairing_dlog_data.PQ.Q);
    copy_point(&mut n_r, &pairing_dlog_data.RS.P);
    copy_point(&mut n_s, &pairing_dlog_data.RS.Q);
    copy_point(&mut n_pq, &pairing_dlog_data.PQ.PmQ);
    copy_point(&mut p_nr, &pairing_dlog_data.diff.PmR);
    copy_point(&mut p_ns, &pairing_dlog_data.diff.PmS);
    copy_point(&mut n_rq, &pairing_dlog_data.diff.RmQ);
    copy_point(&mut n_sq, &pairing_dlog_data.diff.SmQ);

    for _ in 0..(e_full - 1) {
        let in1 = n_pq;
        let in2 = n_p;
        cubical_dbl_add(
            &mut n_pq,
            &mut n_p,
            &in1,
            &in2,
            &pairing_dlog_data.ixQ,
            &pairing_dlog_data.A24,
        );
    }

    for _ in 0..(pairing_dlog_data.e - 1) {
        let in1 = p_nr;
        cubical_add(&mut p_nr, &in1, &n_r, &pairing_dlog_data.ixP);
        let in1 = n_rq;
        let in2 = n_r;
        cubical_dbl_add(
            &mut n_rq,
            &mut n_r,
            &in1,
            &in2,
            &pairing_dlog_data.ixQ,
            &pairing_dlog_data.A24,
        );

        let in1 = p_ns;
        cubical_add(&mut p_ns, &in1, &n_s, &pairing_dlog_data.ixP);
        let in1 = n_sq;
        let in2 = n_s;
        cubical_dbl_add(
            &mut n_sq,
            &mut n_s,
            &in1,
            &in2,
            &pairing_dlog_data.ixQ,
            &pairing_dlog_data.A24,
        );
    }

    let np_const = n_p;
    translate(&mut n_pq, &np_const);
    let nr_const = n_r;
    translate(&mut p_nr, &nr_const);
    translate(&mut n_rq, &nr_const);
    let ns_const = n_s;
    translate(&mut p_ns, &ns_const);
    translate(&mut n_sq, &ns_const);

    let np_const = n_p;
    translate(&mut n_p, &np_const);
    let nq_const = n_q;
    translate(&mut n_q, &nq_const);
    let nr_const = n_r;
    translate(&mut n_r, &nr_const);
    let ns_const = n_s;
    translate(&mut n_s, &ns_const);

    // Reference Tate pairings.
    let mut t0 = EcPoint::zero();
    let mut w1: [Fp2; 5] = [Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    }; 5];
    let mut w2: [Fp2; 5] = w1;

    point_ratio(&mut t0, &n_pq, &n_p, &pairing_dlog_data.PQ.Q);
    fp2_copy(&mut w1[0], &t0.x);
    fp2_copy(&mut w2[0], &t0.z);

    point_ratio(&mut t0, &p_nr, &n_r, &pairing_dlog_data.PQ.P);
    fp2_copy(&mut w1[1], &t0.x);
    fp2_copy(&mut w2[1], &t0.z);

    point_ratio(&mut t0, &n_rq, &n_r, &pairing_dlog_data.PQ.Q);
    fp2_copy(&mut w2[2], &t0.x);
    fp2_copy(&mut w1[2], &t0.z);

    point_ratio(&mut t0, &p_ns, &n_s, &pairing_dlog_data.PQ.P);
    fp2_copy(&mut w1[3], &t0.x);
    fp2_copy(&mut w2[3], &t0.z);

    point_ratio(&mut t0, &n_sq, &n_s, &pairing_dlog_data.PQ.Q);
    fp2_copy(&mut w2[4], &t0.x);
    fp2_copy(&mut w1[4], &t0.z);

    // Batched reduction using projective representation.
    for i in 0..5 {
        let mut frob = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        let mut tmp = frob;
        fp2_copy(&mut tmp, &w1[i]);
        fp2_frob(&mut frob, &w1[i]);
        fp2_mul(&mut w1[i], &w2[i], &frob);
        fp2_frob(&mut frob, &w2[i]);
        fp2_mul(&mut w2[i], &tmp, &frob);
    }

    fp2_batched_inv(&mut w2);
    for i in 0..5 {
        let wv = w1[i];
        fp2_mul(&mut w1[i], &wv, &w2[i]);
    }

    for i in 0..5 {
        let wv = w1[i];
        clear_cofac(&mut w1[i], &wv);
        for _ in 0..e_diff {
            let wv = w1[i];
            fp2_sqr(&mut w1[i], &wv);
        }
    }

    let w1_0 = w1[0];
    fp2_dlog_2e(r2, &w1[1], &w1_0, pairing_dlog_data.e as i32);
    fp2_dlog_2e(r1, &w1[2], &w1_0, pairing_dlog_data.e as i32);
    fp2_dlog_2e(s2, &w1[3], &w1_0, pairing_dlog_data.e as i32);
    fp2_dlog_2e(s1, &w1[4], &w1_0, pairing_dlog_data.e as i32);
}

/// Mirrors `ec_dlog_2_tate`. Given a full E[2^e_full] basis `<P, Q>`
/// and a smaller E[2^e] basis `<R, S>`, computes scalars such that
/// `R = r1*P + r2*Q` and `S = s1*P + s2*Q`.
#[allow(clippy::too_many_arguments)]
pub fn ec_dlog_2_tate(
    r1: &mut [u64],
    r2: &mut [u64],
    s1: &mut [u64],
    s2: &mut [u64],
    pq: &EcBasis,
    rs: &EcBasis,
    curve: &mut EcCurve,
    e: i32,
) {
    ec_curve_normalize_a24(curve);

    let mut pairing_dlog_data = PairingDlogParams::zero();
    pairing_dlog_data.e = e as u32;
    pairing_dlog_data.PQ = *pq;
    pairing_dlog_data.RS = *rs;
    pairing_dlog_data.A24 = curve.A24;

    cubical_normalization_dlog(&mut pairing_dlog_data, curve);
    compute_difference_points(&mut pairing_dlog_data, curve);
    tate_dlog_partial(r1, r2, s1, s2, &pairing_dlog_data);
}
