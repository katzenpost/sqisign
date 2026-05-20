//! SQIsign `ec`: elliptic curve point operations on the Kummer line and
//! Jacobian coordinates, plus degree-2 and degree-4 isogeny construction
//! and evaluation.
//!
//! Mirrors `vendor/the-sqisign/src/ec`. Phase 1 unit 4 (the present
//! delegation) ports the core curve and point arithmetic from
//! `lvlx/ec.c`, the Jacobian arithmetic from `lvlx/ec_jac.c`, and the
//! short isogeny primitives from `lvlx/xeval.c` and `lvlx/xisog.c`. The
//! biextension pairing, basis lifting and recovery, and the long
//! isogeny chains in `basis.c`, `biextension.c`, and `isog_chains.c`
//! are deferred to a follow-up batch.
//!
//! The differential boundary is the raw fp2 limb representation: every
//! ported function replays the recorded C output bit-for-bit on the
//! committed vector battery, identical in spirit to the gf port.
//!
//! ## Defects observed
//!
//! - `xeval_4_singular` and `xisog_4_singular` are declared in
//!   `vendor/the-sqisign/src/ec/ref/include/isog.h` but never defined
//!   anywhere in the upstream tree (no callers either): they are dead
//!   declarations. The port skips them; if upstream ever adds bodies,
//!   we add the boundaries then.

#![forbid(unsafe_code)]
#![allow(non_snake_case)]

use sqisign_gf::{
    fp2_add, fp2_add_one, fp2_copy, fp2_cswap, fp2_inv, fp2_is_equal, fp2_is_one, fp2_is_zero,
    fp2_mul, fp2_neg, fp2_select, fp2_set_one, fp2_set_small, fp2_set_zero, fp2_sqr, fp2_sqrt,
    fp2_sub, fp_add, fp_copy, fp_div3, fp_neg, fp_set_one, fp_sub, Fp, Fp2, NWORDS_FIELD,
};
use sqisign_mp::{mp_shiftr, mp_sub, select_ct, swap_ct};

// ---------------------------------------------------------------------------
// ec_params.c (lvl1)
// ---------------------------------------------------------------------------

/// The power of two dividing `p + 1`. Mirrors the C macro
/// `TORSION_EVEN_POWER` defined in
/// `vendor/the-sqisign/src/precomp/ref/lvl1/include/ec_params.h`.
pub const TORSION_EVEN_POWER: usize = 248;

/// The odd cofactor `(p + 1) / 2^TORSION_EVEN_POWER`, kept as a
/// single-element array to mirror the C declaration
/// `const digit_t p_cofactor_for_2f[1] = {5};` in
/// `vendor/the-sqisign/src/precomp/ref/lvl1/ec_params.c`.
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
