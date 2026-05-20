//! Quaternion algebra and its elements.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion/ref/generic/algebra.c` and
//! the constructors from `finit.c`. The quaternion algebra modelled here
//! is ramified at `p = 3 (mod 4)` and infinity; elements are represented
//! in the basis `(1, i, j, ij)` with `i^2 = -1` and `j^2 = -p`.
//!
//! See the crate root for the list of in-scope and deferred functions.

use num_bigint::BigInt;
use num_traits::{One, Zero};

use crate::dim4::{
    ibz_vec_4_add, ibz_vec_4_content, ibz_vec_4_is_zero, ibz_vec_4_scalar_div,
    ibz_vec_4_scalar_mul, ibz_vec_4_sub,
};
use crate::ibz::{
    ibz_abs, ibz_add, ibz_cmp, ibz_const_zero, ibz_copy_digits as _ibz_copy_digits, ibz_div,
    ibz_gcd, ibz_mul, ibz_neg, ibz_set, ibz_sub, Ibz,
};
use crate::lattice::{quat_lattice_contains, QuatLattice};

/// `quat_alg_t`: the quaternion algebra ramified at `p` (and infinity).
///
/// `p` must be prime and congruent to `3 (mod 4)`. We do not enforce this
/// at construction; the constraint is documented per the C reference.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct QuatAlg {
    /// The ramified prime.
    pub p: Ibz,
}

impl QuatAlg {
    /// `quat_alg_init_set(alg, p)`: construct over the prime `p`.
    pub fn init_set(p: &Ibz) -> Self {
        Self { p: p.clone() }
    }

    /// `quat_alg_init_set_ui(alg, p)`: construct over the prime `p`
    /// passed as a small unsigned integer.
    pub fn init_set_ui(p: u32) -> Self {
        Self {
            p: Ibz(BigInt::from(p)),
        }
    }
}

/// `quat_alg_elem_t`: an element of the quaternion algebra represented as
/// four numerators in the basis `(1, i, j, ij)` over a common
/// denominator.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct QuatAlgElem {
    /// Common denominator. Must be non-zero.
    pub denom: Ibz,
    /// The four numerator coordinates.
    pub coord: [Ibz; 4],
}

impl Default for QuatAlgElem {
    fn default() -> Self {
        Self::new()
    }
}

impl QuatAlgElem {
    /// `quat_alg_elem_init(elem)`: zero coordinates, denominator one.
    pub fn new() -> Self {
        Self {
            denom: Ibz(BigInt::one()),
            coord: [
                Ibz(BigInt::zero()),
                Ibz(BigInt::zero()),
                Ibz(BigInt::zero()),
                Ibz(BigInt::zero()),
            ],
        }
    }
}

/// `quat_alg_coord_mul(res, a, b, alg)`: multiply two coordinate vectors
/// as elements of the algebra in the standard basis.
///
/// The product follows the multiplication table `i^2 = -1`, `j^2 = -p`,
/// `ij = -ji`. The implementation tracks the C reference line-for-line.
pub fn quat_alg_coord_mul(res: &mut [Ibz; 4], a: &[Ibz; 4], b: &[Ibz; 4], alg: &QuatAlg) {
    let mut prod = Ibz::zero();
    let mut sum = [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()];

    // Coordinate 0 (the "1" component).
    ibz_mul(&mut prod, &a[2], &b[2]);
    ibz_sub_clone_into(&mut sum[0], &prod, true);
    ibz_mul(&mut prod, &a[3], &b[3]);
    ibz_sub_clone_into(&mut sum[0], &prod, true);
    // Multiply by p (the i = -1, j^2 = -p convention).
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &sum[0], &alg.p);
    sum[0] = tmp;
    ibz_mul(&mut prod, &a[0], &b[0]);
    ibz_add_clone_into(&mut sum[0], &prod);
    ibz_mul(&mut prod, &a[1], &b[1]);
    ibz_sub_clone_into(&mut sum[0], &prod, true);

    // Coordinate 1 (the "i" component).
    ibz_mul(&mut prod, &a[2], &b[3]);
    ibz_add_clone_into(&mut sum[1], &prod);
    ibz_mul(&mut prod, &a[3], &b[2]);
    ibz_sub_clone_into(&mut sum[1], &prod, true);
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &sum[1], &alg.p);
    sum[1] = tmp;
    ibz_mul(&mut prod, &a[0], &b[1]);
    ibz_add_clone_into(&mut sum[1], &prod);
    ibz_mul(&mut prod, &a[1], &b[0]);
    ibz_add_clone_into(&mut sum[1], &prod);

    // Coordinate 2 (the "j" component).
    ibz_mul(&mut prod, &a[0], &b[2]);
    ibz_add_clone_into(&mut sum[2], &prod);
    ibz_mul(&mut prod, &a[2], &b[0]);
    ibz_add_clone_into(&mut sum[2], &prod);
    ibz_mul(&mut prod, &a[1], &b[3]);
    ibz_sub_clone_into(&mut sum[2], &prod, true);
    ibz_mul(&mut prod, &a[3], &b[1]);
    ibz_add_clone_into(&mut sum[2], &prod);

    // Coordinate 3 (the "ij" component).
    ibz_mul(&mut prod, &a[0], &b[3]);
    ibz_add_clone_into(&mut sum[3], &prod);
    ibz_mul(&mut prod, &a[3], &b[0]);
    ibz_add_clone_into(&mut sum[3], &prod);
    ibz_mul(&mut prod, &a[2], &b[1]);
    ibz_sub_clone_into(&mut sum[3], &prod, true);
    ibz_mul(&mut prod, &a[1], &b[2]);
    ibz_add_clone_into(&mut sum[3], &prod);

    res[0] = sum[0].clone();
    res[1] = sum[1].clone();
    res[2] = sum[2].clone();
    res[3] = sum[3].clone();
}

// Helpers to express the C pattern `sum -= prod` / `sum += prod` without
// repeatedly allocating temporaries through ibz_add/ibz_sub (which require
// distinct &mut and & borrows).
fn ibz_add_clone_into(acc: &mut Ibz, addend: &Ibz) {
    let mut tmp = Ibz::zero();
    ibz_add(&mut tmp, acc, addend);
    *acc = tmp;
}
fn ibz_sub_clone_into(acc: &mut Ibz, subtrahend: &Ibz, negate: bool) {
    let mut tmp = Ibz::zero();
    if negate {
        ibz_sub(&mut tmp, acc, subtrahend);
    } else {
        ibz_add(&mut tmp, acc, subtrahend);
    }
    *acc = tmp;
}

/// `quat_alg_equal_denom(res_a, res_b, a, b)`: rewrite `a` and `b` to a
/// common denominator. The numerators and denominators of the outputs
/// are set so that the represented quaternion values are unchanged.
pub fn quat_alg_equal_denom(
    res_a: &mut QuatAlgElem,
    res_b: &mut QuatAlgElem,
    a: &QuatAlgElem,
    b: &QuatAlgElem,
) {
    let mut gcd = Ibz::zero();
    let mut r = Ibz::zero();
    ibz_gcd(&mut gcd, &a.denom, &b.denom);

    // Temporarily store the reduced denominators in res_a.denom and
    // res_b.denom; the C reference uses the same trick to avoid an
    // extra pair of temporaries.
    ibz_div(&mut res_a.denom, &mut r, &a.denom, &gcd);
    ibz_div(&mut res_b.denom, &mut r, &b.denom, &gcd);

    for i in 0..4 {
        let red_a = res_a.denom.clone();
        let red_b = res_b.denom.clone();
        ibz_mul(&mut res_a.coord[i], &a.coord[i], &red_b);
        ibz_mul(&mut res_b.coord[i], &b.coord[i], &red_a);
    }
    // Now build the common denominator.
    let red_a = res_a.denom.clone();
    let red_b = res_b.denom.clone();
    let mut prod = Ibz::zero();
    ibz_mul(&mut prod, &red_a, &red_b);
    res_b.denom = prod.clone();
    res_a.denom = prod;
    // Multiply by the gcd.
    let g = gcd.clone();
    let mut t = Ibz::zero();
    ibz_mul(&mut t, &res_b.denom, &g);
    res_b.denom = t;
    let mut t = Ibz::zero();
    ibz_mul(&mut t, &res_a.denom, &g);
    res_a.denom = t;
}

/// `quat_alg_mul(res, a, b, alg)`: multiply two algebra elements.
pub fn quat_alg_mul(res: &mut QuatAlgElem, a: &QuatAlgElem, b: &QuatAlgElem, alg: &QuatAlg) {
    let mut denom = Ibz::zero();
    ibz_mul(&mut denom, &a.denom, &b.denom);
    res.denom = denom;
    quat_alg_coord_mul(&mut res.coord, &a.coord, &b.coord, alg);
}

/// `quat_alg_norm(res_num, res_denom, a, alg)`: reduced norm of `a`,
/// returned as a reduced fraction `(res_num / res_denom)` with
/// `res_denom > 0`.
pub fn quat_alg_norm(res_num: &mut Ibz, res_denom: &mut Ibz, a: &QuatAlgElem, alg: &QuatAlg) {
    let mut norm = QuatAlgElem::new();
    quat_alg_conj(&mut norm, a);
    let conj = norm.clone();
    quat_alg_mul(&mut norm, a, &conj, alg);

    let mut g = Ibz::zero();
    ibz_gcd(&mut g, &norm.coord[0], &norm.denom);

    let mut r = Ibz::zero();
    ibz_div(res_num, &mut r, &norm.coord[0], &g);
    ibz_div(res_denom, &mut r, &norm.denom, &g);

    let mut tmp = Ibz::zero();
    ibz_abs(&mut tmp, res_denom);
    *res_denom = tmp;
    let mut tmp = Ibz::zero();
    ibz_abs(&mut tmp, res_num);
    *res_num = tmp;

    debug_assert!(
        crate::ibz::ibz_cmp(res_denom, &crate::ibz::ibz_const_zero()) > 0,
        "quat_alg_norm: denominator must be strictly positive"
    );
}

/// `quat_alg_scalar(elem, numerator, denominator)`: build the scalar
/// algebra element `numerator / denominator` (purely in the `1`
/// component).
pub fn quat_alg_scalar(elem: &mut QuatAlgElem, numerator: &Ibz, denominator: &Ibz) {
    elem.denom = denominator.clone();
    elem.coord[0] = numerator.clone();
    ibz_set(&mut elem.coord[1], 0);
    ibz_set(&mut elem.coord[2], 0);
    ibz_set(&mut elem.coord[3], 0);
}

/// `quat_alg_conj(conj, x)`: the conjugate, negating the `i`, `j`, `ij`
/// coordinates.
pub fn quat_alg_conj(conj: &mut QuatAlgElem, x: &QuatAlgElem) {
    conj.denom = x.denom.clone();
    conj.coord[0] = x.coord[0].clone();
    let mut tmp = Ibz::zero();
    ibz_neg(&mut tmp, &x.coord[1]);
    conj.coord[1] = tmp;
    let mut tmp = Ibz::zero();
    ibz_neg(&mut tmp, &x.coord[2]);
    conj.coord[2] = tmp;
    let mut tmp = Ibz::zero();
    ibz_neg(&mut tmp, &x.coord[3]);
    conj.coord[3] = tmp;
}

/// `quat_alg_elem_set(elem, denom, c0, c1, c2, c3)`: set the four
/// coordinates and denominator from small signed integers, without
/// normalization.
pub fn quat_alg_elem_set(elem: &mut QuatAlgElem, denom: i32, c0: i32, c1: i32, c2: i32, c3: i32) {
    ibz_set(&mut elem.coord[0], c0);
    ibz_set(&mut elem.coord[1], c1);
    ibz_set(&mut elem.coord[2], c2);
    ibz_set(&mut elem.coord[3], c3);
    ibz_set(&mut elem.denom, denom);
}

/// `quat_alg_elem_copy(copy, copied)`: deep-copy an element.
pub fn quat_alg_elem_copy(copy: &mut QuatAlgElem, copied: &QuatAlgElem) {
    copy.denom = copied.denom.clone();
    copy.coord[0] = copied.coord[0].clone();
    copy.coord[1] = copied.coord[1].clone();
    copy.coord[2] = copied.coord[2].clone();
    copy.coord[3] = copied.coord[3].clone();
}

/// `quat_alg_elem_copy_ibz(elem, denom, c0, c1, c2, c3)`: set the
/// coordinates and denominator from big-integer inputs without
/// normalization.
pub fn quat_alg_elem_copy_ibz(
    elem: &mut QuatAlgElem,
    denom: &Ibz,
    c0: &Ibz,
    c1: &Ibz,
    c2: &Ibz,
    c3: &Ibz,
) {
    elem.coord[0] = c0.clone();
    elem.coord[1] = c1.clone();
    elem.coord[2] = c2.clone();
    elem.coord[3] = c3.clone();
    elem.denom = denom.clone();
}

/// `quat_alg_elem_mul_by_scalar(res, scalar, elem)`: scale every
/// coordinate by `scalar` (the denominator is preserved).
pub fn quat_alg_elem_mul_by_scalar(res: &mut QuatAlgElem, scalar: &Ibz, elem: &QuatAlgElem) {
    for i in 0..4 {
        ibz_mul(&mut res.coord[i], &elem.coord[i], scalar);
    }
    res.denom = elem.denom.clone();
}

/// `quat_alg_add(res, a, b)`: sum of two algebra elements, written into
/// `res` with a common denominator.
pub fn quat_alg_add(res: &mut QuatAlgElem, a: &QuatAlgElem, b: &QuatAlgElem) {
    let mut res_a = QuatAlgElem::new();
    let mut res_b = QuatAlgElem::new();
    quat_alg_equal_denom(&mut res_a, &mut res_b, a, b);
    res.denom = res_a.denom.clone();
    ibz_vec_4_add(&mut res.coord, &res_a.coord, &res_b.coord);
}

/// `quat_alg_sub(res, a, b)`.
pub fn quat_alg_sub(res: &mut QuatAlgElem, a: &QuatAlgElem, b: &QuatAlgElem) {
    let mut res_a = QuatAlgElem::new();
    let mut res_b = QuatAlgElem::new();
    quat_alg_equal_denom(&mut res_a, &mut res_b, a, b);
    res.denom = res_a.denom.clone();
    ibz_vec_4_sub(&mut res.coord, &res_a.coord, &res_b.coord);
}

/// `quat_alg_normalize(x)`: divide content of coord and denom by their
/// joint gcd, then sign-flip so the denominator is positive.
pub fn quat_alg_normalize(x: &mut QuatAlgElem) {
    let mut gcd = Ibz::zero();
    ibz_vec_4_content(&mut gcd, &x.coord);
    let mut t = Ibz::zero();
    ibz_gcd(&mut t, &gcd, &x.denom);
    gcd = t;
    let mut q = Ibz::zero();
    let mut r = Ibz::zero();
    ibz_div(&mut q, &mut r, &x.denom, &gcd);
    x.denom = q;
    let cloned_coord = clone_coord(&x.coord);
    let _ = ibz_vec_4_scalar_div(&mut x.coord, &gcd, &cloned_coord);
    // sign = 2 * (0 > cmp(0, denom)) - 1, i.e. -1 if denom > 0 else 1.
    // Wait, that's backwards. Let me re-check: the C is
    //   ibz_set(&sign, 2 * (0 > ibz_cmp(&ibz_const_zero, &(x->denom))) - 1);
    // ibz_cmp(zero, denom) > 0 iff zero > denom iff denom < 0. The outer
    // `0 > ...` makes this: ibz_cmp(zero, denom) < 0, i.e. zero < denom
    // iff denom > 0. So sign = 1 if denom > 0, else -1.
    //
    // Hmm wait again. `0 > cmp(zero, denom)`: that's "cmp result negative",
    // which means zero < denom, i.e. denom > 0. So sign = 2*1 - 1 = 1 when
    // denom > 0, else 2*0 - 1 = -1.
    //
    // We want: if denom < 0, negate both. So sign should be -1 when denom < 0,
    // matching the formula. Correct.
    let cmp = ibz_cmp(&ibz_const_zero(), &x.denom);
    let pred = if 0 > cmp { 1 } else { 0 };
    let mut sign = Ibz::zero();
    ibz_set(&mut sign, 2 * pred - 1);
    let cloned_coord = clone_coord(&x.coord);
    ibz_vec_4_scalar_mul(&mut x.coord, &sign, &cloned_coord);
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &sign, &x.denom);
    x.denom = tmp;
}

/// `quat_alg_elem_is_zero(x)`.
pub fn quat_alg_elem_is_zero(x: &QuatAlgElem) -> i32 {
    ibz_vec_4_is_zero(&x.coord)
}

/// `quat_alg_elem_equal(a, b)`.
pub fn quat_alg_elem_equal(a: &QuatAlgElem, b: &QuatAlgElem) -> i32 {
    let mut diff = QuatAlgElem::new();
    quat_alg_sub(&mut diff, a, b);
    quat_alg_elem_is_zero(&diff)
}

/// `quat_alg_make_primitive(primitive_x, content, x, order)`.
pub fn quat_alg_make_primitive(
    primitive_x: &mut crate::dim4::IbzVec4,
    content: &mut Ibz,
    x: &QuatAlgElem,
    order: &QuatLattice,
) {
    let ok = quat_lattice_contains(Some(primitive_x), order, x);
    assert!(ok != 0, "quat_alg_make_primitive: x must be in order");
    ibz_vec_4_content(content, primitive_x);
    let mut r = Ibz::zero();
    for i in 0..4 {
        let saved = primitive_x[i].clone();
        ibz_div(&mut primitive_x[i], &mut r, &saved, content);
    }
}

fn clone_coord(c: &[Ibz; 4]) -> [Ibz; 4] {
    [c[0].clone(), c[1].clone(), c[2].clone(), c[3].clone()]
}

// Suppress unused-import lints on items used only by intra-crate callers.
#[doc(hidden)]
#[allow(dead_code)]
fn _silence_unused() {
    let _ = _ibz_copy_digits;
    let _ = ibz_abs;
}
