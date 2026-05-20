//! Rational arithmetic over `Ibz` (the `ibq_t` type).
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/lll/rationals.c`.
//! Used solely by the LLL verification path (the L2 core itself works in
//! `dpe_t`); the rationals are exact, with no reductions to lowest terms
//! except via the explicit [`ibq_reduce`].
//!
//! The C type `ibq_t = ibz_t[2]` becomes a Rust struct of two `Ibz` values.
//! Convention: `numerator` is index 0, `denominator` is index 1, matching
//! the C reference's `&((*x)[0])` / `&((*x)[1])` access pattern.

use crate::ibz::{
    ibz_add as ibz_add_fn, ibz_cmp, ibz_const_zero, ibz_div, ibz_gcd, ibz_is_zero, ibz_mod,
    ibz_mul, ibz_neg, Ibz,
};

/// `ibq_t`: a numerator/denominator pair over `Ibz`. The denominator is
/// initialised to 1 by [`Ibq::new`].
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Ibq {
    pub num: Ibz,
    pub den: Ibz,
}

impl Ibq {
    /// `ibq_init`. Default denominator is `1`.
    pub fn new() -> Self {
        Self {
            num: Ibz::zero(),
            den: Ibz::from(1),
        }
    }
}

/// `ibq_mat_4x4_t`.
pub type IbqMat4x4 = [[Ibq; 4]; 4];

/// `ibq_vec_4_t`.
pub type IbqVec4 = [Ibq; 4];

/// Allocate a zero 4-vector of rationals, each initialised to `0/1`.
pub fn ibq_vec_4_new() -> IbqVec4 {
    [Ibq::new(), Ibq::new(), Ibq::new(), Ibq::new()]
}

/// Allocate a 4x4 matrix of rationals.
pub fn ibq_mat_4x4_new() -> IbqMat4x4 {
    [
        ibq_vec_4_new(),
        ibq_vec_4_new(),
        ibq_vec_4_new(),
        ibq_vec_4_new(),
    ]
}

/// `ibq_reduce(x)`: divide numerator and denominator by their gcd.
pub fn ibq_reduce(x: &mut Ibq) {
    let mut g = Ibz::zero();
    ibz_gcd(&mut g, &x.num, &x.den);
    let mut qn = Ibz::zero();
    let mut qd = Ibz::zero();
    let mut r = Ibz::zero();
    ibz_div(&mut qn, &mut r, &x.num, &g);
    debug_assert!(ibz_is_zero(&r) == 1);
    ibz_div(&mut qd, &mut r, &x.den, &g);
    debug_assert!(ibz_is_zero(&r) == 1);
    x.num = qn;
    x.den = qd;
}

/// `ibq_add(sum, a, b)`: `sum = a/b1 + c/d1 = (a*d1 + c*b1)/(b1*d1)`.
pub fn ibq_add(sum: &mut Ibq, a: &Ibq, b: &Ibq) {
    let mut tmp1 = Ibz::zero();
    let mut tmp2 = Ibz::zero();
    ibz_mul(&mut tmp1, &a.num, &b.den);
    ibz_mul(&mut tmp2, &b.num, &a.den);
    let mut n = Ibz::zero();
    ibz_add_fn(&mut n, &tmp1, &tmp2);
    let mut d = Ibz::zero();
    ibz_mul(&mut d, &a.den, &b.den);
    sum.num = n;
    sum.den = d;
}

/// `ibq_neg(neg, x)`: `neg = -x` (denominator kept, numerator negated).
pub fn ibq_neg(neg: &mut Ibq, x: &Ibq) {
    neg.den = x.den.clone();
    let mut n = Ibz::zero();
    ibz_neg(&mut n, &x.num);
    neg.num = n;
}

/// `ibq_sub(diff, a, b)`: `diff = a - b`.
pub fn ibq_sub(diff: &mut Ibq, a: &Ibq, b: &Ibq) {
    let mut neg = Ibq::new();
    ibq_neg(&mut neg, b);
    ibq_add(diff, a, &neg);
}

/// `ibq_abs(abs, x)`: `abs = |x|`. The C reference uses `ibq_cmp` against
/// the negation to pick the larger.
pub fn ibq_abs(abs: &mut Ibq, x: &Ibq) {
    let mut neg = Ibq::new();
    ibq_neg(&mut neg, x);
    if ibq_cmp(x, &neg) < 0 {
        ibq_copy(abs, &neg);
    } else {
        ibq_copy(abs, x);
    }
}

/// `ibq_mul(prod, a, b)`: `prod = a*b` (no gcd reduction).
pub fn ibq_mul(prod: &mut Ibq, a: &Ibq, b: &Ibq) {
    let mut n = Ibz::zero();
    let mut d = Ibz::zero();
    ibz_mul(&mut n, &a.num, &b.num);
    ibz_mul(&mut d, &a.den, &b.den);
    prod.num = n;
    prod.den = d;
}

/// `ibq_inv(inv, x)`: swap numerator and denominator. Returns 0 iff `x` is
/// zero (so the inverse does not exist), 1 otherwise.
pub fn ibq_inv(inv: &mut Ibq, x: &Ibq) -> i32 {
    let res = (ibq_is_zero(x) == 0) as i32;
    if res != 0 {
        inv.num = x.den.clone();
        inv.den = x.num.clone();
    }
    res
}

/// `ibq_cmp(a, b)`: return a positive value if `a > b`, zero if equal, and
/// a negative value if `a < b`. Matches the C convention exactly, including
/// the two-step sign-flip in the upstream code.
pub fn ibq_cmp(a: &Ibq, b: &Ibq) -> i32 {
    let mut x = a.num.clone();
    let mut y = b.num.clone();
    // y *= a.den; x *= b.den
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &y, &a.den);
    y = tmp;
    let mut tmp = Ibz::zero();
    ibz_mul(&mut tmp, &x, &b.den);
    x = tmp;
    // The C reference flips signs once if a.den > 0 and once again if
    // b.den > 0. That double-flip is a no-op for positive denominators
    // (the normal case); if one is negative the cross-multiply flips the
    // inequality once, which matches the rational ordering.
    let zero = ibz_const_zero();
    if ibz_cmp(&a.den, &zero) > 0 {
        let mut t = Ibz::zero();
        ibz_neg(&mut t, &y);
        y = t;
        let mut t = Ibz::zero();
        ibz_neg(&mut t, &x);
        x = t;
    }
    if ibz_cmp(&b.den, &zero) > 0 {
        let mut t = Ibz::zero();
        ibz_neg(&mut t, &y);
        y = t;
        let mut t = Ibz::zero();
        ibz_neg(&mut t, &x);
        x = t;
    }
    ibz_cmp(&x, &y)
}

/// `ibq_is_zero(x)`.
pub fn ibq_is_zero(x: &Ibq) -> i32 {
    ibz_is_zero(&x.num)
}

/// `ibq_is_one(x)`.
pub fn ibq_is_one(x: &Ibq) -> i32 {
    (ibz_cmp(&x.num, &x.den) == 0) as i32
}

/// `ibq_set(q, a, b)`: `q = a / b`. Returns 0 iff `b == 0`.
pub fn ibq_set(q: &mut Ibq, a: &Ibz, b: &Ibz) -> i32 {
    q.num = a.clone();
    q.den = b.clone();
    (ibz_is_zero(b) == 0) as i32
}

/// `ibq_copy(target, value)`.
pub fn ibq_copy(target: &mut Ibq, value: &Ibq) {
    target.num = value.num.clone();
    target.den = value.den.clone();
}

/// `ibq_is_ibz(q)`: is `q.num % q.den == 0`?
pub fn ibq_is_ibz(q: &Ibq) -> i32 {
    let mut r = Ibz::zero();
    ibz_mod(&mut r, &q.num, &q.den);
    ibz_is_zero(&r)
}

/// `ibq_to_ibz(z, q)`: divide if integral. Returns 1 if exact, 0 otherwise.
pub fn ibq_to_ibz(z: &mut Ibz, q: &Ibq) -> i32 {
    let mut r = Ibz::zero();
    ibz_div(z, &mut r, &q.num, &q.den);
    ibz_is_zero(&r)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_pair(a: i32, b: i32) -> Ibq {
        let mut q = Ibq::new();
        ibq_set(&mut q, &Ibz::from(a), &Ibz::from(b));
        q
    }

    #[test]
    fn cmp_positive_denoms() {
        // 1/2 vs 1/3 -> 1/2 > 1/3
        let a = from_pair(1, 2);
        let b = from_pair(1, 3);
        assert!(ibq_cmp(&a, &b) > 0);
        assert!(ibq_cmp(&b, &a) < 0);
        assert_eq!(ibq_cmp(&a, &a), 0);
    }

    #[test]
    fn add_then_reduce() {
        let a = from_pair(1, 2);
        let b = from_pair(1, 3);
        let mut s = Ibq::new();
        ibq_add(&mut s, &a, &b);
        // 1/2 + 1/3 = 5/6
        ibq_reduce(&mut s);
        assert_eq!(s.num.0, 5.into());
        assert_eq!(s.den.0, 6.into());
    }

    #[test]
    fn mul_then_reduce() {
        let a = from_pair(2, 3);
        let b = from_pair(3, 4);
        let mut p = Ibq::new();
        ibq_mul(&mut p, &a, &b);
        ibq_reduce(&mut p);
        // 2/3 * 3/4 = 1/2
        assert_eq!(p.num.0, 1.into());
        assert_eq!(p.den.0, 2.into());
    }

    #[test]
    fn inv_nonzero() {
        let a = from_pair(5, 7);
        let mut inv = Ibq::new();
        assert_eq!(ibq_inv(&mut inv, &a), 1);
        assert_eq!(inv.num.0, 7.into());
        assert_eq!(inv.den.0, 5.into());
    }

    #[test]
    fn inv_zero() {
        let a = from_pair(0, 1);
        let mut inv = Ibq::new();
        assert_eq!(ibq_inv(&mut inv, &a), 0);
    }
}
