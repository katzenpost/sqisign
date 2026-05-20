//! Small-integer helpers (mirrors
//! `the-sqisign/src/quaternion/ref/generic/integers.c`).
//!
//! In scope this batch: `ibz_cornacchia_prime` only. The RNG-driven
//! `ibz_generate_random_prime` is deferred.

use crate::ibz::{
    ibz_add, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero, ibz_div, ibz_gcd, ibz_is_one,
    ibz_is_zero, ibz_mul, ibz_neg, ibz_set, ibz_sqrt, ibz_sqrt_mod_p, ibz_sub, Ibz,
};

/// `ibz_cornacchia_prime(x, y, n, p)`: solve `x^2 + n*y^2 == p` for
/// positive integers `x`, `y` assuming `p` is prime and `-n mod p` is a
/// square.
///
/// Returns 1 on success and writes `x, y`; returns 0 on failure. On
/// failure `x, y` are unspecified (the C reference leaves them in
/// scratch state too).
pub fn ibz_cornacchia_prime(x: &mut Ibz, y: &mut Ibz, n: &Ibz, p: &Ibz) -> i32 {
    // Scratch big-integers. r0/r1/prod are overwritten before being read in
    // the general path; the early-return special cases use only r2.
    let mut r0;
    let mut r1;
    let mut r2 = Ibz::zero();
    let mut a = Ibz::zero();
    let mut prod;

    // p == 2 special case.
    if ibz_cmp(p, &ibz_const_two()) == 0 {
        if ibz_is_one(n) == 1 {
            ibz_set(x, 1);
            ibz_set(y, 1);
            return 1;
        }
        return 0;
    }
    // p == n special case.
    if ibz_cmp(p, n) == 0 {
        ibz_set(x, 0);
        ibz_set(y, 1);
        return 1;
    }

    // Coprimality of p and n.
    ibz_gcd(&mut r2, p, n);
    if ibz_is_one(&r2) != 1 {
        return 0;
    }

    // Compute a square root of -n mod p.
    let mut neg_n = Ibz::zero();
    ibz_neg(&mut neg_n, n);
    if ibz_sqrt_mod_p(&mut r2, &neg_n, p) == 0 {
        return 0;
    }

    // The reference's main loop: while prod >= p, run the Euclidean step.
    // It begins by setting prod = r1 = r0 = p (so the loop body always
    // executes once, mirroring `do { ... } while (prod >= p)`).
    prod = p.clone();
    r1 = p.clone();
    r0 = p.clone();

    while ibz_cmp(&prod, p) >= 0 {
        // a, r0 := r2 / r1 (truncated)
        ibz_div(&mut a, &mut r0, &r2, &r1);
        // prod := r0 * r0
        ibz_mul(&mut prod, &r0, &r0);
        r2 = r1.clone();
        r1 = r0.clone();
    }

    // Test if the result is a solution.
    let mut a_tmp = Ibz::zero();
    ibz_sub(&mut a_tmp, p, &prod);
    let mut r2_tmp = Ibz::zero();
    ibz_div(&mut a, &mut r2_tmp, &a_tmp, n);
    if ibz_is_zero(&r2_tmp) != 1 {
        return 0;
    }
    if ibz_sqrt(y, &a) == 0 {
        return 0;
    }

    // x := r0
    *x = r0.clone();

    // Verify: prod + y*y*n == p.
    let mut yy = Ibz::zero();
    ibz_mul(&mut yy, y, y);
    let mut yyn = Ibz::zero();
    ibz_mul(&mut yyn, &yy, n);
    let new_prod = {
        let mut t = Ibz::zero();
        ibz_add(&mut t, &prod, &yyn);
        t
    };
    if ibz_cmp(&new_prod, p) == 0 {
        // Touch the unused-imports to keep the symbol table mirror tidy
        // (mirrors C's `ibz_const_zero`/`ibz_const_one` exports which are
        // referenced by name even when arithmetic uses literals).
        let _ = ibz_const_zero();
        let _ = ibz_const_one();
        1
    } else {
        0
    }
}
