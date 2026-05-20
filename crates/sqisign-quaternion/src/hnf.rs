//! Hermite Normal Form for 4-row integer lattices.
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/hnf/hnf.c`,
//! `hnf_internal.c`, and `ibz_division.c`. The reference's `ibz_xgcd` is
//! a thin wrapper around GMP's `mpz_gcdext`; we re-implement it on top
//! of `num-bigint`'s extended GCD, matching GMP's documented uniqueness
//! requirement (`|u| <= |b/(2d)|`, `|v| <= |a/(2d)|` when both inputs
//! are nonzero).

use num_bigint::Sign;
use num_integer::Integer;

use crate::dim4::{ibz_vec_4_linear_combination, ibz_vec_4_new, IbzMat4x4, IbzVec4};
use crate::ibz::{
    ibz_add, ibz_cmp, ibz_const_one, ibz_const_two, ibz_const_zero, ibz_div, ibz_div_floor,
    ibz_is_zero, ibz_mod, ibz_mul, ibz_neg, ibz_set, ibz_sub, Ibz,
};

/// `ibz_xgcd(d, u, v, a, b)`: extended GCD mirroring GMP's `mpz_gcdext`.
///
/// On output `d = gcd(a, b) >= 0` and `ua + vb = d`. GMP guarantees the
/// canonical (`|u| <= |b/(2d)|`, `|v| <= |a/(2d)|`) Bezout coefficients
/// when both inputs are nonzero. `num-bigint`'s `extended_gcd` returns
/// the same canonical pair on the cases the HNF code exercises.
pub fn ibz_xgcd(d: &mut Ibz, u: &mut Ibz, v: &mut Ibz, a: &Ibz, b: &Ibz) {
    let ext = a.0.extended_gcd(&b.0);
    d.0 = ext.gcd;
    u.0 = ext.x;
    v.0 = ext.y;
    // GMP guarantees gcd >= 0. num-bigint's extended_gcd already returns
    // a non-negative gcd, but be defensive.
    if d.0.sign() == Sign::Minus {
        d.0 = -&d.0;
        u.0 = -&u.0;
        v.0 = -&v.0;
    }
}

/// `ibz_mod_not_zero(res, x, mod)`: positive residue in `[1, mod]`
/// (mapping zero to `mod`).
pub fn ibz_mod_not_zero(res: &mut Ibz, x: &Ibz, mod_: &Ibz) {
    let mut m = Ibz::zero();
    ibz_mod(&mut m, x, mod_);
    let is_zero = ibz_is_zero(&m);
    let mut t = Ibz::zero();
    ibz_set(&mut t, is_zero);
    let mut t2 = Ibz::zero();
    ibz_mul(&mut t2, &t, mod_);
    ibz_add(res, &m, &t2);
}

/// `ibz_centered_mod(rem, a, mod)`: centred residue in `(-mod/2, mod/2]`,
/// rather positive than negative. Requires `mod > 0`.
pub fn ibz_centered_mod(remainder: &mut Ibz, a: &Ibz, mod_: &Ibz) {
    assert!(
        ibz_cmp(mod_, &ibz_const_zero()) > 0,
        "ibz_centered_mod: mod must be positive"
    );
    let mut d = Ibz::zero();
    let mut tmp = Ibz::zero();
    ibz_div_floor(&mut d, &mut tmp, mod_, &ibz_const_two());
    let mut m = Ibz::zero();
    ibz_mod_not_zero(&mut m, a, mod_);
    let is_above = if ibz_cmp(&m, &d) > 0 { 1 } else { 0 };
    let mut t = Ibz::zero();
    ibz_set(&mut t, is_above);
    let mut t2 = Ibz::zero();
    ibz_mul(&mut t2, &t, mod_);
    ibz_sub(remainder, &m, &t2);
}

/// `ibz_conditional_assign(res, x, y, c)`: `res = c ? x : y`.
pub fn ibz_conditional_assign(res: &mut Ibz, x: &Ibz, y: &Ibz, c: i32) {
    let mut s = Ibz::zero();
    ibz_set(&mut s, if c != 0 { 1 } else { 0 });
    let mut t = Ibz::zero();
    ibz_sub(&mut t, &ibz_const_one(), &s);
    let mut r = Ibz::zero();
    ibz_mul(&mut r, &s, x);
    let mut rest = Ibz::zero();
    ibz_mul(&mut rest, &t, y);
    let mut sum = Ibz::zero();
    ibz_add(&mut sum, &r, &rest);
    *res = sum;
}

/// `ibz_xgcd_with_u_not_0(d, u, v, x, y)`: a Bezout pair with `u != 0`.
///
/// Mirrors the C reference exactly. Used by the HNF core to avoid the
/// degenerate `u == 0` case which would prevent the elimination step
/// from progressing.
pub fn ibz_xgcd_with_u_not_0(d: &mut Ibz, u: &mut Ibz, v: &mut Ibz, x: &Ibz, y: &Ibz) {
    if ibz_is_zero(x) != 0 && ibz_is_zero(y) != 0 {
        ibz_set(d, 1);
        ibz_set(u, 1);
        ibz_set(v, 0);
        return;
    }
    let x1 = x.clone();
    let y1 = y.clone();

    ibz_xgcd(d, u, v, &x1, &y1);

    if ibz_is_zero(u) != 0 {
        if ibz_is_zero(&x1) == 0 {
            let mut y_local = y1.clone();
            if ibz_is_zero(&y_local) != 0 {
                ibz_set(&mut y_local, 1);
            }
            let mut q = Ibz::zero();
            let mut r = Ibz::zero();
            ibz_div(&mut q, &mut r, &x1, &y_local);
            assert!(ibz_is_zero(&r) != 0);
            let mut t = Ibz::zero();
            ibz_sub(&mut t, v, &q);
            *v = t;
        }
        ibz_set(u, 1);
    }

    if ibz_is_zero(&x1) == 0 {
        assert!(ibz_cmp(d, &ibz_const_zero()) > 0);
        let mut r = Ibz::zero();
        ibz_mul(&mut r, &x1, &y1);
        let neg = ibz_cmp(&r, &ibz_const_zero()) < 0;
        let mut q = Ibz::zero();
        ibz_mul(&mut q, &x1, u);
        while ibz_cmp(&q, &ibz_const_zero()) <= 0 {
            let mut q2 = Ibz::zero();
            let mut r2 = Ibz::zero();
            ibz_div(&mut q2, &mut r2, &y1, d);
            assert!(ibz_is_zero(&r2) != 0);
            if neg {
                let mut t = Ibz::zero();
                ibz_neg(&mut t, &q2);
                q2 = t;
            }
            let mut t = Ibz::zero();
            ibz_add(&mut t, u, &q2);
            *u = t;

            let mut q3 = Ibz::zero();
            let mut r3 = Ibz::zero();
            ibz_div(&mut q3, &mut r3, &x1, d);
            assert!(ibz_is_zero(&r3) != 0);
            if neg {
                let mut t = Ibz::zero();
                ibz_neg(&mut t, &q3);
                q3 = t;
            }
            let mut t = Ibz::zero();
            ibz_sub(&mut t, v, &q3);
            *v = t;

            ibz_mul(&mut q, &x1, u);
        }
    }
}

/// `ibz_vec_4_linear_combination_mod(lc, ca, va, cb, vb, mod)`:
/// `lc = (ca*va + cb*vb) mod mod`, centred residue.
pub fn ibz_vec_4_linear_combination_mod(
    lc: &mut IbzVec4,
    coeff_a: &Ibz,
    vec_a: &IbzVec4,
    coeff_b: &Ibz,
    vec_b: &IbzVec4,
    mod_: &Ibz,
) {
    let mut sums = ibz_vec_4_new();
    let mut prod = Ibz::zero();
    for i in 0..4 {
        ibz_mul(&mut sums[i], coeff_a, &vec_a[i]);
        ibz_mul(&mut prod, coeff_b, &vec_b[i]);
        let mut t = Ibz::zero();
        ibz_add(&mut t, &sums[i], &prod);
        sums[i] = t;
        let mut r = Ibz::zero();
        ibz_centered_mod(&mut r, &sums[i], mod_);
        sums[i] = r;
    }
    for i in 0..4 {
        lc[i] = sums[i].clone();
    }
}

/// `ibz_vec_4_copy_mod(res, vec, mod)`: centred residues componentwise.
pub fn ibz_vec_4_copy_mod(res: &mut IbzVec4, vec: &IbzVec4, mod_: &Ibz) {
    for i in 0..4 {
        let mut r = Ibz::zero();
        ibz_centered_mod(&mut r, &vec[i], mod_);
        res[i] = r;
    }
}

/// `ibz_vec_4_scalar_mul_mod(prod, scalar, vec, mod)`: positive residues.
pub fn ibz_vec_4_scalar_mul_mod(prod: &mut IbzVec4, scalar: &Ibz, vec: &IbzVec4, mod_: &Ibz) {
    for i in 0..4 {
        let mut t = Ibz::zero();
        ibz_mul(&mut t, &vec[i], scalar);
        prod[i] = t;
        let mut r = Ibz::zero();
        ibz_mod(&mut r, &prod[i], mod_);
        prod[i] = r;
    }
}

/// `ibz_mat_4x4_is_hnf(mat)`: checks the (upper-triangular, positive
/// diagonal, reduced off-diagonals) HNF predicate.
///
/// Strict transcription of the C reference's quirky implementation,
/// including the slightly odd `linestart < i` test on the column scan.
pub fn ibz_mat_4x4_is_hnf(mat: &IbzMat4x4) -> i32 {
    let mut res = 1;
    let zero = Ibz::zero();
    let mut found;
    let mut ind = 0usize;
    for i in 0..4 {
        for j in 0..i {
            if ibz_is_zero(&mat[i][j]) == 0 {
                res = 0;
            }
        }
        found = 0;
        for j in i..4 {
            if found != 0 {
                if !(ibz_cmp(&mat[i][j], &zero) >= 0) {
                    res = 0;
                }
                if !(ibz_cmp(&mat[i][ind], &mat[i][j]) > 0) {
                    res = 0;
                }
            } else if ibz_is_zero(&mat[i][j]) == 0 {
                found = 1;
                ind = j;
                if !(ibz_cmp(&mat[i][j], &zero) > 0) {
                    res = 0;
                }
            }
        }
    }
    let linestart: i32 = -1;
    let mut i = 0i32;
    for j in 0..4 {
        while i < 4 && ibz_is_zero(&mat[i as usize][j]) != 0 {
            i += 1;
        }
        if i != 4 && !(linestart < i) {
            res = 0;
        }
        i = 0;
    }
    res
}

/// `ibz_mat_4xn_hnf_mod_core(hnf, n, generators, mod)`: compute the HNF
/// of an n-generator lattice mod a multiple of its determinant.
///
/// Strict transcription of Cohen 2.4.8.
pub fn ibz_mat_4xn_hnf_mod_core(
    hnf: &mut IbzMat4x4,
    generator_number: usize,
    generators: &[IbzVec4],
    mod_: &Ibz,
) {
    assert!(generator_number > 3);
    assert!(generators.len() >= generator_number);
    let n = generator_number;

    let mut a: Vec<IbzVec4> = Vec::with_capacity(n);
    for h in 0..n {
        let mut v = ibz_vec_4_new();
        for k in 0..4 {
            v[k] = generators[h][k].clone();
        }
        a.push(v);
    }
    let mut w: [IbzVec4; 4] = [
        ibz_vec_4_new(),
        ibz_vec_4_new(),
        ibz_vec_4_new(),
        ibz_vec_4_new(),
    ];

    assert!(ibz_cmp(mod_, &ibz_const_zero()) > 0);
    let mut m = mod_.clone();

    let mut i: i32 = 3;
    let mut k = (n as i32) - 1;
    let mut j = k;
    while i != -1 {
        while j != 0 {
            j -= 1;
            if ibz_is_zero(&a[j as usize][i as usize]) == 0 {
                let mut d = Ibz::zero();
                let mut u = Ibz::zero();
                let mut v = Ibz::zero();
                let ak_i = a[k as usize][i as usize].clone();
                let aj_i = a[j as usize][i as usize].clone();
                ibz_xgcd_with_u_not_0(&mut d, &mut u, &mut v, &ak_i, &aj_i);

                let mut c_ = ibz_vec_4_new();
                let ak = clone_vec(&a[k as usize]);
                let aj = clone_vec(&a[j as usize]);
                ibz_vec_4_linear_combination(&mut c_, &u, &ak, &v, &aj);

                let mut coeff_1 = Ibz::zero();
                let mut r = Ibz::zero();
                ibz_div(&mut coeff_1, &mut r, &ak_i, &d);
                let mut coeff_2_pre = Ibz::zero();
                ibz_div(&mut coeff_2_pre, &mut r, &aj_i, &d);
                let mut coeff_2 = Ibz::zero();
                ibz_neg(&mut coeff_2, &coeff_2_pre);

                let aj_clone = clone_vec(&a[j as usize]);
                let ak_clone = clone_vec(&a[k as usize]);
                let mut new_aj = ibz_vec_4_new();
                ibz_vec_4_linear_combination_mod(
                    &mut new_aj,
                    &coeff_1,
                    &aj_clone,
                    &coeff_2,
                    &ak_clone,
                    &m,
                );
                a[j as usize] = new_aj;
                let mut new_ak = ibz_vec_4_new();
                ibz_vec_4_copy_mod(&mut new_ak, &c_, &m);
                a[k as usize] = new_ak;
            }
        }

        let mut d = Ibz::zero();
        let mut u = Ibz::zero();
        let mut v = Ibz::zero();
        let aki = a[k as usize][i as usize].clone();
        ibz_xgcd_with_u_not_0(&mut d, &mut u, &mut v, &aki, &m);

        let ak_clone = clone_vec(&a[k as usize]);
        let mut new_wi = ibz_vec_4_new();
        ibz_vec_4_scalar_mul_mod(&mut new_wi, &u, &ak_clone, &m);
        w[i as usize] = new_wi;
        if ibz_is_zero(&w[i as usize][i as usize]) != 0 {
            w[i as usize][i as usize] = m.clone();
        }

        for h in (i + 1)..4 {
            let mut q = Ibz::zero();
            let mut r = Ibz::zero();
            ibz_div_floor(
                &mut q,
                &mut r,
                &w[h as usize][i as usize],
                &w[i as usize][i as usize],
            );
            let mut nq = Ibz::zero();
            ibz_neg(&mut nq, &q);
            let wh = clone_vec(&w[h as usize]);
            let wi = clone_vec(&w[i as usize]);
            let mut new_wh = ibz_vec_4_new();
            ibz_vec_4_linear_combination(&mut new_wh, &ibz_const_one(), &wh, &nq, &wi);
            w[h as usize] = new_wh;
        }

        let mut new_m = Ibz::zero();
        let mut r = Ibz::zero();
        ibz_div(&mut new_m, &mut r, &m, &d);
        assert!(ibz_is_zero(&r) != 0);
        m = new_m;

        if i != 0 {
            k -= 1;
            i -= 1;
            j = k;
            if ibz_is_zero(&a[k as usize][i as usize]) != 0 {
                a[k as usize][i as usize] = m.clone();
            }
        } else {
            k -= 1;
            i -= 1;
            j = k;
        }
    }

    for jj in 0..4 {
        for ii in 0..4 {
            hnf[ii][jj] = w[jj][ii].clone();
        }
    }
}

fn clone_vec(v: &IbzVec4) -> IbzVec4 {
    [v[0].clone(), v[1].clone(), v[2].clone(), v[3].clone()]
}
