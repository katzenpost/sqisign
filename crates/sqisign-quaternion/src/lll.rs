//! L2 (Nguyen-Stehlé) LLL reduction.
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/lll/l2.c`. Reduces
//! a 4x4 lattice basis in place, updating the Gram matrix to reflect the
//! new basis. The implementation is dimension-4 specific because the
//! upstream code is, and the quaternion module relies on that.
//!
//! ## Constants
//!
//! `DELTABAR = 0.995`, `ETABAR = 0.505`. Both match `lll_internals.h`.
//!
//! ## Non-uniqueness
//!
//! Many bases are LLL-reduced for the same input lattice: the algorithm
//! makes deterministic swaps but the output basis is not canonical in any
//! intrinsic way. The differential test against the C reference therefore
//! checks two properties side by side:
//!
//!  * structural validity (Lovász + size-reduction conditions), via
//!    [`crate::lll_verification::quat_lll_verify`];
//!  * bit-exact agreement with the C reference's `dpe` arithmetic, which
//!    is achievable in pure IEEE 754 except for the `dpe_set_z` quirk
//!    that we faithfully mirror in [`crate::dpe::dpe_set_z`].

use crate::algebra::QuatAlg;
use crate::dim4::{ibz_mat_4x4_copy, IbzMat4x4};
use crate::dpe::{
    dpe_cmp, dpe_cmp_d, dpe_div, dpe_get_z, dpe_mul, dpe_round, dpe_set, dpe_set_d, dpe_set_z,
    dpe_sub, Dpe,
};
use crate::ibz::{ibz_mul, ibz_sub, Ibz};
use crate::lattice::{quat_lattice_gram, QuatLattice};

/// `DELTABAR` in `lll_internals.h`.
pub const DELTABAR: f64 = 0.995;
/// `ETABAR` in `lll_internals.h`.
pub const ETABAR: f64 = 0.505;
/// `DELTA_NUM` / `DELTA_DENOM` rational form of `DELTABAR` minus an `ε`.
pub const DELTA_NUM: i32 = 99;
pub const DELTA_DENOM: i32 = 100;
/// `EPSILON_NUM` / `EPSILON_DENOM` for the `eta = 1/2 + ε` parameter.
pub const EPSILON_NUM: i32 = 1;
pub const EPSILON_DENOM: i32 = 100;

/// Helper: index the lower triangle, swapping `(i, j)` if needed so we
/// always read `g[max][min]`. Matches the `SYM` macro in `l2.c`.
#[inline]
fn sym(g: &IbzMat4x4, i: usize, j: usize) -> &Ibz {
    if i < j {
        &g[j][i]
    } else {
        &g[i][j]
    }
}

#[inline]
fn sym_mut(g: &mut IbzMat4x4, i: usize, j: usize) -> &mut Ibz {
    if i < j {
        &mut g[j][i]
    } else {
        &mut g[i][j]
    }
}

/// `quat_lll_core(G, basis)`: in-place L2 reduction. The Gram matrix `G`
/// is updated to reflect the new basis. The input Gram matrix's upper
/// triangle is overwritten by a symmetric copy of the lower triangle at
/// the end of the call (matching the C reference's last loop).
pub fn quat_lll_core(g: &mut IbzMat4x4, basis: &mut IbzMat4x4) {
    let mut delta_bar = Dpe::zero();
    dpe_set_d(&mut delta_bar, DELTABAR);

    // dpe scratch for Gram-Schmidt and Lovász.
    let mut r = [[Dpe::zero(); 4]; 4];
    let mut u = [[Dpe::zero(); 4]; 4];
    let mut lovasz = [Dpe::zero(); 4];

    let mut xf = Dpe::zero();
    let mut tmp_f = Dpe::zero();
    let mut x = Ibz::zero();
    let mut tmp_i = Ibz::zero();

    // Initialize r[0][0] = G[0][0]
    dpe_set_z(&mut r[0][0], &g[0][0]);

    let mut kappa: usize = 1;
    while kappa < 4 {
        // Size-reduce b_kappa
        let mut done = false;
        while !done {
            // Recompute the kappa-th row of the Cholesky factorization.
            for j in 0..=kappa {
                let mut r_kj = Dpe::zero();
                dpe_set_z(&mut r_kj, &g[kappa][j]);
                for k in 0..j {
                    let mut tmp = Dpe::zero();
                    dpe_mul(&mut tmp, &r[kappa][k], &u[j][k]);
                    let mut new_r = Dpe::zero();
                    dpe_sub(&mut new_r, &r_kj, &tmp);
                    r_kj = new_r;
                }
                r[kappa][j] = r_kj;
                if j < kappa {
                    let mut new_u = Dpe::zero();
                    dpe_div(&mut new_u, &r[kappa][j], &r[j][j]);
                    u[kappa][j] = new_u;
                }
            }

            done = true;
            // Size-reduce (descending i)
            for i in (0..kappa).rev() {
                if dpe_cmp_d(&u[kappa][i], ETABAR) > 0 || dpe_cmp_d(&u[kappa][i], -ETABAR) < 0 {
                    done = false;
                    dpe_set(&mut xf, &u[kappa][i]);
                    let xf_in = xf;
                    dpe_round(&mut xf, &xf_in);
                    dpe_get_z(&mut x, &xf);
                    // Update basis: b_kappa -= X * b_i
                    for j in 0..4 {
                        ibz_mul(&mut tmp_i, &x, &basis[j][i]);
                        let old = basis[j][kappa].clone();
                        ibz_sub(&mut basis[j][kappa], &old, &tmp_i);
                    }
                    // Update lower half of the Gram matrix.
                    // <b_kappa, b_kappa> -= X * <b_kappa, b_i>
                    ibz_mul(&mut tmp_i, &x, &g[kappa][i]);
                    let old = g[kappa][kappa].clone();
                    ibz_sub(&mut g[kappa][kappa], &old, &tmp_i);
                    // For each j: <b_kappa, b_j> -= X * <b_i, b_j>
                    for j in 0..4 {
                        let sym_ij = sym(g, i, j).clone();
                        ibz_mul(&mut tmp_i, &x, &sym_ij);
                        let cell = sym_mut(g, kappa, j);
                        let old = cell.clone();
                        ibz_sub(cell, &old, &tmp_i);
                    }
                    // Update u[kappa][j] for j < i
                    for j in 0..i {
                        let mut prod = Dpe::zero();
                        dpe_mul(&mut prod, &xf, &u[i][j]);
                        let mut new_u = Dpe::zero();
                        dpe_sub(&mut new_u, &u[kappa][j], &prod);
                        u[kappa][j] = new_u;
                    }
                }
            }
        }

        // Check Lovász' conditions.
        dpe_set_z(&mut lovasz[0], &g[kappa][kappa]);
        for i in 1..kappa {
            dpe_mul(&mut tmp_f, &u[kappa][i - 1], &r[kappa][i - 1]);
            let prev = lovasz[i - 1];
            let mut new_lov = Dpe::zero();
            dpe_sub(&mut new_lov, &prev, &tmp_f);
            lovasz[i] = new_lov;
        }
        let mut swap = kappa;
        while swap > 0 {
            dpe_mul(&mut tmp_f, &delta_bar, &r[swap - 1][swap - 1]);
            if dpe_cmp(&tmp_f, &lovasz[swap - 1]) < 0 {
                break;
            }
            swap -= 1;
        }

        if kappa != swap {
            // Insert b_kappa before b_swap in the basis and lower half Gram.
            let mut j = kappa;
            while j > swap {
                for i in 0..4 {
                    basis[i].swap(j, j - 1);
                    if i == j - 1 {
                        // swap diagonal entries g[i][i] and g[j][j]
                        // g[i][i] is g[j-1][j-1], g[j][j] is g[j][j]
                        // The C reference does ibz_swap(&G[i][i], &G[j][j]).
                        // For i == j-1, this swaps g[j-1][j-1] <-> g[j][j].
                        let a = g[j - 1][j - 1].clone();
                        let b = g[j][j].clone();
                        g[j - 1][j - 1] = b;
                        g[j][j] = a;
                    } else if i != j {
                        // swap SYM(G, i, j) with SYM(G, i, j-1)
                        let a = sym(g, i, j).clone();
                        let b = sym(g, i, j - 1).clone();
                        *sym_mut(g, i, j) = b;
                        *sym_mut(g, i, j - 1) = a;
                    }
                }
                j -= 1;
            }
            // Copy u[kappa][i] and r[kappa][i] into row `swap`, for i < swap.
            for i in 0..swap {
                let uki = u[kappa][i];
                let rki = r[kappa][i];
                u[swap][i] = uki;
                r[swap][i] = rki;
            }
            // r[swap][swap] <- lovasz[swap]
            let ls = lovasz[swap];
            r[swap][swap] = ls;
            kappa = swap;
        }

        kappa += 1;
    }

    // Fill in the upper half of the Gram matrix.
    for i in 0..4 {
        for j in (i + 1)..4 {
            g[i][j] = g[j][i].clone();
        }
    }
}

/// `quat_lattice_lll(red, lattice, alg)`: compute the lattice Gram matrix
/// and run [`quat_lll_core`] on a copy of the basis, returning the reduced
/// basis in `red`. Always returns `0` (the C reference does the same; the
/// return code is reserved for future error reporting).
pub fn quat_lattice_lll(red: &mut IbzMat4x4, lattice: &QuatLattice, alg: &QuatAlg) -> i32 {
    let mut g = crate::dim4::ibz_mat_4x4_new();
    quat_lattice_gram(&mut g, lattice, alg);
    ibz_mat_4x4_copy(red, &lattice.basis);
    quat_lll_core(&mut g, red);
    0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dim4::{ibz_mat_4x4_new, IbzMat4x4};
    use crate::ibz::ibz_set;

    fn identity() -> IbzMat4x4 {
        let mut m = ibz_mat_4x4_new();
        for i in 0..4 {
            ibz_set(&mut m[i][i], 1);
        }
        m
    }

    #[test]
    fn lll_on_identity_is_idempotent() {
        // For the identity Gram matrix and identity basis, no reduction
        // is necessary; the basis must come back unchanged (and the Gram
        // matrix likewise).
        let mut g = identity();
        let mut basis = identity();
        let g_before = g.clone();
        let basis_before = basis.clone();
        quat_lll_core(&mut g, &mut basis);
        assert_eq!(g, g_before, "Gram changed on identity input");
        assert_eq!(basis, basis_before, "basis changed on identity input");
    }
}
