//! Left ideals and orders of quaternion algebras.
//!
//! Mirrors `the-sqisign/src/quaternion/ref/generic/ideal.c`.
//! `quat_lideal_t` carries a lattice, a norm, and a borrowed reference
//! to its parent (maximal) order.

use crate::algebra::{quat_alg_norm, QuatAlg, QuatAlgElem};
use crate::dim4::{
    ibz_mat_4x4_identity, ibz_mat_4x4_mul, ibz_mat_4x4_new, ibz_mat_4x4_scalar_mul,
    ibz_mat_4x4_transpose, IbzMat4x4,
};
use crate::ibz::{
    ibz_cmp, ibz_const_two, ibz_div, ibz_divides, ibz_is_one, ibz_is_zero, ibz_mul, ibz_sqrt, Ibz,
};
use crate::lattice::{
    quat_lattice_add, quat_lattice_alg_elem_mul, quat_lattice_conjugate_without_hnf,
    quat_lattice_equal, quat_lattice_gram, quat_lattice_index, quat_lattice_intersect,
    quat_lattice_mul, quat_lattice_reduce_denom, QuatLattice,
};

/// `quat_left_ideal_t`: a left ideal of a maximal order, with its norm.
///
/// We omit the borrowed `parent_order` pointer field from the C ref;
/// every public function that needs the parent order takes it as an
/// explicit argument or expects the caller to keep it alongside the
/// ideal. This avoids lifetime parameters polluting the type at the
/// price of a small amount of bookkeeping on the caller's side.
#[derive(Clone, Debug, Default)]
pub struct QuatLeftIdeal {
    pub lattice: QuatLattice,
    pub norm: Ibz,
}

impl QuatLeftIdeal {
    pub fn new() -> Self {
        Self {
            lattice: QuatLattice::new(),
            norm: Ibz::zero(),
        }
    }
}

/// `quat_lideal_norm(lideal, parent_order)`: compute and set the norm.
pub fn quat_lideal_norm(lideal: &mut QuatLeftIdeal, parent_order: &QuatLattice) {
    quat_lattice_index(&mut lideal.norm, &lideal.lattice, parent_order);
    let saved = lideal.norm.clone();
    let ok = ibz_sqrt(&mut lideal.norm, &saved);
    assert!(ok != 0, "quat_lideal_norm: index must be a perfect square");
}

/// `quat_lideal_copy(copy, copied)`.
pub fn quat_lideal_copy(copy: &mut QuatLeftIdeal, copied: &QuatLeftIdeal) {
    copy.norm = copied.norm.clone();
    copy.lattice = copied.lattice.clone();
}

/// `quat_lideal_create_principal(lideal, x, order, alg)`.
pub fn quat_lideal_create_principal(
    lideal: &mut QuatLeftIdeal,
    x: &QuatAlgElem,
    order: &QuatLattice,
    alg: &QuatAlg,
) {
    quat_lattice_alg_elem_mul(&mut lideal.lattice, order, x, alg);
    let cloned = lideal.lattice.clone();
    quat_lattice_reduce_denom(&mut lideal.lattice, &cloned);
    let mut norm_n = Ibz::zero();
    let mut norm_d = Ibz::zero();
    quat_alg_norm(&mut norm_n, &mut norm_d, x, alg);
    assert!(ibz_is_one(&norm_d) != 0);
    lideal.norm = norm_n;
}

/// `quat_lideal_create(lideal, x, N, order, alg)`.
pub fn quat_lideal_create(
    lideal: &mut QuatLeftIdeal,
    x: &QuatAlgElem,
    big_n: &Ibz,
    order: &QuatLattice,
    alg: &QuatAlg,
) {
    let mut on = QuatLattice::new();
    quat_lideal_create_principal(lideal, x, order, alg);
    ibz_mat_4x4_scalar_mul(&mut on.basis, big_n, &order.basis);
    on.denom = order.denom.clone();
    let lat_clone = lideal.lattice.clone();
    quat_lattice_add(&mut lideal.lattice, &lat_clone, &on);
    quat_lideal_norm(lideal, order);
}

/// `quat_lideal_mul(product, lideal, alpha, alg, parent_order)`.
pub fn quat_lideal_mul(
    product: &mut QuatLeftIdeal,
    lideal: &QuatLeftIdeal,
    alpha: &QuatAlgElem,
    alg: &QuatAlg,
) {
    quat_lattice_alg_elem_mul(&mut product.lattice, &lideal.lattice, alpha, alg);
    let mut norm = Ibz::zero();
    let mut norm_d = Ibz::zero();
    quat_alg_norm(&mut norm, &mut norm_d, alpha, alg);
    ibz_mul(&mut product.norm, &lideal.norm, &norm);
    assert!(ibz_divides(&product.norm, &norm_d) != 0);
    let saved = product.norm.clone();
    let mut r = Ibz::zero();
    ibz_div(&mut product.norm, &mut r, &saved, &norm_d);
}

/// `quat_lideal_add(sum, I1, I2, alg, parent_order)`.
pub fn quat_lideal_add(
    sum: &mut QuatLeftIdeal,
    i1: &QuatLeftIdeal,
    i2: &QuatLeftIdeal,
    parent_order: &QuatLattice,
) {
    quat_lattice_add(&mut sum.lattice, &i1.lattice, &i2.lattice);
    quat_lideal_norm(sum, parent_order);
}

/// `quat_lideal_inter(inter, I1, I2, parent_order)`.
pub fn quat_lideal_inter(
    inter: &mut QuatLeftIdeal,
    i1: &QuatLeftIdeal,
    i2: &QuatLeftIdeal,
    parent_order: &QuatLattice,
) {
    quat_lattice_intersect(&mut inter.lattice, &i1.lattice, &i2.lattice);
    quat_lideal_norm(inter, parent_order);
}

/// `quat_lideal_equals(I1, I2, alg)`: equality assuming same parent order.
pub fn quat_lideal_equals(i1: &QuatLeftIdeal, i2: &QuatLeftIdeal) -> i32 {
    let norm_eq = ibz_cmp(&i1.norm, &i2.norm) == 0;
    let lat_eq = quat_lattice_equal(&i1.lattice, &i2.lattice) != 0;
    if norm_eq && lat_eq {
        1
    } else {
        0
    }
}

/// `quat_lideal_inverse_lattice_without_hnf(inv, lideal, alg)`.
pub fn quat_lideal_inverse_lattice_without_hnf(inv: &mut QuatLattice, lideal: &QuatLeftIdeal) {
    quat_lattice_conjugate_without_hnf(inv, &lideal.lattice);
    let saved = inv.denom.clone();
    ibz_mul(&mut inv.denom, &saved, &lideal.norm);
}

/// `quat_lideal_right_transporter(trans, l1, l2, alg)`.
pub fn quat_lideal_right_transporter(
    trans: &mut QuatLattice,
    lideal1: &QuatLeftIdeal,
    lideal2: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    let mut inv = QuatLattice::new();
    quat_lideal_inverse_lattice_without_hnf(&mut inv, lideal1);
    quat_lattice_mul(trans, &inv, &lideal2.lattice, alg);
}

/// `quat_lideal_right_order(order, lideal, alg)`.
pub fn quat_lideal_right_order(order: &mut QuatLattice, lideal: &QuatLeftIdeal, alg: &QuatAlg) {
    quat_lideal_right_transporter(order, lideal, lideal, alg);
}

/// `quat_lideal_class_gram(G, lideal, alg)`.
pub fn quat_lideal_class_gram(g: &mut IbzMat4x4, lideal: &QuatLeftIdeal, alg: &QuatAlg) {
    quat_lattice_gram(g, &lideal.lattice, alg);
    let mut divisor = Ibz::zero();
    ibz_mul(&mut divisor, &lideal.lattice.denom, &lideal.lattice.denom);
    let cloned = divisor.clone();
    ibz_mul(&mut divisor, &cloned, &lideal.norm);

    let mut rmd = Ibz::zero();
    for i in 0..4 {
        for j in 0..=i {
            let saved = g[i][j].clone();
            ibz_div(&mut g[i][j], &mut rmd, &saved, &divisor);
            assert!(ibz_is_zero(&rmd) != 0);
        }
    }
    for i in 0..4 {
        for j in 0..i {
            g[j][i] = g[i][j].clone();
        }
    }
}

/// `quat_lideal_conjugate_without_hnf(conj, new_parent_order, lideal, alg)`.
pub fn quat_lideal_conjugate_without_hnf(
    conj: &mut QuatLeftIdeal,
    new_parent_order: &mut QuatLattice,
    lideal: &QuatLeftIdeal,
    alg: &QuatAlg,
) {
    quat_lideal_right_order(new_parent_order, lideal, alg);
    quat_lattice_conjugate_without_hnf(&mut conj.lattice, &lideal.lattice);
    conj.norm = lideal.norm.clone();
}

/// `quat_lideal_generator(gen, lideal, alg)`: search for a non-scalar
/// generator. Deterministic by construction (no RNG). Mirrors the C
/// reference's nested loop verbatim.
///
/// Returns 1 on success, 0 on failure (in practice the C code loops
/// forever on failure; we cap the outer norm at 100 to avoid hangs in
/// pathological cases the differential vectors do not exercise).
pub fn quat_lideal_generator(gen: &mut QuatAlgElem, lideal: &QuatLeftIdeal, alg: &QuatAlg) -> i32 {
    use crate::dim4::{ibz_mat_4x4_eval, ibz_vec_4_content, ibz_vec_4_set};
    use crate::ibz::{ibz_const_one, ibz_gcd};
    let mut vec = [Ibz::zero(), Ibz::zero(), Ibz::zero(), Ibz::zero()];
    let mut gcd_v = Ibz::zero();
    let mut norm_int = Ibz::zero();
    let mut norm_denom = Ibz::zero();
    let mut q = Ibz::zero();
    let mut r = Ibz::zero();
    let mut int_norm = 0i32;
    while int_norm < 100 {
        int_norm += 1;
        for a in -int_norm..=int_norm {
            for b in -(int_norm - a.abs())..=(int_norm - a.abs()) {
                for c in -(int_norm - a.abs() - b.abs())..=(int_norm - a.abs() - b.abs()) {
                    let d = int_norm - a.abs() - b.abs() - c.abs();
                    ibz_vec_4_set(&mut vec, a, b, c, d);
                    ibz_vec_4_content(&mut gcd_v, &vec);
                    if ibz_is_one(&gcd_v) != 0 {
                        ibz_mat_4x4_eval(&mut gen.coord, &lideal.lattice.basis, &vec);
                        gen.denom = lideal.lattice.denom.clone();
                        quat_alg_norm(&mut norm_int, &mut norm_denom, gen, alg);
                        assert!(ibz_is_one(&norm_denom) != 0);
                        ibz_div(&mut q, &mut r, &norm_int, &lideal.norm);
                        assert!(ibz_is_zero(&r) != 0);
                        let mut g2 = Ibz::zero();
                        ibz_gcd(&mut g2, &lideal.norm, &q);
                        let found = if ibz_cmp(&g2, &ibz_const_one()) == 0 {
                            1
                        } else {
                            0
                        };
                        if found != 0 {
                            return 1;
                        }
                    }
                }
            }
        }
    }
    0
}

/// `quat_order_discriminant(disc, order, alg)`.
pub fn quat_order_discriminant(disc: &mut Ibz, order: &QuatLattice, alg: &QuatAlg) -> i32 {
    let mut det = Ibz::zero();
    let mut sqr = Ibz::zero();
    let mut div = Ibz::zero();
    let mut transposed = ibz_mat_4x4_new();
    let mut norm = ibz_mat_4x4_new();
    let mut prod = ibz_mat_4x4_new();
    ibz_mat_4x4_transpose(&mut transposed, &order.basis);
    ibz_mat_4x4_identity(&mut norm);
    norm[2][2] = alg.p.clone();
    norm[3][3] = alg.p.clone();
    let norm_clone = clone_mat(&norm);
    ibz_mat_4x4_scalar_mul(&mut norm, &ibz_const_two(), &norm_clone);
    ibz_mat_4x4_mul(&mut prod, &transposed, &norm);
    let prod_clone = clone_mat(&prod);
    ibz_mat_4x4_mul(&mut prod, &prod_clone, &order.basis);
    let _ = crate::dim4::ibz_mat_4x4_inv_with_det_as_denom(None, Some(&mut det), &prod);
    ibz_mul(&mut div, &order.denom, &order.denom);
    let div_clone = div.clone();
    ibz_mul(&mut div, &div_clone, &div_clone);
    let div_clone = div.clone();
    ibz_mul(&mut div, &div_clone, &div_clone);
    let mut div_rem = Ibz::zero();
    ibz_div(&mut sqr, &mut div_rem, &det, &div);
    let mut ok = ibz_is_zero(&div_rem);
    let s = ibz_sqrt(disc, &sqr);
    ok &= s;
    ok
}

/// `quat_order_is_maximal(order, alg)`.
pub fn quat_order_is_maximal(order: &QuatLattice, alg: &QuatAlg) -> i32 {
    let mut disc = Ibz::zero();
    quat_order_discriminant(&mut disc, order, alg);
    if ibz_cmp(&disc, &alg.p) == 0 {
        1
    } else {
        0
    }
}

fn clone_mat(m: &IbzMat4x4) -> IbzMat4x4 {
    let mut out = ibz_mat_4x4_new();
    for i in 0..4 {
        for j in 0..4 {
            out[i][j] = m[i][j].clone();
        }
    }
    out
}
