//! SQIsign `quaternion`: orders and ideals.
//!
//! Mirrors `vendor/the-sqisign/src/quaternion`. The verification subset is
//! ported in **Phase 2, unit 6**; the remainder in **Phase 3, unit 10**.
//!
//! # Phase 2, unit 6: foundations
//!
//! This batch ports the four foundational source files:
//! `intbig.c` (`ibz_t` wrappers around GMP `mpz_t`), `integers.c` (small
//! integer helpers and the deterministic Cornacchia primitive),
//! `algebra.c` (quaternion algebra and elements), and `finit.c` (constructors
//! and destructors for the algebra structs).
//!
//! ## Differential boundary
//!
//! Unlike the `mp`, `gf`, `ec`, and `hd` modules whose `ibz_t`-free
//! word-array boundaries permit raw-memory byte equality with the C
//! reference, `ibz_t` is an `mpz_t` with library-specific limb sizes and
//! layouts. Equivalence here is at the **canonical value level**: the C
//! harness emits each `ibz_t` as a 1-byte sign tag (0x00 non-negative,
//! 0x01 negative) followed by a 4-byte little-endian length prefix N then
//! N bytes of magnitude in big-endian (the natural shape of
//! `mpz_export`); the Rust port deserializes with
//! [`num_bigint::BigInt::from_signed_bytes_be`] composed with the sign tag,
//! computes via `num-bigint`, and reserializes the same way for byte
//! comparison.
//!
//! ## Cryptographic note
//!
//! These arithmetic paths operate on **public** lattice data. Signature
//! secrets do not flow through here: they flow through the `gf` and `ec`
//! layers, which have already been ported with explicit constant-time
//! primitives. Hence `num-bigint`'s non-constant-time semantics are
//! acceptable for the quaternion module. See `~/sqisign-port-notes.md` for
//! the full architectural decision record.
//!
//! ## In scope this batch
//!
//! `intbig.c`: the deterministic arithmetic, predicate, conversion, and
//! number-theory primitives that have no RNG dependency. RNG-driven
//! `ibz_rand_*` are deferred (they would be RNG-dependent vectors, which
//! is a separate harness pattern). String I/O `ibz_print`,
//! `ibz_convert_to_str`, `ibz_set_from_str` are deferred as they have no
//! cryptographic role.
//!
//! `integers.c`: `ibz_cornacchia_prime` only. `ibz_generate_random_prime`
//! depends on RNG output and is deferred.
//!
//! `algebra.c`: the pure-`ibz_t` element primitives: `quat_alg_coord_mul`,
//! `quat_alg_mul`, `quat_alg_conj`, `quat_alg_scalar`, `quat_alg_norm`,
//! `quat_alg_equal_denom`, plus the trivial `*_set`, `*_copy`,
//! `*_mul_by_scalar` helpers. `quat_alg_add`, `quat_alg_sub`,
//! `quat_alg_normalize`, `quat_alg_elem_equal`, `quat_alg_elem_is_zero`,
//! and `quat_alg_make_primitive` are **deferred** because they depend on
//! `ibz_vec_4_*` helpers defined in `dim4.c` (out of scope this batch)
//! and `quat_lattice_contains` defined in `lattice.c` (out of scope this
//! batch).
//!
//! `finit.c`: constructors and destructors. In Rust these collapse into
//! `Default` impls and `Drop`; no differential vectors are needed because
//! they expose no observable computation beyond zero-initialization.

#![forbid(unsafe_code)]
// The quaternion module is a strict transcription of the C reference
// (intbig.c, integers.c, algebra.c, finit.c, dim2.c, dim4.c, lattice.c,
// ideal.c, normeq.c, hnf.c). We preserve the C control flow and the
// quirks of the reference's mpz_t-on-mini-gmp arithmetic; clippy's
// idiomatic-Rust lints would obscure that mirror.
#![allow(clippy::needless_range_loop)]
#![allow(clippy::nonminimal_bool)]
#![allow(clippy::manual_memcpy)]
#![allow(clippy::needless_borrows_for_generic_args)]

pub mod algebra;
pub mod dim2;
pub mod dim4;
pub mod dpe;
pub mod hnf;
pub mod ibz;
pub mod ideal;
pub mod integers;
pub mod lat_ball;
pub mod lattice;
pub mod lll;
pub mod lll_applications;
pub mod lll_verification;
pub mod normeq;
pub mod rationals;

pub use algebra::{
    quat_alg_add, quat_alg_conj, quat_alg_coord_mul, quat_alg_elem_copy, quat_alg_elem_copy_ibz,
    quat_alg_elem_equal, quat_alg_elem_is_zero, quat_alg_elem_mul_by_scalar, quat_alg_elem_set,
    quat_alg_equal_denom, quat_alg_make_primitive, quat_alg_mul, quat_alg_norm, quat_alg_normalize,
    quat_alg_scalar, quat_alg_sub, QuatAlg, QuatAlgElem,
};
pub use dim2::{
    ibz_2x2_mul_mod, ibz_mat_2x2_add, ibz_mat_2x2_copy, ibz_mat_2x2_det_from_ibz, ibz_mat_2x2_eval,
    ibz_mat_2x2_inv_mod, ibz_mat_2x2_new, ibz_mat_2x2_set, ibz_vec_2_new, ibz_vec_2_set, IbzMat2x2,
    IbzVec2,
};
pub use dim4::{
    ibz_inv_dim4_make_coeff_mpm, ibz_inv_dim4_make_coeff_pmp, ibz_mat_4x4_copy, ibz_mat_4x4_equal,
    ibz_mat_4x4_eval, ibz_mat_4x4_eval_t, ibz_mat_4x4_gcd, ibz_mat_4x4_identity,
    ibz_mat_4x4_inv_with_det_as_denom, ibz_mat_4x4_is_identity, ibz_mat_4x4_mul,
    ibz_mat_4x4_negate, ibz_mat_4x4_new, ibz_mat_4x4_scalar_div, ibz_mat_4x4_scalar_mul,
    ibz_mat_4x4_transpose, ibz_mat_4x4_zero, ibz_vec_4_add, ibz_vec_4_content, ibz_vec_4_copy,
    ibz_vec_4_copy_ibz, ibz_vec_4_is_zero, ibz_vec_4_linear_combination, ibz_vec_4_negate,
    ibz_vec_4_new, ibz_vec_4_scalar_div, ibz_vec_4_scalar_mul, ibz_vec_4_set, ibz_vec_4_sub,
    quat_qf_eval, IbzMat4x4, IbzVec4,
};
pub use dpe::{
    dpe_abs, dpe_add, dpe_cmp, dpe_cmp_d, dpe_div, dpe_get_d, dpe_get_z, dpe_mul, dpe_neg,
    dpe_round, dpe_set, dpe_set_d, dpe_set_si, dpe_set_ui, dpe_set_z, dpe_sqrt, dpe_sub,
    dpe_zero_p, Dpe,
};
pub use hnf::{
    ibz_centered_mod, ibz_conditional_assign, ibz_mat_4x4_is_hnf, ibz_mat_4xn_hnf_mod_core,
    ibz_mod_not_zero, ibz_vec_4_copy_mod, ibz_vec_4_linear_combination_mod,
    ibz_vec_4_scalar_mul_mod, ibz_xgcd, ibz_xgcd_with_u_not_0,
};
pub use ibz::{
    ibz_abs, ibz_add, ibz_bitsize, ibz_cmp, ibz_cmp_int32, ibz_const_one, ibz_const_three,
    ibz_const_two, ibz_const_zero, ibz_copy_digits, ibz_div, ibz_div_2exp, ibz_div_floor,
    ibz_divides, ibz_gcd, ibz_get, ibz_invmod, ibz_is_even, ibz_is_odd, ibz_is_one, ibz_is_zero,
    ibz_legendre, ibz_mod, ibz_mod_ui, ibz_mul, ibz_neg, ibz_pow, ibz_pow_mod, ibz_probab_prime,
    ibz_set, ibz_size_in_base, ibz_sqrt, ibz_sqrt_floor, ibz_sqrt_mod_p, ibz_sub, ibz_to_digits,
    ibz_two_adic, Ibz,
};
pub use ideal::{
    quat_lideal_add, quat_lideal_class_gram, quat_lideal_conjugate_without_hnf, quat_lideal_copy,
    quat_lideal_create, quat_lideal_create_principal, quat_lideal_equals, quat_lideal_generator,
    quat_lideal_inter, quat_lideal_inverse_lattice_without_hnf, quat_lideal_mul, quat_lideal_norm,
    quat_lideal_right_order, quat_lideal_right_transporter, quat_order_discriminant,
    quat_order_is_maximal, QuatLeftIdeal,
};
pub use integers::ibz_cornacchia_prime;
pub use lat_ball::quat_lattice_bound_parallelogram;
pub use lattice::{
    quat_lattice_add, quat_lattice_alg_elem_mul, quat_lattice_conjugate_without_hnf,
    quat_lattice_contains, quat_lattice_dual_without_hnf, quat_lattice_equal, quat_lattice_gram,
    quat_lattice_hnf, quat_lattice_inclusion, quat_lattice_index, quat_lattice_intersect,
    quat_lattice_mat_alg_coord_mul_without_hnf, quat_lattice_mul, quat_lattice_reduce_denom,
    QuatLattice,
};
pub use lll::{
    quat_lattice_lll, quat_lll_core, DELTABAR, DELTA_DENOM, DELTA_NUM, EPSILON_DENOM, EPSILON_NUM,
    ETABAR,
};
pub use lll_applications::{quat_lideal_lideal_mul_reduced, quat_lideal_reduce_basis};
pub use lll_verification::{
    ibq_vec_4_copy_ibz, quat_lll_bilinear, quat_lll_gram_schmidt_transposed_with_ibq,
    quat_lll_set_ibq_parameters, quat_lll_verify,
};
pub use normeq::{
    quat_change_to_O0_basis, quat_lattice_O0_set, quat_lattice_O0_set_extremal,
    quat_order_elem_create, QuatPExtremalMaximalOrder,
};
pub use rationals::{
    ibq_abs, ibq_add, ibq_cmp, ibq_copy, ibq_inv, ibq_is_ibz, ibq_is_one, ibq_is_zero,
    ibq_mat_4x4_new, ibq_mul, ibq_neg, ibq_reduce, ibq_set, ibq_sub, ibq_to_ibz, ibq_vec_4_new,
    Ibq, IbqMat4x4, IbqVec4,
};
