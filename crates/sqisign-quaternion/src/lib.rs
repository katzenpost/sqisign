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

mod algebra;
mod ibz;
mod integers;

pub use algebra::{
    quat_alg_conj, quat_alg_coord_mul, quat_alg_elem_copy, quat_alg_elem_copy_ibz,
    quat_alg_elem_mul_by_scalar, quat_alg_elem_set, quat_alg_equal_denom, quat_alg_mul,
    quat_alg_norm, quat_alg_scalar, QuatAlg, QuatAlgElem,
};
pub use ibz::{
    ibz_abs, ibz_add, ibz_bitsize, ibz_cmp, ibz_cmp_int32, ibz_const_one, ibz_const_three,
    ibz_const_two, ibz_const_zero, ibz_copy_digits, ibz_div, ibz_div_2exp, ibz_div_floor,
    ibz_divides, ibz_gcd, ibz_get, ibz_invmod, ibz_is_even, ibz_is_odd, ibz_is_one, ibz_is_zero,
    ibz_legendre, ibz_mod, ibz_mod_ui, ibz_mul, ibz_neg, ibz_pow, ibz_pow_mod, ibz_probab_prime,
    ibz_set, ibz_size_in_base, ibz_sqrt, ibz_sqrt_floor, ibz_sqrt_mod_p, ibz_sub, ibz_to_digits,
    ibz_two_adic, Ibz,
};
pub use integers::ibz_cornacchia_prime;
