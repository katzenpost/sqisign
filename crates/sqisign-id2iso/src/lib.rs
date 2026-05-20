//! SQIsign `id2iso`: ideal-to-isogeny translation.
//!
//! Mirrors `vendor/the-sqisign/src/id2iso`. Ported in **Phase 3, unit 11**.
//!
//! # Scope
//!
//! Two source files, ported in full:
//! * `id2iso.c` (~338 lines) - kernel/ideal translation helpers.
//! * `dim2id2iso.c` (~1172 lines) - the Clapotis-style `find_uv` enumerator
//!   and the dimension-two Kani-diagram ideal-to-isogeny driver.
//!
//! # Determinism and differential boundaries
//!
//! The deterministic helpers carry per-boundary differential batteries
//! recorded by the `cdump` harness against the pinned C reference:
//!
//! * [`ec_biscalar_mul_ibz_vec`]
//! * [`id2iso_ideal_to_kernel_dlogs_even`]
//! * [`matrix_application_even_basis`]
//! * [`endomorphism_application_even_basis`]
//! * [`id2iso_kernel_dlogs_to_ideal_even`]
//! * [`change_of_basis_matrix_tate`] and [`change_of_basis_matrix_tate_invert`]
//!
//! The remaining entry points sit downstream of [`quat_represent_integer`],
//! which is RNG-driven and itself **deferred** in the quaternion crate
//! (see `normeq.rs`). Those entry points (`find_uv`,
//! `fixed_degree_isogeny_and_eval`, `dim2id2iso_ideal_to_isogeny_clapotis`,
//! `dim2id2iso_arbitrary_isogeny_evaluation`) are ported in full and
//! compile, but their differential vectors are deferred until the
//! quaternion RNG path is wired:
//!
//! * `find_uv` is deterministic *given* an input ideal, but the C reference
//!   test constructs that input via `quat_represent_integer`. Without a
//!   handcrafted reproducer (the search space is sensitive to ideal class,
//!   not arbitrary handcrafted norms), per-boundary vectors are deferred.
//! * `fixed_degree_isogeny_and_eval` and the two `dim2id2iso_*` drivers
//!   call `quat_represent_integer` directly. They are RNG-driven by
//!   construction.

#![forbid(unsafe_code)]
// Strict transcription of the C reference: preserve loop and naming
// idioms so a side-by-side review reads as a 1:1 port.
#![allow(non_snake_case)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::needless_borrows_for_generic_args)]

pub mod dim2id2iso;
pub mod id2iso;

pub use dim2id2iso::{
    dim2id2iso_arbitrary_isogeny_evaluation, dim2id2iso_ideal_to_isogeny_clapotis, find_uv,
    fixed_degree_isogeny_and_eval,
};
pub use id2iso::{
    change_of_basis_matrix_tate, change_of_basis_matrix_tate_invert, ec_biscalar_mul_ibz_vec,
    endomorphism_application_even_basis, id2iso_ideal_to_kernel_dlogs_even,
    id2iso_kernel_dlogs_to_ideal_even, matrix_application_even_basis,
};
