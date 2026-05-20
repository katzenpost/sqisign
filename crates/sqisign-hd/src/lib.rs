//! SQIsign `hd`: (2,2)-isogenies in the theta model.
//!
//! Mirrors `vendor/the-sqisign/src/hd`. The module operates on abelian
//! surfaces (dimension 2): points are pairs of [`EcPoint`]s, curves are
//! pairs of [`EcCurve`]s, and isogenies are computed in the theta model.
//! All arithmetic flows through [`sqisign_gf`] and [`sqisign_ec`].
//!
//! Ported sources (Phase 2 unit 8):
//!
//! - `hd/ref/lvlx/hd.c` (top-level couple-point and Jacobian helpers).
//! - `hd/ref/lvlx/theta_structure.c` (theta-point precomputation,
//!   pointwise square, Hadamard transform, doubling).
//! - `hd/ref/lvlx/theta_isogenies.c` (the (2,2)-isogeny chain and its
//!   building blocks: gluing, intermediate steps, splitting).
//! - `precomp/ref/lvl1/hd_splitting_transforms.c` (the precomputed
//!   splitting and normalisation transforms used during splitting).
//!
//! ## Differential boundaries exposed for testing
//!
//! The public boundaries are the three chain entry points and a handful
//! of trivially exercised helpers. The intermediate static helpers in
//! `theta_isogenies.c` are not exposed as separate boundaries since they
//! are deeply intertwined with `_theta_chain_compute_impl`; their
//! correctness is witnessed by the chain vectors.
//!
//! ## Deferred boundaries
//!
//! - `theta_chain_compute_and_eval_randomized` is **not ported in this
//!   unit**. The reference calls `randombytes` via `sample_random_index`
//!   to choose one of six normalisation matrices for the splitting step.
//!   That makes the boundary caller-visible non-deterministic, so a
//!   differential vector battery is impossible without a caller-supplied
//!   RNG handle. The current `sqisign-common` crate exposes the CTR-DRBG
//!   only via [`sqisign_common::CtrDrbg::fill`] (no top-level
//!   `randombytes`), so wiring this faithfully requires a follow-up that
//!   thread-locals or argument-passes the DRBG. The deterministic
//!   variants ([`theta_chain_compute_and_eval`] and
//!   [`theta_chain_compute_and_eval_verify`]) cover all paths the
//!   reference reaches with `randomize=false`.
//!
//! ## Ambiguous decisions
//!
//! - `theta_chain_compute_impl` stack-allocates two variable-length arrays
//!   (`uint16_t todo[space]` and `theta_couple_jac_point_t jacQ1[space]`)
//!   plus an opaque `theta_point_t pts[numP ? numP : 1]`. Rust has no
//!   safe VLAs, so we use heap [`Vec`] allocations of the same lengths.
//!   The boundary is unchanged (no per-element observable behaviour).
//! - `test_point_order_twof` and `test_jac_order_twof` are `static inline`
//!   helpers in `ec.h` not yet ported into [`sqisign_ec`]. They appear
//!   only inside debug-only sanity prints in `gluing_compute`, with no
//!   effect on the boundary output. We inline tiny equivalents below so
//!   the port stays self-contained without expanding the `ec` surface.

#![forbid(unsafe_code)]
#![allow(non_snake_case)]
// The reference uses explicit `for (j = 0; j < numP; ++j) { ... pts[j] ... }`
// loops in the chain-impl translation; reshaping to iterators makes the
// transcription harder to audit against the C side, so we keep them as
// index loops.
#![allow(clippy::needless_range_loop)]

use sqisign_ec::{
    add, copy_curve, copy_point, dbl, dblw, ec_curve_init, ec_dbl, ec_dbl_iter, ec_is_equal,
    ec_is_zero, jac_from_ws, jac_to_ws, jac_to_xz, jac_to_xz_add_components, lift_basis,
    AddComponents, EcBasis, EcCurve, EcPoint, JacPoint,
};
use sqisign_gf::{
    fp2_add, fp2_batched_inv, fp2_copy, fp2_is_equal, fp2_is_zero, fp2_mul, fp2_neg, fp2_select,
    fp2_set_one, fp2_set_zero, fp2_sqr, fp2_sqrt, fp2_sub, Fp2, NWORDS_FIELD,
};

// ---------------------------------------------------------------------------
// hd.h: type definitions for dimension 2
// ---------------------------------------------------------------------------

/// Exponent of additional torsion required when the kernel is supplied
/// in `E[2^(n+HD_extra_torsion)]`. Mirrors `#define HD_extra_torsion 2`.
pub const HD_EXTRA_TORSION: u32 = 2;

/// Couple point on an elliptic product, in `(X : Z)` coordinates.
/// Mirrors `theta_couple_point_t {ec_point_t P1; ec_point_t P2;}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaCouplePoint {
    pub p1: EcPoint,
    pub p2: EcPoint,
}

impl ThetaCouplePoint {
    pub const fn zero() -> Self {
        ThetaCouplePoint {
            p1: EcPoint::zero(),
            p2: EcPoint::zero(),
        }
    }
}

/// Triple of couple points `(T1, T2, T1 - T2)`. Mirrors
/// `theta_kernel_couple_points_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaKernelCouplePoints {
    pub t1: ThetaCouplePoint,
    pub t2: ThetaCouplePoint,
    pub t1m2: ThetaCouplePoint,
}

impl ThetaKernelCouplePoints {
    pub const fn zero() -> Self {
        ThetaKernelCouplePoints {
            t1: ThetaCouplePoint::zero(),
            t2: ThetaCouplePoint::zero(),
            t1m2: ThetaCouplePoint::zero(),
        }
    }
}

/// Couple point on an elliptic product, in `(X : Y : Z)` Jacobian
/// coordinates. Mirrors `theta_couple_jac_point_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaCoupleJacPoint {
    pub p1: JacPoint,
    pub p2: JacPoint,
}

impl ThetaCoupleJacPoint {
    pub const fn zero() -> Self {
        ThetaCoupleJacPoint {
            p1: JacPoint::zero(),
            p2: JacPoint::zero(),
        }
    }
}

/// An elliptic product `E1 x E2`. Mirrors `theta_couple_curve_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaCoupleCurve {
    pub e1: EcCurve,
    pub e2: EcCurve,
}

impl ThetaCoupleCurve {
    pub const fn zero() -> Self {
        ThetaCoupleCurve {
            e1: EcCurve::zero(),
            e2: EcCurve::zero(),
        }
    }
}

/// A theta point: four-tuple of fp2 coordinates. Mirrors
/// `theta_point_t {fp2_t x; fp2_t y; fp2_t z; fp2_t t;}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaPoint {
    pub x: Fp2,
    pub y: Fp2,
    pub z: Fp2,
    pub t: Fp2,
}

impl ThetaPoint {
    pub const fn zero() -> Self {
        let z = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        ThetaPoint {
            x: z,
            y: z,
            z,
            t: z,
        }
    }
}

/// Compact theta point with repeating components. Mirrors
/// `theta_point_compact_t {fp2_t x; fp2_t y;}`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaPointCompact {
    pub x: Fp2,
    pub y: Fp2,
}

impl ThetaPointCompact {
    pub const fn zero() -> Self {
        let z = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        ThetaPointCompact { x: z, y: z }
    }
}

/// Theta structure: null point plus the eight precomputed factors used
/// in doubling and (2,2)-isogenies. Mirrors `theta_structure_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ThetaStructure {
    pub null_point: ThetaPoint,
    pub precomputation: bool,
    pub xyz_big_0: Fp2,
    pub yzt_big_0: Fp2,
    pub xzt_big_0: Fp2,
    pub xyt_big_0: Fp2,
    pub xyz_0: Fp2,
    pub yzt_0: Fp2,
    pub xzt_0: Fp2,
    pub xyt_0: Fp2,
}

impl ThetaStructure {
    pub const fn zero() -> Self {
        let z = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        ThetaStructure {
            null_point: ThetaPoint::zero(),
            precomputation: false,
            xyz_big_0: z,
            yzt_big_0: z,
            xzt_big_0: z,
            xyt_big_0: z,
            xyz_0: z,
            yzt_0: z,
            xzt_0: z,
            xyt_0: z,
        }
    }
}

/// 2x2 translation matrix used while computing a compatible theta
/// structure during gluing. Mirrors `translation_matrix_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TranslationMatrix {
    pub g00: Fp2,
    pub g01: Fp2,
    pub g10: Fp2,
    pub g11: Fp2,
}

impl TranslationMatrix {
    pub const fn zero() -> Self {
        let z = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        TranslationMatrix {
            g00: z,
            g01: z,
            g10: z,
            g11: z,
        }
    }
}

/// 4x4 matrix used for basis changes. Mirrors `basis_change_matrix_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct BasisChangeMatrix {
    pub m: [[Fp2; 4]; 4],
}

impl BasisChangeMatrix {
    pub const fn zero() -> Self {
        let z = Fp2 {
            re: [0u64; NWORDS_FIELD],
            im: [0u64; NWORDS_FIELD],
        };
        BasisChangeMatrix { m: [[z; 4]; 4] }
    }
}

/// State of a gluing (2,2)-isogeny step. Mirrors `theta_gluing_t`.
#[derive(Clone, Copy, Debug)]
pub struct ThetaGluing {
    pub domain: ThetaCoupleCurve,
    pub xy_k1_8: ThetaCoupleJacPoint,
    pub image_k1_8: ThetaPointCompact,
    pub mat: BasisChangeMatrix,
    pub precomputation: ThetaPoint,
    pub codomain: ThetaPoint,
}

impl ThetaGluing {
    pub const fn zero() -> Self {
        ThetaGluing {
            domain: ThetaCoupleCurve::zero(),
            xy_k1_8: ThetaCoupleJacPoint::zero(),
            image_k1_8: ThetaPointCompact::zero(),
            mat: BasisChangeMatrix::zero(),
            precomputation: ThetaPoint::zero(),
            codomain: ThetaPoint::zero(),
        }
    }
}

/// State of a standard (2,2)-isogeny step. Mirrors `theta_isogeny_t`.
#[derive(Clone, Copy, Debug)]
pub struct ThetaIsogeny {
    pub t1_8: ThetaPoint,
    pub t2_8: ThetaPoint,
    pub hadamard_bool_1: bool,
    pub hadamard_bool_2: bool,
    pub domain: ThetaStructure,
    pub precomputation: ThetaPoint,
    pub codomain: ThetaStructure,
}

impl ThetaIsogeny {
    pub const fn zero() -> Self {
        ThetaIsogeny {
            t1_8: ThetaPoint::zero(),
            t2_8: ThetaPoint::zero(),
            hadamard_bool_1: false,
            hadamard_bool_2: false,
            domain: ThetaStructure::zero(),
            precomputation: ThetaPoint::zero(),
            codomain: ThetaStructure::zero(),
        }
    }
}

/// State of the final splitting isomorphism. Mirrors
/// `theta_splitting_t`.
#[derive(Clone, Copy, Debug)]
pub struct ThetaSplitting {
    pub mat: BasisChangeMatrix,
    pub b: ThetaStructure,
}

impl ThetaSplitting {
    pub const fn zero() -> Self {
        ThetaSplitting {
            mat: BasisChangeMatrix::zero(),
            b: ThetaStructure::zero(),
        }
    }
}

// ---------------------------------------------------------------------------
// hd_splitting_transforms.c (lvl1 constants used by splitting)
// ---------------------------------------------------------------------------

/// Indices of the ten even theta-characteristic pairs. Mirrors the lvl1
/// constant `EVEN_INDEX`.
pub const EVEN_INDEX: [[u32; 2]; 10] = [
    [0, 0],
    [0, 1],
    [0, 2],
    [0, 3],
    [1, 0],
    [1, 2],
    [2, 0],
    [2, 1],
    [3, 0],
    [3, 3],
];

/// Character evaluation table `CHI_EVAL[a][i] in {-1, +1}`. Mirrors the
/// lvl1 constant `CHI_EVAL`; the reference uses `int` and relies on
/// `(int >> 1)` of `-1` producing `-1` (so cast to `uint32_t` yields
/// `0xFFFFFFFF`), which we preserve.
pub const CHI_EVAL: [[i32; 4]; 4] = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]];

const FP2_ZERO_IDX: u8 = 0;
const FP2_ONE_IDX: u8 = 1;
const FP2_I_IDX: u8 = 2;
const FP2_MINUS_ONE_IDX: u8 = 3;
const FP2_MINUS_I_IDX: u8 = 4;

/// `FP2_CONSTANTS[k]` is `{0, 1, i, -1, -i}` represented as a level-1
/// Montgomery fp2. The reference's lvl1 source spells these limb arrays
/// out for each radix; we construct them from `fp2_set_zero`,
/// `fp2_set_one`, and `fp2_neg`, which produces byte-identical limbs at
/// our pinned 64-bit, non-Broadwell build.
fn fp2_constants() -> [Fp2; 5] {
    let mut zero = Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    };
    fp2_set_zero(&mut zero);
    let mut one = zero;
    fp2_set_one(&mut one);
    let mut i_const = zero;
    // i = 0 + 1 * i in Fp2; the imaginary part is Montgomery one.
    i_const.im = one.re;
    let mut minus_one = zero;
    fp2_neg(&mut minus_one, &one);
    let mut minus_i = zero;
    fp2_neg(&mut minus_i, &i_const);
    [zero, one, i_const, minus_one, minus_i]
}

/// Precomputed basis-change matrix: each cell is an index into
/// `FP2_CONSTANTS`. Mirrors `precomp_basis_change_matrix_t`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PrecompBasisChangeMatrix {
    pub m: [[u8; 4]; 4],
}

/// `SPLITTING_TRANSFORMS[i]` for `i in 0..10`: the ten precomputed
/// matrices `splitting_compute` selects between. Mirrors the lvl1
/// constant of the same name.
pub const SPLITTING_TRANSFORMS: [PrecompBasisChangeMatrix; 10] = [
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_I_IDX, FP2_ONE_IDX, FP2_I_IDX],
            [FP2_ONE_IDX, FP2_MINUS_I_IDX, FP2_MINUS_ONE_IDX, FP2_I_IDX],
            [FP2_ONE_IDX, FP2_I_IDX, FP2_MINUS_ONE_IDX, FP2_MINUS_I_IDX],
            [FP2_MINUS_ONE_IDX, FP2_I_IDX, FP2_MINUS_ONE_IDX, FP2_I_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_MINUS_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_MINUS_ONE_IDX, FP2_ZERO_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_MINUS_ONE_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
            [
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_ONE_IDX,
            ],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX],
        ],
    },
];

/// `NORMALIZATION_TRANSFORMS[i]` for `i in 0..6`: the six matrices the
/// randomised splitting variant chooses between. Mirrors the lvl1
/// constant of the same name.
pub const NORMALIZATION_TRANSFORMS: [PrecompBasisChangeMatrix; 6] = [
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX],
            [FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX],
            [FP2_ZERO_IDX, FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
            [FP2_ONE_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX, FP2_ZERO_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
            ],
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
            [
                FP2_MINUS_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_ONE_IDX,
            ],
            [
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
                FP2_MINUS_ONE_IDX,
                FP2_ONE_IDX,
            ],
            [FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX, FP2_ONE_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_MINUS_ONE_IDX, FP2_I_IDX, FP2_I_IDX, FP2_ONE_IDX],
            [FP2_I_IDX, FP2_MINUS_ONE_IDX, FP2_ONE_IDX, FP2_I_IDX],
            [FP2_I_IDX, FP2_ONE_IDX, FP2_MINUS_ONE_IDX, FP2_I_IDX],
            [FP2_ONE_IDX, FP2_I_IDX, FP2_I_IDX, FP2_MINUS_ONE_IDX],
        ],
    },
    PrecompBasisChangeMatrix {
        m: [
            [FP2_ONE_IDX, FP2_I_IDX, FP2_I_IDX, FP2_MINUS_ONE_IDX],
            [FP2_I_IDX, FP2_ONE_IDX, FP2_MINUS_ONE_IDX, FP2_I_IDX],
            [FP2_I_IDX, FP2_MINUS_ONE_IDX, FP2_ONE_IDX, FP2_I_IDX],
            [FP2_MINUS_ONE_IDX, FP2_I_IDX, FP2_I_IDX, FP2_ONE_IDX],
        ],
    },
];

// ---------------------------------------------------------------------------
// hd.c: couple-point and Jacobian-couple helpers
// ---------------------------------------------------------------------------

#[inline]
fn fp2_zero() -> Fp2 {
    Fp2 {
        re: [0u64; NWORDS_FIELD],
        im: [0u64; NWORDS_FIELD],
    }
}

/// Mirrors `double_couple_point`: doubles both legs in parallel on the
/// elliptic product.
pub fn double_couple_point(
    out: &mut ThetaCouplePoint,
    inp: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) {
    let in_p1 = inp.p1;
    let in_p2 = inp.p2;
    ec_dbl(&mut out.p1, &in_p1, &e12.e1);
    ec_dbl(&mut out.p2, &in_p2, &e12.e2);
}

/// Mirrors `double_couple_point_iter`: doubles `n` times.
pub fn double_couple_point_iter(
    out: &mut ThetaCouplePoint,
    n: u32,
    inp: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) {
    if n == 0 {
        *out = *inp;
    } else {
        double_couple_point(out, inp, e12);
        for _ in 0..(n - 1) {
            let cur = *out;
            double_couple_point(out, &cur, e12);
        }
    }
}

/// Mirrors `add_couple_jac_points`: parallel `ADD` on both legs.
pub fn add_couple_jac_points(
    out: &mut ThetaCoupleJacPoint,
    t1: &ThetaCoupleJacPoint,
    t2: &ThetaCoupleJacPoint,
    e12: &ThetaCoupleCurve,
) {
    let t1p1 = t1.p1;
    let t1p2 = t1.p2;
    let t2p1 = t2.p1;
    let t2p2 = t2.p2;
    add(&mut out.p1, &t1p1, &t2p1, &e12.e1);
    add(&mut out.p2, &t1p2, &t2p2, &e12.e2);
}

/// Mirrors `double_couple_jac_point`: parallel `DBL` on both legs.
pub fn double_couple_jac_point(
    out: &mut ThetaCoupleJacPoint,
    inp: &ThetaCoupleJacPoint,
    e12: &ThetaCoupleCurve,
) {
    let in_p1 = inp.p1;
    let in_p2 = inp.p2;
    dbl(&mut out.p1, &in_p1, &e12.e1);
    dbl(&mut out.p2, &in_p2, &e12.e2);
}

/// Mirrors `double_couple_jac_point_iter`: doubles `n` times, using the
/// modified Jacobian (Weierstrass-scratch) form for `n >= 2`.
pub fn double_couple_jac_point_iter(
    out: &mut ThetaCoupleJacPoint,
    n: u32,
    inp: &ThetaCoupleJacPoint,
    e12: &ThetaCoupleCurve,
) {
    if n == 0 {
        *out = *inp;
    } else if n == 1 {
        double_couple_jac_point(out, inp, e12);
    } else {
        let mut a1 = fp2_zero();
        let mut a2 = fp2_zero();
        let mut t1 = fp2_zero();
        let mut t2 = fp2_zero();

        jac_to_ws(&mut out.p1, &mut t1, &mut a1, &inp.p1, &e12.e1);
        jac_to_ws(&mut out.p2, &mut t2, &mut a2, &inp.p2, &e12.e2);

        let in_p1 = out.p1;
        let in_t1 = t1;
        dblw(&mut out.p1, &mut t1, &in_p1, &in_t1);
        let in_p2 = out.p2;
        let in_t2 = t2;
        dblw(&mut out.p2, &mut t2, &in_p2, &in_t2);

        for _ in 0..(n - 1) {
            let in_p1 = out.p1;
            let in_t1 = t1;
            dblw(&mut out.p1, &mut t1, &in_p1, &in_t1);
            let in_p2 = out.p2;
            let in_t2 = t2;
            dblw(&mut out.p2, &mut t2, &in_p2, &in_t2);
        }

        let in_p1 = out.p1;
        jac_from_ws(&mut out.p1, &in_p1, &a1, &e12.e1);
        let in_p2 = out.p2;
        jac_from_ws(&mut out.p2, &in_p2, &a2, &e12.e2);
    }
}

/// Mirrors `couple_jac_to_xz`: forgetful map `(X : Y : Z) -> (X : Z)`
/// on both legs.
pub fn couple_jac_to_xz(p: &mut ThetaCouplePoint, xy_p: &ThetaCoupleJacPoint) {
    jac_to_xz(&mut p.p1, &xy_p.p1);
    jac_to_xz(&mut p.p2, &xy_p.p2);
}

/// Mirrors `copy_bases_to_kernel`: lays out `(B1, B2)` as
/// `(T1, T2, T1 - T2)` on the elliptic product.
pub fn copy_bases_to_kernel(ker: &mut ThetaKernelCouplePoints, b1: &EcBasis, b2: &EcBasis) {
    copy_point(&mut ker.t1.p1, &b1.P);
    copy_point(&mut ker.t2.p1, &b1.Q);
    copy_point(&mut ker.t1m2.p1, &b1.PmQ);

    copy_point(&mut ker.t1.p2, &b2.P);
    copy_point(&mut ker.t2.p2, &b2.Q);
    copy_point(&mut ker.t1m2.p2, &b2.PmQ);
}

/// Local port of the `static inline test_point_order_twof` in `ec.h`.
/// Used only for debug-only sanity checks inside `gluing_compute`.
fn test_point_order_twof(p: &EcPoint, e: &EcCurve, t: i32) -> i32 {
    let mut test = *p;
    let mut curve = EcCurve::zero();
    copy_curve(&mut curve, e);

    if ec_is_zero(&test) != 0 {
        return 0;
    }
    let cur = test;
    ec_dbl_iter(&mut test, t - 1, &cur, &mut curve);
    if ec_is_zero(&test) != 0 {
        return 0;
    }
    let cur = test;
    ec_dbl(&mut test, &cur, &curve);
    ec_is_zero(&test) as i32
}

/// Local port of the `static inline test_jac_order_twof` in `ec.h`.
/// Used only for debug-only sanity checks inside `gluing_compute`.
fn test_jac_order_twof(p: &JacPoint, e: &EcCurve, t: i32) -> i32 {
    let mut test = *p;
    if fp2_is_zero(&test.z) != 0 {
        return 0;
    }
    for _ in 0..(t - 1) {
        let cur = test;
        dbl(&mut test, &cur, e);
    }
    if fp2_is_zero(&test.z) != 0 {
        return 0;
    }
    let cur = test;
    dbl(&mut test, &cur, e);
    fp2_is_zero(&test.z) as i32
}

/// Mirrors the inline `test_couple_point_order_twof` in `hd.h`: AND of
/// the per-leg `test_point_order_twof` checks.
pub fn test_couple_point_order_twof(t: &ThetaCouplePoint, e: &ThetaCoupleCurve, twof: i32) -> i32 {
    let check_p1 = test_point_order_twof(&t.p1, &e.e1, twof);
    let check_p2 = test_point_order_twof(&t.p2, &e.e2, twof);
    check_p1 & check_p2
}

// ---------------------------------------------------------------------------
// theta_structure.c: theta-point primitives
// ---------------------------------------------------------------------------

/// Mirrors `hadamard`: in = (x,y,z,t),
/// out = (x+y+z+t, x-y+z-t, x+y-z-t, x-y-z+t).
pub fn hadamard(out: &mut ThetaPoint, inp: &ThetaPoint) {
    let mut t1 = fp2_zero();
    let mut t2 = fp2_zero();
    let mut t3 = fp2_zero();
    let mut t4 = fp2_zero();
    fp2_add(&mut t1, &inp.x, &inp.y);
    fp2_sub(&mut t2, &inp.x, &inp.y);
    fp2_add(&mut t3, &inp.z, &inp.t);
    fp2_sub(&mut t4, &inp.z, &inp.t);
    fp2_add(&mut out.x, &t1, &t3);
    fp2_add(&mut out.y, &t2, &t4);
    fp2_sub(&mut out.z, &t1, &t3);
    fp2_sub(&mut out.t, &t2, &t4);
}

/// Mirrors `pointwise_square`: squares each coordinate.
pub fn pointwise_square(out: &mut ThetaPoint, inp: &ThetaPoint) {
    fp2_sqr(&mut out.x, &inp.x);
    fp2_sqr(&mut out.y, &inp.y);
    fp2_sqr(&mut out.z, &inp.z);
    fp2_sqr(&mut out.t, &inp.t);
}

/// Mirrors `to_squared_theta`: pointwise square then Hadamard.
pub fn to_squared_theta(out: &mut ThetaPoint, inp: &ThetaPoint) {
    pointwise_square(out, inp);
    let cur = *out;
    hadamard(out, &cur);
}

/// Mirrors `theta_precomputation`: caches the eight projective factors
/// used in subsequent doublings and (2,2)-isogenies. Idempotent.
pub fn theta_precomputation(a: &mut ThetaStructure) {
    if a.precomputation {
        return;
    }

    let mut a_dual = ThetaPoint::zero();
    to_squared_theta(&mut a_dual, &a.null_point);

    let mut t1 = fp2_zero();
    let mut t2 = fp2_zero();
    fp2_mul(&mut t1, &a_dual.x, &a_dual.y);
    fp2_mul(&mut t2, &a_dual.z, &a_dual.t);
    fp2_mul(&mut a.xyz_big_0, &t1, &a_dual.z);
    fp2_mul(&mut a.xyt_big_0, &t1, &a_dual.t);
    fp2_mul(&mut a.yzt_big_0, &t2, &a_dual.y);
    fp2_mul(&mut a.xzt_big_0, &t2, &a_dual.x);

    let np = a.null_point;
    fp2_mul(&mut t1, &np.x, &np.y);
    fp2_mul(&mut t2, &np.z, &np.t);
    fp2_mul(&mut a.xyz_0, &t1, &np.z);
    fp2_mul(&mut a.xyt_0, &t1, &np.t);
    fp2_mul(&mut a.yzt_0, &t2, &np.y);
    fp2_mul(&mut a.xzt_0, &t2, &np.x);

    a.precomputation = true;
}

/// Mirrors `double_point`: doubles `in` on the theta structure `A`.
pub fn double_point(out: &mut ThetaPoint, a: &mut ThetaStructure, inp: &ThetaPoint) {
    to_squared_theta(out, inp);
    let cur = *out;
    fp2_sqr(&mut out.x, &cur.x);
    fp2_sqr(&mut out.y, &cur.y);
    fp2_sqr(&mut out.z, &cur.z);
    fp2_sqr(&mut out.t, &cur.t);

    if !a.precomputation {
        theta_precomputation(a);
    }
    let cur = *out;
    fp2_mul(&mut out.x, &cur.x, &a.yzt_big_0);
    fp2_mul(&mut out.y, &cur.y, &a.xzt_big_0);
    fp2_mul(&mut out.z, &cur.z, &a.xyt_big_0);
    fp2_mul(&mut out.t, &cur.t, &a.xyz_big_0);

    let cur = *out;
    hadamard(out, &cur);

    let cur = *out;
    fp2_mul(&mut out.x, &cur.x, &a.yzt_0);
    fp2_mul(&mut out.y, &cur.y, &a.xzt_0);
    fp2_mul(&mut out.z, &cur.z, &a.xyt_0);
    fp2_mul(&mut out.t, &cur.t, &a.xyz_0);
}

/// Mirrors `double_iter`: doubles `exp` times.
pub fn double_iter(out: &mut ThetaPoint, a: &mut ThetaStructure, inp: &ThetaPoint, exp: i32) {
    if exp == 0 {
        *out = *inp;
    } else {
        double_point(out, a, inp);
        for _ in 1..exp {
            let cur = *out;
            double_point(out, a, &cur);
        }
    }
}

/// Mirrors `is_product_theta_point`: returns `0xFFFFFFFF` if
/// `P.x * P.t == P.y * P.z`, else zero.
pub fn is_product_theta_point(p: &ThetaPoint) -> u32 {
    let mut t1 = fp2_zero();
    let mut t2 = fp2_zero();
    fp2_mul(&mut t1, &p.x, &p.t);
    fp2_mul(&mut t2, &p.y, &p.z);
    fp2_is_equal(&t1, &t2)
}

// ---------------------------------------------------------------------------
// theta_isogenies.c: helpers (static in C, internal here)
// ---------------------------------------------------------------------------

fn select_base_change_matrix(
    m: &mut BasisChangeMatrix,
    m1: &BasisChangeMatrix,
    m2: &PrecompBasisChangeMatrix,
    option: u32,
) {
    let consts = fp2_constants();
    for i in 0..4 {
        for j in 0..4 {
            fp2_select(
                &mut m.m[i][j],
                &m1.m[i][j],
                &consts[m2.m[i][j] as usize],
                option,
            );
        }
    }
}

fn set_base_change_matrix_from_precomp(res: &mut BasisChangeMatrix, m: &PrecompBasisChangeMatrix) {
    let consts = fp2_constants();
    for i in 0..4 {
        for j in 0..4 {
            res.m[i][j] = consts[m.m[i][j] as usize];
        }
    }
}

fn choose_index_theta_point(res: &mut Fp2, ind: i32, t: &ThetaPoint) {
    let src = match ind.rem_euclid(4) {
        0 => &t.x,
        1 => &t.y,
        2 => &t.z,
        3 => &t.t,
        _ => unreachable!(),
    };
    fp2_copy(res, src);
}

fn apply_isomorphism_general(
    res: &mut ThetaPoint,
    m: &BasisChangeMatrix,
    p: &ThetaPoint,
    pt_not_zero: bool,
) {
    let mut x1 = fp2_zero();
    let mut temp = ThetaPoint::zero();

    fp2_mul(&mut temp.x, &p.x, &m.m[0][0]);
    fp2_mul(&mut x1, &p.y, &m.m[0][1]);
    let cur = temp.x;
    fp2_add(&mut temp.x, &cur, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[0][2]);
    let cur = temp.x;
    fp2_add(&mut temp.x, &cur, &x1);

    fp2_mul(&mut temp.y, &p.x, &m.m[1][0]);
    fp2_mul(&mut x1, &p.y, &m.m[1][1]);
    let cur = temp.y;
    fp2_add(&mut temp.y, &cur, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[1][2]);
    let cur = temp.y;
    fp2_add(&mut temp.y, &cur, &x1);

    fp2_mul(&mut temp.z, &p.x, &m.m[2][0]);
    fp2_mul(&mut x1, &p.y, &m.m[2][1]);
    let cur = temp.z;
    fp2_add(&mut temp.z, &cur, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[2][2]);
    let cur = temp.z;
    fp2_add(&mut temp.z, &cur, &x1);

    fp2_mul(&mut temp.t, &p.x, &m.m[3][0]);
    fp2_mul(&mut x1, &p.y, &m.m[3][1]);
    let cur = temp.t;
    fp2_add(&mut temp.t, &cur, &x1);
    fp2_mul(&mut x1, &p.z, &m.m[3][2]);
    let cur = temp.t;
    fp2_add(&mut temp.t, &cur, &x1);

    if pt_not_zero {
        fp2_mul(&mut x1, &p.t, &m.m[0][3]);
        let cur = temp.x;
        fp2_add(&mut temp.x, &cur, &x1);

        fp2_mul(&mut x1, &p.t, &m.m[1][3]);
        let cur = temp.y;
        fp2_add(&mut temp.y, &cur, &x1);

        fp2_mul(&mut x1, &p.t, &m.m[2][3]);
        let cur = temp.z;
        fp2_add(&mut temp.z, &cur, &x1);

        fp2_mul(&mut x1, &p.t, &m.m[3][3]);
        let cur = temp.t;
        fp2_add(&mut temp.t, &cur, &x1);
    }

    fp2_copy(&mut res.x, &temp.x);
    fp2_copy(&mut res.y, &temp.y);
    fp2_copy(&mut res.z, &temp.z);
    fp2_copy(&mut res.t, &temp.t);
}

fn apply_isomorphism(res: &mut ThetaPoint, m: &BasisChangeMatrix, p: &ThetaPoint) {
    apply_isomorphism_general(res, m, p, true);
}

fn base_change_matrix_multiplication(
    res: &mut BasisChangeMatrix,
    m1: &BasisChangeMatrix,
    m2: &BasisChangeMatrix,
) {
    let mut tmp = BasisChangeMatrix::zero();
    let mut sum = fp2_zero();
    for i in 0..4 {
        for j in 0..4 {
            fp2_set_zero(&mut sum);
            for k in 0..4 {
                let mut m_ik = m1.m[i][k];
                let m_kj = m2.m[k][j];
                let cur = m_ik;
                fp2_mul(&mut m_ik, &cur, &m_kj);
                let cur = sum;
                fp2_add(&mut sum, &cur, &m_ik);
            }
            tmp.m[i][j] = sum;
        }
    }
    *res = tmp;
}

fn base_change(out: &mut ThetaPoint, phi: &ThetaGluing, t: &ThetaCouplePoint) {
    let mut null_point = ThetaPoint::zero();
    fp2_mul(&mut null_point.x, &t.p1.x, &t.p2.x);
    fp2_mul(&mut null_point.y, &t.p1.x, &t.p2.z);
    fp2_mul(&mut null_point.z, &t.p2.x, &t.p1.z);
    fp2_mul(&mut null_point.t, &t.p1.z, &t.p2.z);
    apply_isomorphism(out, &phi.mat, &null_point);
}

fn action_by_translation_z_and_det(z_inv: &mut Fp2, det_inv: &mut Fp2, p4: &EcPoint, p2: &EcPoint) {
    fp2_copy(z_inv, &p4.z);
    let mut tmp = fp2_zero();
    fp2_mul(det_inv, &p4.x, &p2.z);
    fp2_mul(&mut tmp, &p4.z, &p2.x);
    let cur = *det_inv;
    fp2_sub(det_inv, &cur, &tmp);
}

fn action_by_translation_compute_matrix(
    g: &mut TranslationMatrix,
    p4: &EcPoint,
    p2: &EcPoint,
    z_inv: &Fp2,
    det_inv: &Fp2,
) {
    let mut tmp = fp2_zero();

    fp2_mul(&mut tmp, &p4.x, z_inv);
    fp2_mul(&mut g.g10, &p4.x, &p2.x);
    let cur = g.g10;
    fp2_mul(&mut g.g10, &cur, det_inv);
    let cur = g.g10;
    fp2_sub(&mut g.g10, &cur, &tmp);

    fp2_mul(&mut g.g11, &p2.x, det_inv);
    let cur = g.g11;
    fp2_mul(&mut g.g11, &cur, &p4.z);

    fp2_neg(&mut g.g00, &g.g11);

    fp2_mul(&mut g.g01, &p2.z, det_inv);
    let cur = g.g01;
    fp2_mul(&mut g.g01, &cur, &p4.z);
    let cur = g.g01;
    fp2_neg(&mut g.g01, &cur);
}

fn verify_two_torsion(
    k1_2: &ThetaCouplePoint,
    k2_2: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) -> i32 {
    if (ec_is_zero(&k1_2.p1) | ec_is_zero(&k1_2.p2) | ec_is_zero(&k2_2.p1) | ec_is_zero(&k2_2.p2))
        != 0
    {
        return 0;
    }

    if (ec_is_equal(&k1_2.p1, &k2_2.p1) | ec_is_equal(&k1_2.p2, &k2_2.p2)) != 0 {
        return 0;
    }

    let mut o1 = ThetaCouplePoint::zero();
    let mut o2 = ThetaCouplePoint::zero();
    double_couple_point(&mut o1, k1_2, e12);
    double_couple_point(&mut o2, k2_2, e12);
    if (ec_is_zero(&o1.p1) & ec_is_zero(&o1.p2) & ec_is_zero(&o2.p1) & ec_is_zero(&o2.p2)) == 0 {
        return 0;
    }

    1
}

fn action_by_translation(
    gi: &mut [TranslationMatrix; 4],
    k1_4: &ThetaCouplePoint,
    k2_4: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) -> i32 {
    let mut k1_2 = ThetaCouplePoint::zero();
    let mut k2_2 = ThetaCouplePoint::zero();
    double_couple_point(&mut k1_2, k1_4, e12);
    double_couple_point(&mut k2_2, k2_4, e12);

    if verify_two_torsion(&k1_2, &k2_2, e12) == 0 {
        return 0;
    }

    let mut inverses = [fp2_zero(); 8];
    let (a, b) = inverses.split_at_mut(4);
    action_by_translation_z_and_det(&mut a[0], &mut b[0], &k1_4.p1, &k1_2.p1);
    action_by_translation_z_and_det(&mut a[1], &mut b[1], &k1_4.p2, &k1_2.p2);
    action_by_translation_z_and_det(&mut a[2], &mut b[2], &k2_4.p1, &k2_2.p1);
    action_by_translation_z_and_det(&mut a[3], &mut b[3], &k2_4.p2, &k2_2.p2);

    fp2_batched_inv(&mut inverses);
    if fp2_is_zero(&inverses[0]) != 0 {
        return 0;
    }

    action_by_translation_compute_matrix(
        &mut gi[0],
        &k1_4.p1,
        &k1_2.p1,
        &inverses[0],
        &inverses[4],
    );
    action_by_translation_compute_matrix(
        &mut gi[1],
        &k1_4.p2,
        &k1_2.p2,
        &inverses[1],
        &inverses[5],
    );
    action_by_translation_compute_matrix(
        &mut gi[2],
        &k2_4.p1,
        &k2_2.p1,
        &inverses[2],
        &inverses[6],
    );
    action_by_translation_compute_matrix(
        &mut gi[3],
        &k2_4.p2,
        &k2_2.p2,
        &inverses[3],
        &inverses[7],
    );

    1
}

fn gluing_change_of_basis(
    m: &mut BasisChangeMatrix,
    k1_4: &ThetaCouplePoint,
    k2_4: &ThetaCouplePoint,
    e12: &ThetaCoupleCurve,
) -> i32 {
    let mut gi = [TranslationMatrix::zero(); 4];
    if action_by_translation(&mut gi, k1_4, k2_4, e12) == 0 {
        return 0;
    }

    let mut t001 = fp2_zero();
    let mut t101 = fp2_zero();
    let mut t002 = fp2_zero();
    let mut t102 = fp2_zero();
    let mut tmp = fp2_zero();

    fp2_mul(&mut t001, &gi[0].g00, &gi[2].g00);
    fp2_mul(&mut tmp, &gi[0].g01, &gi[2].g10);
    let cur = t001;
    fp2_add(&mut t001, &cur, &tmp);

    fp2_mul(&mut t101, &gi[0].g10, &gi[2].g00);
    fp2_mul(&mut tmp, &gi[0].g11, &gi[2].g10);
    let cur = t101;
    fp2_add(&mut t101, &cur, &tmp);

    fp2_mul(&mut t002, &gi[1].g00, &gi[3].g00);
    fp2_mul(&mut tmp, &gi[1].g01, &gi[3].g10);
    let cur = t002;
    fp2_add(&mut t002, &cur, &tmp);

    fp2_mul(&mut t102, &gi[1].g10, &gi[3].g00);
    fp2_mul(&mut tmp, &gi[1].g11, &gi[3].g10);
    let cur = t102;
    fp2_add(&mut t102, &cur, &tmp);

    fp2_set_one(&mut m.m[0][0]);
    fp2_mul(&mut tmp, &t001, &t002);
    let cur = m.m[0][0];
    fp2_add(&mut m.m[0][0], &cur, &tmp);
    fp2_mul(&mut tmp, &gi[2].g00, &gi[3].g00);
    let cur = m.m[0][0];
    fp2_add(&mut m.m[0][0], &cur, &tmp);
    fp2_mul(&mut tmp, &gi[0].g00, &gi[1].g00);
    let cur = m.m[0][0];
    fp2_add(&mut m.m[0][0], &cur, &tmp);

    fp2_mul(&mut m.m[0][1], &t001, &t102);
    fp2_mul(&mut tmp, &gi[2].g00, &gi[3].g10);
    let cur = m.m[0][1];
    fp2_add(&mut m.m[0][1], &cur, &tmp);
    fp2_mul(&mut tmp, &gi[0].g00, &gi[1].g10);
    let cur = m.m[0][1];
    fp2_add(&mut m.m[0][1], &cur, &tmp);

    fp2_mul(&mut m.m[0][2], &t101, &t002);
    fp2_mul(&mut tmp, &gi[2].g10, &gi[3].g00);
    let cur = m.m[0][2];
    fp2_add(&mut m.m[0][2], &cur, &tmp);
    fp2_mul(&mut tmp, &gi[0].g10, &gi[1].g00);
    let cur = m.m[0][2];
    fp2_add(&mut m.m[0][2], &cur, &tmp);

    fp2_mul(&mut m.m[0][3], &t101, &t102);
    fp2_mul(&mut tmp, &gi[2].g10, &gi[3].g10);
    let cur = m.m[0][3];
    fp2_add(&mut m.m[0][3], &cur, &tmp);
    fp2_mul(&mut tmp, &gi[0].g10, &gi[1].g10);
    let cur = m.m[0][3];
    fp2_add(&mut m.m[0][3], &cur, &tmp);

    let m01 = m.m[0][1];
    let m00 = m.m[0][0];
    fp2_mul(&mut tmp, &gi[3].g01, &m01);
    fp2_mul(&mut m.m[1][0], &gi[3].g00, &m00);
    let cur = m.m[1][0];
    fp2_add(&mut m.m[1][0], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[3].g11, &m01);
    fp2_mul(&mut m.m[1][1], &gi[3].g10, &m00);
    let cur = m.m[1][1];
    fp2_add(&mut m.m[1][1], &cur, &tmp);

    let m03 = m.m[0][3];
    let m02 = m.m[0][2];
    fp2_mul(&mut tmp, &gi[3].g01, &m03);
    fp2_mul(&mut m.m[1][2], &gi[3].g00, &m02);
    let cur = m.m[1][2];
    fp2_add(&mut m.m[1][2], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[3].g11, &m03);
    fp2_mul(&mut m.m[1][3], &gi[3].g10, &m02);
    let cur = m.m[1][3];
    fp2_add(&mut m.m[1][3], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g01, &m02);
    fp2_mul(&mut m.m[2][0], &gi[0].g00, &m00);
    let cur = m.m[2][0];
    fp2_add(&mut m.m[2][0], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g01, &m03);
    fp2_mul(&mut m.m[2][1], &gi[0].g00, &m01);
    let cur = m.m[2][1];
    fp2_add(&mut m.m[2][1], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &m02);
    fp2_mul(&mut m.m[2][2], &gi[0].g10, &m00);
    let cur = m.m[2][2];
    fp2_add(&mut m.m[2][2], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &m03);
    fp2_mul(&mut m.m[2][3], &gi[0].g10, &m01);
    let cur = m.m[2][3];
    fp2_add(&mut m.m[2][3], &cur, &tmp);

    let m12 = m.m[1][2];
    let m13 = m.m[1][3];
    let m10 = m.m[1][0];
    let m11 = m.m[1][1];
    fp2_mul(&mut tmp, &gi[0].g01, &m12);
    fp2_mul(&mut m.m[3][0], &gi[0].g00, &m10);
    let cur = m.m[3][0];
    fp2_add(&mut m.m[3][0], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g01, &m13);
    fp2_mul(&mut m.m[3][1], &gi[0].g00, &m11);
    let cur = m.m[3][1];
    fp2_add(&mut m.m[3][1], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &m12);
    fp2_mul(&mut m.m[3][2], &gi[0].g10, &m10);
    let cur = m.m[3][2];
    fp2_add(&mut m.m[3][2], &cur, &tmp);

    fp2_mul(&mut tmp, &gi[0].g11, &m13);
    fp2_mul(&mut m.m[3][3], &gi[0].g10, &m11);
    let cur = m.m[3][3];
    fp2_add(&mut m.m[3][3], &cur, &tmp);

    1
}

fn gluing_compute(
    out: &mut ThetaGluing,
    e12: &ThetaCoupleCurve,
    xy_k1_8: &ThetaCoupleJacPoint,
    xy_k2_8: &ThetaCoupleJacPoint,
    verify: bool,
) -> i32 {
    // The reference performs `assert`-style sanity checks in NDEBUG-off
    // builds; we keep them as discardable diagnostics. They affect no
    // boundary output.
    let _ = (
        test_jac_order_twof(&xy_k1_8.p1, &e12.e1, 3),
        test_jac_order_twof(&xy_k2_8.p1, &e12.e1, 3),
        test_jac_order_twof(&xy_k1_8.p2, &e12.e2, 3),
        test_jac_order_twof(&xy_k2_8.p2, &e12.e2, 3),
    );

    out.xy_k1_8 = *xy_k1_8;
    out.domain = *e12;

    let mut xy_k1_4 = ThetaCoupleJacPoint::zero();
    let mut xy_k2_4 = ThetaCoupleJacPoint::zero();
    double_couple_jac_point(&mut xy_k1_4, xy_k1_8, e12);
    double_couple_jac_point(&mut xy_k2_4, xy_k2_8, e12);

    let mut k1_8 = ThetaCouplePoint::zero();
    let mut k2_8 = ThetaCouplePoint::zero();
    let mut k1_4 = ThetaCouplePoint::zero();
    let mut k2_4 = ThetaCouplePoint::zero();
    couple_jac_to_xz(&mut k1_8, xy_k1_8);
    couple_jac_to_xz(&mut k2_8, xy_k2_8);
    couple_jac_to_xz(&mut k1_4, &xy_k1_4);
    couple_jac_to_xz(&mut k2_4, &xy_k2_4);

    if gluing_change_of_basis(&mut out.mat, &k1_4, &k2_4, e12) == 0 {
        return 0;
    }

    let mut tt1 = ThetaPoint::zero();
    let mut tt2 = ThetaPoint::zero();
    base_change(&mut tt1, out, &k1_8);
    base_change(&mut tt2, out, &k2_8);

    let cur = tt1;
    to_squared_theta(&mut tt1, &cur);
    let cur = tt2;
    to_squared_theta(&mut tt2, &cur);

    if (fp2_is_zero(&tt1.t) & fp2_is_zero(&tt2.t)) == 0 {
        return 0;
    }
    if (fp2_is_zero(&tt1.x)
        | fp2_is_zero(&tt2.x)
        | fp2_is_zero(&tt1.y)
        | fp2_is_zero(&tt2.z)
        | fp2_is_zero(&tt1.z))
        != 0
    {
        return 0;
    }

    fp2_mul(&mut out.codomain.x, &tt1.x, &tt2.x);
    fp2_mul(&mut out.codomain.y, &tt1.y, &tt2.x);
    fp2_mul(&mut out.codomain.z, &tt1.x, &tt2.z);
    fp2_set_zero(&mut out.codomain.t);

    fp2_mul(&mut out.precomputation.x, &tt1.y, &tt2.z);
    let cd_z = out.codomain.z;
    fp2_copy(&mut out.precomputation.y, &cd_z);
    let cd_y = out.codomain.y;
    fp2_copy(&mut out.precomputation.z, &cd_y);
    fp2_set_zero(&mut out.precomputation.t);

    fp2_mul(&mut out.image_k1_8.x, &tt1.x, &out.precomputation.x);
    fp2_mul(&mut out.image_k1_8.y, &tt1.z, &out.precomputation.z);

    if verify {
        let mut t1 = fp2_zero();
        let mut t2 = fp2_zero();
        fp2_mul(&mut t1, &tt1.y, &out.precomputation.y);
        if fp2_is_equal(&out.image_k1_8.x, &t1) == 0 {
            return 0;
        }
        fp2_mul(&mut t1, &tt2.x, &out.precomputation.x);
        fp2_mul(&mut t2, &tt2.z, &out.precomputation.z);
        if fp2_is_equal(&t2, &t1) == 0 {
            return 0;
        }
    }

    let cur = out.codomain;
    hadamard(&mut out.codomain, &cur);
    1
}

fn gluing_eval_point(image: &mut ThetaPoint, p: &ThetaCoupleJacPoint, phi: &ThetaGluing) {
    let mut t1 = ThetaPoint::zero();
    let mut t2 = ThetaPoint::zero();
    let mut add_comp1 = AddComponents::zero();
    let mut add_comp2 = AddComponents::zero();

    jac_to_xz_add_components(&mut add_comp1, &p.p1, &phi.xy_k1_8.p1, &phi.domain.e1);
    jac_to_xz_add_components(&mut add_comp2, &p.p2, &phi.xy_k1_8.p2, &phi.domain.e2);

    fp2_mul(&mut t1.x, &add_comp1.u, &add_comp2.u);
    fp2_mul(&mut t2.t, &add_comp1.v, &add_comp2.v);
    let cur = t1.x;
    fp2_add(&mut t1.x, &cur, &t2.t);
    fp2_mul(&mut t1.y, &add_comp1.u, &add_comp2.w);
    fp2_mul(&mut t1.z, &add_comp1.w, &add_comp2.u);
    fp2_mul(&mut t1.t, &add_comp1.w, &add_comp2.w);
    fp2_add(&mut t2.x, &add_comp1.u, &add_comp1.v);
    fp2_add(&mut t2.y, &add_comp2.u, &add_comp2.v);
    let cur = t2.x;
    fp2_mul(&mut t2.x, &cur, &t2.y);
    let cur = t2.x;
    fp2_sub(&mut t2.x, &cur, &t1.x);
    fp2_mul(&mut t2.y, &add_comp1.v, &add_comp2.w);
    fp2_mul(&mut t2.z, &add_comp1.w, &add_comp2.v);
    fp2_set_zero(&mut t2.t);

    let cur = t1;
    apply_isomorphism_general(&mut t1, &phi.mat, &cur, true);
    let cur = t2;
    apply_isomorphism_general(&mut t2, &phi.mat, &cur, false);
    let cur = t1;
    pointwise_square(&mut t1, &cur);
    let cur = t2;
    pointwise_square(&mut t2, &cur);

    let cur = t1.x;
    fp2_sub(&mut t1.x, &cur, &t2.x);
    let cur = t1.y;
    fp2_sub(&mut t1.y, &cur, &t2.y);
    let cur = t1.z;
    fp2_sub(&mut t1.z, &cur, &t2.z);
    let cur = t1.t;
    fp2_sub(&mut t1.t, &cur, &t2.t);
    let cur = t1;
    hadamard(&mut t1, &cur);

    fp2_mul(&mut image.x, &t1.x, &phi.image_k1_8.y);
    fp2_mul(&mut image.y, &t1.y, &phi.image_k1_8.y);
    fp2_mul(&mut image.z, &t1.z, &phi.image_k1_8.x);
    fp2_mul(&mut image.t, &t1.t, &phi.image_k1_8.x);

    let cur = *image;
    hadamard(image, &cur);
}

fn gluing_eval_point_special_case(
    image: &mut ThetaPoint,
    p: &ThetaCouplePoint,
    phi: &ThetaGluing,
) -> i32 {
    let mut t = ThetaPoint::zero();
    base_change(&mut t, phi, p);

    let cur = t;
    to_squared_theta(&mut t, &cur);

    if fp2_is_zero(&t.t) == 0 {
        return 0;
    }

    fp2_mul(&mut image.x, &t.x, &phi.precomputation.x);
    fp2_mul(&mut image.y, &t.y, &phi.precomputation.y);
    fp2_mul(&mut image.z, &t.z, &phi.precomputation.z);
    fp2_set_zero(&mut image.t);

    let cur = *image;
    hadamard(image, &cur);
    1
}

fn gluing_eval_basis(
    image1: &mut ThetaPoint,
    image2: &mut ThetaPoint,
    xy_t1: &ThetaCoupleJacPoint,
    xy_t2: &ThetaCoupleJacPoint,
    phi: &ThetaGluing,
) {
    gluing_eval_point(image1, xy_t1, phi);
    gluing_eval_point(image2, xy_t2, phi);
}

fn theta_isogeny_compute(
    out: &mut ThetaIsogeny,
    a: &ThetaStructure,
    t1_8: &ThetaPoint,
    t2_8: &ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
    verify: bool,
) -> i32 {
    out.hadamard_bool_1 = hadamard_bool_1;
    out.hadamard_bool_2 = hadamard_bool_2;
    out.domain = *a;
    out.t1_8 = *t1_8;
    out.t2_8 = *t2_8;
    out.codomain.precomputation = false;

    let mut tt1 = ThetaPoint::zero();
    let mut tt2 = ThetaPoint::zero();

    if hadamard_bool_1 {
        hadamard(&mut tt1, t1_8);
        let cur = tt1;
        to_squared_theta(&mut tt1, &cur);
        hadamard(&mut tt2, t2_8);
        let cur = tt2;
        to_squared_theta(&mut tt2, &cur);
    } else {
        to_squared_theta(&mut tt1, t1_8);
        to_squared_theta(&mut tt2, t2_8);
    }

    let mut t1 = fp2_zero();
    let mut t2 = fp2_zero();

    if (fp2_is_zero(&tt2.x)
        | fp2_is_zero(&tt2.y)
        | fp2_is_zero(&tt2.z)
        | fp2_is_zero(&tt2.t)
        | fp2_is_zero(&tt1.x)
        | fp2_is_zero(&tt1.y))
        != 0
    {
        return 0;
    }

    fp2_mul(&mut t1, &tt1.x, &tt2.y);
    fp2_mul(&mut t2, &tt1.y, &tt2.x);
    fp2_mul(&mut out.codomain.null_point.x, &tt2.x, &t1);
    fp2_mul(&mut out.codomain.null_point.y, &tt2.y, &t2);
    fp2_mul(&mut out.codomain.null_point.z, &tt2.z, &t1);
    fp2_mul(&mut out.codomain.null_point.t, &tt2.t, &t2);

    let mut t3 = fp2_zero();
    fp2_mul(&mut t3, &tt2.z, &tt2.t);
    fp2_mul(&mut out.precomputation.x, &t3, &tt1.y);
    fp2_mul(&mut out.precomputation.y, &t3, &tt1.x);
    let cd_t = out.codomain.null_point.t;
    fp2_copy(&mut out.precomputation.z, &cd_t);
    let cd_z = out.codomain.null_point.z;
    fp2_copy(&mut out.precomputation.t, &cd_z);

    if verify {
        fp2_mul(&mut t1, &tt1.x, &out.precomputation.x);
        fp2_mul(&mut t2, &tt1.y, &out.precomputation.y);
        if fp2_is_equal(&t1, &t2) == 0 {
            return 0;
        }
        fp2_mul(&mut t1, &tt1.z, &out.precomputation.z);
        fp2_mul(&mut t2, &tt1.t, &out.precomputation.t);
        if fp2_is_equal(&t1, &t2) == 0 {
            return 0;
        }
        fp2_mul(&mut t1, &tt2.x, &out.precomputation.x);
        fp2_mul(&mut t2, &tt2.z, &out.precomputation.z);
        if fp2_is_equal(&t1, &t2) == 0 {
            return 0;
        }
        fp2_mul(&mut t1, &tt2.y, &out.precomputation.y);
        fp2_mul(&mut t2, &tt2.t, &out.precomputation.t);
        if fp2_is_equal(&t1, &t2) == 0 {
            return 0;
        }
    }

    if hadamard_bool_2 {
        let cur = out.codomain.null_point;
        hadamard(&mut out.codomain.null_point, &cur);
    }
    1
}

fn theta_isogeny_compute_4(
    out: &mut ThetaIsogeny,
    a: &ThetaStructure,
    t1_4: &ThetaPoint,
    t2_4: &ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
) {
    out.hadamard_bool_1 = hadamard_bool_1;
    out.hadamard_bool_2 = hadamard_bool_2;
    out.domain = *a;
    out.t1_8 = *t1_4;
    out.t2_8 = *t2_4;
    out.codomain.precomputation = false;

    let mut tt1 = ThetaPoint::zero();
    let mut tt2 = ThetaPoint::zero();

    if hadamard_bool_1 {
        hadamard(&mut tt1, t1_4);
        let cur = tt1;
        to_squared_theta(&mut tt1, &cur);
        hadamard(&mut tt2, &a.null_point);
        let cur = tt2;
        to_squared_theta(&mut tt2, &cur);
    } else {
        to_squared_theta(&mut tt1, t1_4);
        to_squared_theta(&mut tt2, &a.null_point);
    }

    let mut sqaabb = fp2_zero();
    let mut sqaacc = fp2_zero();
    fp2_mul(&mut sqaabb, &tt2.x, &tt2.y);
    fp2_mul(&mut sqaacc, &tt2.x, &tt2.z);
    fp2_sqrt(&mut sqaabb);
    fp2_sqrt(&mut sqaacc);

    fp2_mul(&mut out.codomain.null_point.y, &sqaabb, &sqaacc);
    fp2_mul(
        &mut out.precomputation.t,
        &out.codomain.null_point.y,
        &tt1.z,
    );
    let cur = out.codomain.null_point.y;
    fp2_mul(&mut out.codomain.null_point.y, &cur, &tt1.x);

    fp2_mul(&mut out.codomain.null_point.t, &tt1.z, &sqaabb);
    let cur = out.codomain.null_point.t;
    fp2_mul(&mut out.codomain.null_point.t, &cur, &tt2.x);

    fp2_mul(&mut out.codomain.null_point.x, &tt1.x, &tt2.x);
    fp2_mul(
        &mut out.codomain.null_point.z,
        &out.codomain.null_point.x,
        &tt2.z,
    );
    let cur = out.codomain.null_point.x;
    fp2_mul(&mut out.codomain.null_point.x, &cur, &sqaacc);

    fp2_mul(&mut out.precomputation.x, &tt1.x, &tt2.t);
    fp2_mul(&mut out.precomputation.z, &out.precomputation.x, &tt2.y);
    let cur = out.precomputation.x;
    fp2_mul(&mut out.precomputation.x, &cur, &tt2.z);
    fp2_mul(&mut out.precomputation.y, &out.precomputation.x, &sqaabb);
    let cur = out.precomputation.x;
    fp2_mul(&mut out.precomputation.x, &cur, &tt2.y);
    let cur = out.precomputation.z;
    fp2_mul(&mut out.precomputation.z, &cur, &sqaacc);
    let cur = out.precomputation.t;
    fp2_mul(&mut out.precomputation.t, &cur, &tt2.y);

    if hadamard_bool_2 {
        let cur = out.codomain.null_point;
        hadamard(&mut out.codomain.null_point, &cur);
    }
}

fn theta_isogeny_compute_2(
    out: &mut ThetaIsogeny,
    a: &ThetaStructure,
    t1_2: &ThetaPoint,
    t2_2: &ThetaPoint,
    hadamard_bool_1: bool,
    hadamard_bool_2: bool,
) {
    out.hadamard_bool_1 = hadamard_bool_1;
    out.hadamard_bool_2 = hadamard_bool_2;
    out.domain = *a;
    out.t1_8 = *t1_2;
    out.t2_8 = *t2_2;
    out.codomain.precomputation = false;

    let mut tt2 = ThetaPoint::zero();

    if hadamard_bool_1 {
        hadamard(&mut tt2, &a.null_point);
        let cur = tt2;
        to_squared_theta(&mut tt2, &cur);
    } else {
        to_squared_theta(&mut tt2, &a.null_point);
    }

    fp2_copy(&mut out.codomain.null_point.x, &tt2.x);
    fp2_mul(&mut out.codomain.null_point.y, &tt2.x, &tt2.y);
    fp2_mul(&mut out.codomain.null_point.z, &tt2.x, &tt2.z);
    fp2_mul(&mut out.codomain.null_point.t, &tt2.x, &tt2.t);
    fp2_sqrt(&mut out.codomain.null_point.y);
    fp2_sqrt(&mut out.codomain.null_point.z);
    fp2_sqrt(&mut out.codomain.null_point.t);

    fp2_mul(&mut out.precomputation.x, &tt2.z, &tt2.t);
    fp2_mul(
        &mut out.precomputation.y,
        &out.precomputation.x,
        &out.codomain.null_point.y,
    );
    let cur = out.precomputation.x;
    fp2_mul(&mut out.precomputation.x, &cur, &tt2.y);
    fp2_mul(
        &mut out.precomputation.z,
        &tt2.t,
        &out.codomain.null_point.z,
    );
    let cur = out.precomputation.z;
    fp2_mul(&mut out.precomputation.z, &cur, &tt2.y);
    fp2_mul(
        &mut out.precomputation.t,
        &tt2.z,
        &out.codomain.null_point.t,
    );
    let cur = out.precomputation.t;
    fp2_mul(&mut out.precomputation.t, &cur, &tt2.y);

    if hadamard_bool_2 {
        let cur = out.codomain.null_point;
        hadamard(&mut out.codomain.null_point, &cur);
    }
}

fn theta_isogeny_eval(out: &mut ThetaPoint, phi: &ThetaIsogeny, p: &ThetaPoint) {
    if phi.hadamard_bool_1 {
        hadamard(out, p);
        let cur = *out;
        to_squared_theta(out, &cur);
    } else {
        to_squared_theta(out, p);
    }
    let cur = *out;
    fp2_mul(&mut out.x, &cur.x, &phi.precomputation.x);
    fp2_mul(&mut out.y, &cur.y, &phi.precomputation.y);
    fp2_mul(&mut out.z, &cur.z, &phi.precomputation.z);
    fp2_mul(&mut out.t, &cur.t, &phi.precomputation.t);

    if phi.hadamard_bool_2 {
        let cur = *out;
        hadamard(out, &cur);
    }
}

fn splitting_compute(
    out: &mut ThetaSplitting,
    a: &ThetaStructure,
    zero_index: i32,
    randomize: bool,
) -> bool {
    let mut count: u32 = 0;
    let mut u_cst = fp2_zero();
    let mut t1 = fp2_zero();
    let mut t2 = fp2_zero();

    // Zero the matrix.
    out.mat = BasisChangeMatrix::zero();

    for i in 0..10usize {
        fp2_set_zero(&mut u_cst);
        for t in 0..4i32 {
            choose_index_theta_point(&mut t2, t, &a.null_point);
            choose_index_theta_point(&mut t1, t ^ EVEN_INDEX[i][1] as i32, &a.null_point);

            let cur = t1;
            fp2_mul(&mut t1, &cur, &t2);
            // CHI_EVAL[a][t] in {-1, +1}; cast via i32 -> uint32 then `>> 1`
            // yields 0 (positive) or 0xFFFFFFFF (negative two's-complement).
            let chi = CHI_EVAL[EVEN_INDEX[i][0] as usize][t as usize];
            let ctl: u32 = (chi as u32) >> 1;
            debug_assert!(ctl == 0 || ctl == 0xFFFF_FFFF);

            fp2_neg(&mut t2, &t1);
            let cur = t1;
            fp2_select(&mut t1, &cur, &t2, ctl);

            let cur = u_cst;
            fp2_add(&mut u_cst, &cur, &t1);
        }

        let ctl = fp2_is_zero(&u_cst);
        count = count.wrapping_sub(ctl);
        let cur = out.mat;
        select_base_change_matrix(&mut out.mat, &cur, &SPLITTING_TRANSFORMS[i], ctl);
        if zero_index != -1 && i as i32 == zero_index && ctl == 0 {
            return false;
        }
    }

    if randomize {
        // Deferred in this unit: see crate-level docs for the
        // `theta_chain_compute_and_eval_randomized` note. The two
        // referenced helpers (`set_base_change_matrix_from_precomp`,
        // `base_change_matrix_multiplication`) are kept so the
        // non-random path can call them locally if needed and so the
        // randomised pathway is a small extension.
        let _ = NORMALIZATION_TRANSFORMS;
        let _ = (
            set_base_change_matrix_from_precomp as fn(_, _),
            base_change_matrix_multiplication as fn(_, _, _),
        );
        unimplemented!(
            "theta_chain_compute_and_eval_randomized: requires caller-supplied RNG (deferred)"
        );
    }

    apply_isomorphism(&mut out.b.null_point, &out.mat, &a.null_point);

    count == 1
}

fn theta_product_structure_to_elliptic_product(
    e12: &mut ThetaCoupleCurve,
    a: &ThetaStructure,
) -> i32 {
    let mut xx = fp2_zero();
    let mut yy = fp2_zero();

    if is_product_theta_point(&a.null_point) == 0 {
        return 0;
    }

    ec_curve_init(&mut e12.e1);
    ec_curve_init(&mut e12.e2);

    if (fp2_is_zero(&a.null_point.x) | fp2_is_zero(&a.null_point.y) | fp2_is_zero(&a.null_point.z))
        != 0
    {
        return 0;
    }

    fp2_sqr(&mut xx, &a.null_point.x);
    fp2_sqr(&mut yy, &a.null_point.y);
    let cur = xx;
    fp2_sqr(&mut xx, &cur);
    let cur = yy;
    fp2_sqr(&mut yy, &cur);

    fp2_add(&mut e12.e2.A, &xx, &yy);
    fp2_sub(&mut e12.e2.C, &xx, &yy);
    let cur = e12.e2.A;
    fp2_add(&mut e12.e2.A, &cur, &cur);
    let cur = e12.e2.A;
    fp2_neg(&mut e12.e2.A, &cur);

    fp2_sqr(&mut xx, &a.null_point.x);
    fp2_sqr(&mut yy, &a.null_point.z);
    let cur = xx;
    fp2_sqr(&mut xx, &cur);
    let cur = yy;
    fp2_sqr(&mut yy, &cur);

    fp2_add(&mut e12.e1.A, &xx, &yy);
    fp2_sub(&mut e12.e1.C, &xx, &yy);
    let cur = e12.e1.A;
    fp2_add(&mut e12.e1.A, &cur, &cur);
    let cur = e12.e1.A;
    fp2_neg(&mut e12.e1.A, &cur);

    if (fp2_is_zero(&e12.e1.C) | fp2_is_zero(&e12.e2.C)) != 0 {
        return 0;
    }

    1
}

fn theta_point_to_montgomery_point(
    p12: &mut ThetaCouplePoint,
    p: &ThetaPoint,
    a: &ThetaStructure,
) -> i32 {
    if is_product_theta_point(p) == 0 {
        return 0;
    }

    // The C source uses `const fp2_t *` aliasing through P; we mirror by
    // copying the selected coordinates and reading from there.
    let (mut x_sel, mut z_sel) = (p.x, p.y);
    if (fp2_is_zero(&x_sel) & fp2_is_zero(&z_sel)) != 0 {
        x_sel = p.z;
        z_sel = p.t;
    }
    if (fp2_is_zero(&x_sel) & fp2_is_zero(&z_sel)) != 0 {
        return 0;
    }
    let mut temp = fp2_zero();
    fp2_mul(&mut p12.p2.x, &a.null_point.y, &x_sel);
    fp2_mul(&mut temp, &a.null_point.x, &z_sel);
    fp2_sub(&mut p12.p2.z, &temp, &p12.p2.x);
    let cur = p12.p2.x;
    fp2_add(&mut p12.p2.x, &cur, &temp);

    let (mut x_sel, mut z_sel) = (p.x, p.z);
    if (fp2_is_zero(&x_sel) & fp2_is_zero(&z_sel)) != 0 {
        x_sel = p.y;
        z_sel = p.t;
    }
    fp2_mul(&mut p12.p1.x, &a.null_point.z, &x_sel);
    fp2_mul(&mut temp, &a.null_point.x, &z_sel);
    fp2_sub(&mut p12.p1.z, &temp, &p12.p1.x);
    let cur = p12.p1.x;
    fp2_add(&mut p12.p1.x, &cur, &temp);
    1
}

#[allow(clippy::too_many_arguments)]
fn theta_chain_compute_impl(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
    verify: bool,
    randomize: bool,
) -> i32 {
    let num_p = p12.len();
    let mut theta = ThetaStructure::zero();

    let mut xy_t1 = ThetaCoupleJacPoint::zero();
    let mut xy_t2 = ThetaCoupleJacPoint::zero();

    let mut bas1 = EcBasis {
        P: ker.t1.p1,
        Q: ker.t2.p1,
        PmQ: ker.t1m2.p1,
    };
    let mut bas2 = EcBasis {
        P: ker.t1.p2,
        Q: ker.t2.p2,
        PmQ: ker.t1m2.p2,
    };
    if lift_basis(&mut xy_t1.p1, &mut xy_t2.p1, &mut bas1, &mut e12.e1) == 0 {
        return 0;
    }
    if lift_basis(&mut xy_t1.p2, &mut xy_t2.p2, &mut bas2, &mut e12.e2) == 0 {
        return 0;
    }

    let extra: u32 = HD_EXTRA_TORSION * (extra_torsion as u32);
    debug_assert!(extra == 0 || extra == 2);

    let mut pts: Vec<ThetaPoint> = vec![ThetaPoint::zero(); num_p];

    // `space = ceil(log2(n)) + 1`; mirrors the C `space` counter.
    let mut space: usize = 1;
    {
        let mut i: u32 = 1;
        while i < n {
            i *= 2;
            space += 1;
        }
    }

    let mut todo: Vec<u16> = vec![0u16; space];
    todo[0] = (n - 2 + extra) as u16;

    let mut current: i32 = 0;

    let mut jac_q1: Vec<ThetaCoupleJacPoint> = vec![ThetaCoupleJacPoint::zero(); space];
    let mut jac_q2: Vec<ThetaCoupleJacPoint> = vec![ThetaCoupleJacPoint::zero(); space];
    jac_q1[0] = xy_t1;
    jac_q2[0] = xy_t2;

    while todo[current as usize] != 1 {
        debug_assert!(todo[current as usize] >= 2);
        current += 1;
        debug_assert!((current as usize) < space);
        let prev = todo[(current - 1) as usize];
        let num_dbls: u32 = if prev >= 16 {
            (prev / 2) as u32
        } else {
            (prev - 1) as u32
        };
        debug_assert!(num_dbls != 0 && num_dbls < prev as u32);
        let prev_q1 = jac_q1[(current - 1) as usize];
        let prev_q2 = jac_q2[(current - 1) as usize];
        double_couple_jac_point_iter(&mut jac_q1[current as usize], num_dbls, &prev_q1, e12);
        double_couple_jac_point_iter(&mut jac_q2[current as usize], num_dbls, &prev_q2, e12);
        todo[current as usize] = (prev as u32 - num_dbls) as u16;
    }

    let mut theta_q1: Vec<ThetaPoint> = vec![ThetaPoint::zero(); space];
    let mut theta_q2: Vec<ThetaPoint> = vec![ThetaPoint::zero(); space];

    let mut first_step = ThetaGluing::zero();
    {
        debug_assert_eq!(todo[current as usize], 1);
        if gluing_compute(
            &mut first_step,
            e12,
            &jac_q1[current as usize],
            &jac_q2[current as usize],
            verify,
        ) == 0
        {
            return 0;
        }

        for j in 0..num_p {
            debug_assert!(ec_is_zero(&p12[j].p1) != 0 || ec_is_zero(&p12[j].p2) != 0);
            if gluing_eval_point_special_case(&mut pts[j], &p12[j], &first_step) == 0 {
                return 0;
            }
        }

        for j in 0..(current as usize) {
            let q1 = jac_q1[j];
            let q2 = jac_q2[j];
            gluing_eval_basis(&mut theta_q1[j], &mut theta_q2[j], &q1, &q2, &first_step);
            todo[j] -= 1;
        }

        current -= 1;
    }

    theta.null_point = first_step.codomain;
    theta.precomputation = false;
    theta_precomputation(&mut theta);

    let mut step = ThetaIsogeny::zero();

    let mut i: u32 = 1;
    while current >= 0 && todo[current as usize] != 0 {
        debug_assert!((current as usize) < space);
        while todo[current as usize] != 1 {
            debug_assert!(todo[current as usize] >= 2);
            current += 1;
            debug_assert!((current as usize) < space);
            let prev = todo[(current - 1) as usize];
            let num_dbls: i32 = (prev / 2) as i32;
            debug_assert!(num_dbls != 0 && (num_dbls as u16) < prev);
            let prev_q1 = theta_q1[(current - 1) as usize];
            let prev_q2 = theta_q2[(current - 1) as usize];
            double_iter(
                &mut theta_q1[current as usize],
                &mut theta,
                &prev_q1,
                num_dbls,
            );
            double_iter(
                &mut theta_q2[current as usize],
                &mut theta,
                &prev_q2,
                num_dbls,
            );
            todo[current as usize] = (prev as i32 - num_dbls) as u16;
        }

        let ret = if i == n - 2 {
            theta_isogeny_compute(
                &mut step,
                &theta,
                &theta_q1[current as usize],
                &theta_q2[current as usize],
                false,
                false,
                verify,
            )
        } else if i == n - 1 {
            theta_isogeny_compute(
                &mut step,
                &theta,
                &theta_q1[current as usize],
                &theta_q2[current as usize],
                true,
                false,
                false,
            )
        } else {
            theta_isogeny_compute(
                &mut step,
                &theta,
                &theta_q1[current as usize],
                &theta_q2[current as usize],
                false,
                true,
                verify,
            )
        };
        if ret == 0 {
            return 0;
        }

        for j in 0..num_p {
            let cur = pts[j];
            theta_isogeny_eval(&mut pts[j], &step, &cur);
        }

        theta = step.codomain;

        debug_assert_eq!(todo[current as usize], 1);
        for j in 0..(current as usize) {
            let cur1 = theta_q1[j];
            let cur2 = theta_q2[j];
            theta_isogeny_eval(&mut theta_q1[j], &step, &cur1);
            theta_isogeny_eval(&mut theta_q2[j], &step, &cur2);
            debug_assert!(todo[j] != 0);
            todo[j] -= 1;
        }

        current -= 1;
        i += 1;
    }

    debug_assert_eq!(current, -1);

    if !extra_torsion {
        if n >= 3 {
            let cur1 = theta_q1[0];
            let cur2 = theta_q2[0];
            theta_isogeny_eval(&mut theta_q1[0], &step, &cur1);
            theta_isogeny_eval(&mut theta_q2[0], &step, &cur2);
        }

        theta_isogeny_compute_4(&mut step, &theta, &theta_q1[0], &theta_q2[0], false, false);
        for j in 0..num_p {
            let cur = pts[j];
            theta_isogeny_eval(&mut pts[j], &step, &cur);
        }
        theta = step.codomain;
        let cur1 = theta_q1[0];
        let cur2 = theta_q2[0];
        theta_isogeny_eval(&mut theta_q1[0], &step, &cur1);
        theta_isogeny_eval(&mut theta_q2[0], &step, &cur2);

        theta_isogeny_compute_2(&mut step, &theta, &theta_q1[0], &theta_q2[0], true, false);
        for j in 0..num_p {
            let cur = pts[j];
            theta_isogeny_eval(&mut pts[j], &step, &cur);
        }
        theta = step.codomain;
    }

    let mut last_step = ThetaSplitting::zero();
    let is_split = splitting_compute(
        &mut last_step,
        &theta,
        if extra_torsion { 8 } else { -1 },
        randomize,
    );
    if !is_split {
        return 0;
    }

    if theta_product_structure_to_elliptic_product(e34, &last_step.b) == 0 {
        return 0;
    }

    for j in 0..num_p {
        let cur = pts[j];
        apply_isomorphism(&mut pts[j], &last_step.mat, &cur);
        if theta_point_to_montgomery_point(&mut p12[j], &pts[j], &last_step.b) == 0 {
            return 0;
        }
    }

    1
}

/// Mirrors `theta_chain_compute_and_eval`: compute a (2,2)-isogeny
/// chain of length `n` between elliptic products in the theta model
/// and evaluate at every point of `p12`. Returns 1 on success, 0 on
/// failure.
pub fn theta_chain_compute_and_eval(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
) -> i32 {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, e34, p12, false, false)
}

/// Mirrors `theta_chain_compute_and_eval_verify`: like
/// [`theta_chain_compute_and_eval`] but with extra isotropy checks on
/// the kernel during the gluing and intermediate steps.
pub fn theta_chain_compute_and_eval_verify(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
) -> i32 {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, e34, p12, true, false)
}

/// Mirrors `theta_chain_compute_and_eval_randomized`: as
/// [`theta_chain_compute_and_eval`] but draws a random normalisation
/// matrix during the splitting step. **Deferred** in this unit; see
/// the crate-level docs.
pub fn theta_chain_compute_and_eval_randomized(
    n: u32,
    e12: &mut ThetaCoupleCurve,
    ker: &ThetaKernelCouplePoints,
    extra_torsion: bool,
    e34: &mut ThetaCoupleCurve,
    p12: &mut [ThetaCouplePoint],
) -> i32 {
    theta_chain_compute_impl(n, e12, ker, extra_torsion, e34, p12, false, true)
}
