//! Property tests for `fp_set_zero`.
//!
//! `fp_set_zero` is the reference's `modzer`: a plain five-limb
//! zero-fill, `x[i] = 0` for `i` in `0..5`. No `prop`, no `2p`
//! correction, no reduction. Unlike `fp_add`/`fp_sub`/`fp_neg` it is not
//! constrained to the redundant-mod-`p` notion of equality: the output
//! is the **bit-exact canonical all-zero representative** on every
//! input, canonical or not.
//!
//! Two sound raw-limb properties:
//!
//! 1. **All-zero.** For arbitrary five-limb destination pre-fill,
//!    `fp_set_zero(out)` leaves `out == [0; 5]`, bit for bit. Verified
//!    against the full 1012-vector C-derived battery
//!    (`fp_set_zero_vectors.rs` pins this as a `count == total`
//!    assertion).
//! 2. **Idempotent.** Calling `fp_set_zero` twice on the same buffer
//!    leaves it bit-exact equal to a single call: zero overwritten with
//!    zero is zero. A second call is therefore observationally
//!    indistinguishable from a single call, which guards against any
//!    port that quietly accumulated state across invocations.
//!
//! Pre-fills are drawn from arbitrary `u64` limb tuples (a setter
//! accepts non-canonical destinations by construction), and the
//! all-zero pre-fill is itself a sound input that exercises the
//! coincidence between input and output without being privileged.

use proptest::prelude::*;
use sqisign_gf::{fp_set_zero, Fp, NWORDS_FIELD};

const ZERO: Fp = [0u64; NWORDS_FIELD];

proptest! {
    // (1) All-zero: out == [0; 5], bit-exact, for arbitrary five-limb
    // destination pre-fill.
    #[test]
    fn produces_all_zero(
        prefill in proptest::array::uniform5(any::<u64>()),
    ) {
        let mut out: Fp = prefill;
        fp_set_zero(&mut out);
        prop_assert_eq!(out, ZERO);
    }

    // (2) Idempotent: a second call leaves the buffer bit-exact
    // unchanged from the first, regardless of the original pre-fill.
    #[test]
    fn idempotent(
        prefill in proptest::array::uniform5(any::<u64>()),
    ) {
        let mut out: Fp = prefill;
        fp_set_zero(&mut out);
        let after_first = out;
        fp_set_zero(&mut out);
        prop_assert_eq!(after_first, ZERO);
        prop_assert_eq!(out, ZERO);
        prop_assert_eq!(out, after_first);
    }
}
