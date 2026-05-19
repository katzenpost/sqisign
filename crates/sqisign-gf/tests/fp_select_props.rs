//! Property tests for `fp_select`: branchless constant-time conditional
//! select with a documented `ctl in {0x00000000, 0xFFFFFFFF}` contract.
//!
//! These mirror the properties that pin `mp::select_ct` for the u64
//! mask case, narrowed to the two declared endpoints `fp_select`
//! accepts. `ctl` values outside those two endpoints are undefined per
//! the reference and are not exercised.
//!
//! Three properties:
//!
//! 1. **`ctl == 0` selects `a0` bit for bit.** For arbitrary five-limb
//!    `a0`, `a1`, and destination pre-fill, `fp_select(d, a0, a1, 0)`
//!    leaves `d == a0`, limb for limb.
//! 2. **`ctl == 0xFFFFFFFF` selects `a1` bit for bit.** Symmetrically,
//!    for arbitrary five-limb inputs and pre-fill,
//!    `fp_select(d, a0, a1, 0xFFFFFFFF)` leaves `d == a1`, limb for
//!    limb.
//! 3. **Per-bit blend at the two endpoints.** At each declared endpoint
//!    the output equals the per-bit blend `(a0 & !cw) | (a1 & cw)`,
//!    where `cw` is `ctl` sign-extended to `u64` (i.e. `0` and
//!    `u64::MAX` at the two endpoints). This is the exact form the
//!    `mp::select_ct` property suite pins for the u64 mask case; here
//!    it pins that `fp_select`'s output at both endpoints is the same
//!    branchless blend as `select_ct` (just constrained to the two
//!    declared `cw` values), so the bit-for-bit oracle correspondence
//!    with the reference's `d[i] = a0[i] ^ (cw & (a0[i] ^ a1[i]))` is
//!    indistinguishable from the equivalent expression.

use proptest::prelude::*;
use sqisign_gf::{fp_select, Fp, NWORDS_FIELD};

fn uniform5() -> impl Strategy<Value = Fp> {
    (
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
        any::<u64>(),
    )
        .prop_map(|(a, b, c, d, e)| [a, b, c, d, e])
}

fn sign_extend_ctl(ctl: u32) -> u64 {
    // Reproduces the C cast chain `digit_t cw = (int32_t)ctl;` exactly:
    // bit-preserving u32 -> i32, then sign-extending widening i32 -> i64,
    // then bit-cast to u64. See the lib.rs note on `fp_select`.
    (ctl as i32) as u64
}

proptest! {
    // (1) ctl == 0 selects a0 bit for bit.
    #[test]
    fn ctl_zero_selects_a0(a0 in uniform5(), a1 in uniform5(), prefill in uniform5()) {
        let mut d: Fp = prefill;
        fp_select(&mut d, &a0, &a1, 0x00000000);
        prop_assert_eq!(d, a0);
        for i in 0..NWORDS_FIELD {
            prop_assert_eq!(d[i], a0[i]);
        }
    }

    // (2) ctl == 0xFFFFFFFF selects a1 bit for bit.
    #[test]
    fn ctl_ones_selects_a1(a0 in uniform5(), a1 in uniform5(), prefill in uniform5()) {
        let mut d: Fp = prefill;
        fp_select(&mut d, &a0, &a1, 0xFFFFFFFF);
        prop_assert_eq!(d, a1);
        for i in 0..NWORDS_FIELD {
            prop_assert_eq!(d[i], a1[i]);
        }
    }

    // (3) Per-bit blend at both endpoints. cw == 0 at ctl == 0,
    // cw == u64::MAX at ctl == 0xFFFFFFFF (sign extension).
    #[test]
    fn endpoint_outputs_are_bit_blend(
        a0 in uniform5(),
        a1 in uniform5(),
        prefill_zero in uniform5(),
        prefill_ones in uniform5(),
    ) {
        for ctl in [0x00000000u32, 0xFFFFFFFFu32] {
            let cw = sign_extend_ctl(ctl);
            let prefill = if ctl == 0 { prefill_zero } else { prefill_ones };
            let mut d: Fp = prefill;
            fp_select(&mut d, &a0, &a1, ctl);
            for i in 0..NWORDS_FIELD {
                prop_assert_eq!(d[i], (a0[i] & !cw) | (a1[i] & cw));
            }
        }
        // And the sign-extension itself, the load-bearing subtlety, is
        // exactly the two endpoint values.
        prop_assert_eq!(sign_extend_ctl(0x00000000), 0u64);
        prop_assert_eq!(sign_extend_ctl(0xFFFFFFFF), u64::MAX);
    }
}
