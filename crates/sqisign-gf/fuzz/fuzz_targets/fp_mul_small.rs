//! Invariant fuzz target for `sqisign_gf::fp_mul_small`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_mul_small` on a fuzzing host.
//!
//! `fp_mul_small` is the reference's `modmli(*a, (int)val, *x)`: narrow
//! `val` to `int32`, build the Montgomery representative of the
//! narrowed-and-sign-extended integer via `modint` into a five-limb
//! scratch, then `modmul(a, scratch, x)` into the destination. Three
//! sound raw-limb invariants survive the bit-exact form into property
//! space and are exercised by the harness:
//!
//! 1. `val == 0` produces the canonical all-zero limb vector for
//!    arbitrary `a`. The Montgomery image of zero is zero, and every
//!    column sum in `modmul(a, [0; 5])` accumulates only zero partial
//!    products, so the final per-limb writes are all zero.
//! 2. `val == 1` produces the same redundant representative as
//!    `fp_mul(a, &MONTGOMERY_ONE)`, the literal expansion of the
//!    `modmli` two-line body. Since `modint(1)` writes the Montgomery
//!    ONE bit pattern (pinned by the `nres(positional 1) ==
//!    MONTGOMERY_ONE` unit test in lib.rs), this is the structural
//!    identity the cross-oracle for property (3) generalises.
//! 3. Cross-oracle: `fp_mul_small(out, a, val) == fp_mul(a,
//!    &fp_set_small_temp(val))` bit-for-bit, for arbitrary `a` and
//!    `val`. The boundary is literally that two-step chain (build the
//!    Montgomery image of the integer scalar, multiply), so the
//!    equivalence is a structural tautology of the port.
//! 4. Structural carry-propagation invariant: limbs 0..=3 of the
//!    output are below `2^51` (limb 4 is left fully unmasked by the
//!    underlying `modmul`'s final truncation, no `& mask`, and is not
//!    constrained). Inherited from `fp_mul`'s structural invariant.
//!
//! The fuzz input is split into a five-limb `a` plus a single u32
//! `val`, both consumed verbatim. Linking the reference's
//! `fp_p5248_64.c` for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_mul, fp_mul_small, fp_set_small, Fp, NWORDS_FIELD};

fn fp(bytes: &[u8]) -> Fp {
    let mut n = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks(8).take(NWORDS_FIELD).enumerate() {
        let mut w = [0u8; 8];
        w[..chunk.len()].copy_from_slice(chunk);
        n[i] = u64::from_le_bytes(w);
    }
    n
}

fuzz_target!(|data: &[u8]| {
    // 5 * 8 (a) + 4 (val as u32) = 44 bytes minimum.
    if data.len() < NWORDS_FIELD * 8 + 4 {
        return;
    }
    let a = fp(&data[..NWORDS_FIELD * 8]);
    let val_bytes: [u8; 4] = data[NWORDS_FIELD * 8..NWORDS_FIELD * 8 + 4]
        .try_into()
        .unwrap();
    let val = u32::from_le_bytes(val_bytes);

    // (1) val == 0 -> canonical all-zero limb vector for arbitrary a.
    let mut zero_out: Fp = [0u64; NWORDS_FIELD];
    fp_mul_small(&mut zero_out, &a, 0);
    assert_eq!(
        zero_out,
        [0u64; NWORDS_FIELD],
        "fp_mul_small(_, a, 0) must produce the canonical all-zero limb vector"
    );

    // (2) val == 1 -> same redundant representative as fp_mul(a, &MONTGOMERY_ONE).
    let mut one_direct: Fp = [0u64; NWORDS_FIELD];
    fp_mul_small(&mut one_direct, &a, 1);
    let mut one_montgomery: Fp = [0u64; NWORDS_FIELD];
    fp_set_small(&mut one_montgomery, 1);
    let mut one_via_chain: Fp = [0u64; NWORDS_FIELD];
    fp_mul(&mut one_via_chain, &a, &one_montgomery);
    assert_eq!(
        one_direct, one_via_chain,
        "fp_mul_small(_, a, 1) must equal fp_mul(a, MONTGOMERY_ONE) bit-exact"
    );

    // (3) Cross-oracle: fp_mul_small(a, val) == fp_mul(a, fp_set_small(val)).
    // Sign-extend val into the u64 fp_set_small accepts so the
    // narrowed-int32 image is reproduced exactly.
    let val_sx = val as i32 as u64;
    let mut scratch: Fp = [0u64; NWORDS_FIELD];
    fp_set_small(&mut scratch, val_sx);
    let mut via_chain: Fp = [0u64; NWORDS_FIELD];
    fp_mul(&mut via_chain, &a, &scratch);
    let mut direct: Fp = [0u64; NWORDS_FIELD];
    fp_mul_small(&mut direct, &a, val);
    assert_eq!(
        direct, via_chain,
        "fp_mul_small(_, a, val) must equal fp_mul(a, fp_set_small(val as i32 as u64)) bit-exact"
    );

    // (4) Structural carry-propagation invariant on limbs 0..=3.
    for (k, &limb) in direct.iter().take(4).enumerate() {
        assert!(
            limb < (1u64 << 51),
            "limb {k} of fp_mul_small output not reduced below 2^51"
        );
    }
});
