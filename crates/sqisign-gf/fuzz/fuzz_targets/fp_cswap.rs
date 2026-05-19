//! Invariant fuzz target for `sqisign_gf::fp_cswap`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_cswap` on a fuzzing host.
//!
//! `fp_cswap` is a branchless constant-time conditional swap that
//! consults only the LSB of `ctl`. The asserted invariants are the
//! ones that transfer from the bit-exact form: at `ctl & 1 == 0`,
//! both operands are unchanged limb for limb; at `ctl & 1 == 1`, the
//! two operands are exchanged limb for limb; and a double cswap with
//! the same `ctl` returns both operands to their pre-call values.
//! The fuzz input drives both fp operands and the full 32-bit `ctl`
//! (every bit pattern, not just the canonical endpoints), so any
//! port that consumed more than the LSB would diverge on the first
//! property. Linking the reference's `fp_p5248_64.c` (lvl1) for
//! byte-equality against C is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_cswap, Fp, NWORDS_FIELD};

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
    // Two fp operands plus a 4-byte ctl.
    let limbs = NWORDS_FIELD * 8;
    if data.len() < 2 * limbs + 4 {
        return;
    }
    let a = fp(&data[..limbs]);
    let b = fp(&data[limbs..2 * limbs]);
    let ctl = u32::from_le_bytes(data[2 * limbs..2 * limbs + 4].try_into().unwrap());

    // Single call: LSB partitions into no-op and swap.
    let mut a_mut: Fp = a;
    let mut b_mut: Fp = b;
    fp_cswap(&mut a_mut, &mut b_mut, ctl);
    if (ctl & 1) == 0 {
        assert_eq!(a_mut, a, "ctl & 1 == 0 must leave a unchanged");
        assert_eq!(b_mut, b, "ctl & 1 == 0 must leave b unchanged");
        for i in 0..NWORDS_FIELD {
            assert_eq!(a_mut[i], a[i]);
            assert_eq!(b_mut[i], b[i]);
        }
    } else {
        assert_eq!(a_mut, b, "ctl & 1 == 1 must place b's limbs into a");
        assert_eq!(b_mut, a, "ctl & 1 == 1 must place a's limbs into b");
        for i in 0..NWORDS_FIELD {
            assert_eq!(a_mut[i], b[i]);
            assert_eq!(b_mut[i], a[i]);
        }
    }

    // Involution: a second cswap with the same ctl returns both to
    // their pre-call values, regardless of whether the first call
    // swapped or was a no-op.
    fp_cswap(&mut a_mut, &mut b_mut, ctl);
    assert_eq!(a_mut, a, "double cswap with same ctl must restore a");
    assert_eq!(b_mut, b, "double cswap with same ctl must restore b");

    // LSB-only contract: clearing all but the LSB must produce the
    // identical post-call state, on the original inputs.
    let mut a_lsb: Fp = a;
    let mut b_lsb: Fp = b;
    fp_cswap(&mut a_lsb, &mut b_lsb, ctl & 1);
    let mut a_full: Fp = a;
    let mut b_full: Fp = b;
    fp_cswap(&mut a_full, &mut b_full, ctl);
    assert_eq!(
        a_lsb, a_full,
        "fp_cswap must depend on only the LSB of ctl (a)"
    );
    assert_eq!(
        b_lsb, b_full,
        "fp_cswap must depend on only the LSB of ctl (b)"
    );
});
