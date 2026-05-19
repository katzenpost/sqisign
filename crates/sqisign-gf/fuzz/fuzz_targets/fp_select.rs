//! Invariant fuzz target for `sqisign_gf::fp_select`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_select` on a fuzzing host.
//!
//! `fp_select` is a branchless constant-time conditional select with a
//! documented `ctl in {0x00000000, 0xFFFFFFFF}` contract. The asserted
//! invariants are the ones that transfer from the bit-exact form,
//! narrowed to the two declared endpoints: at `ctl == 0` the
//! destination equals `a0` limb for limb; at `ctl == 0xFFFFFFFF` it
//! equals `a1` limb for limb; and at both endpoints the output equals
//! the per-bit blend `(a0 & !cw) | (a1 & cw)`, where `cw` is `ctl`
//! sign-extended to `u64`. The fuzz input drives the destination
//! pre-fill alongside `a0` and `a1` so a no-op or partial-write port
//! would leave a visible residue. `ctl` values outside the two
//! declared endpoints are undefined per the reference's contract and
//! are deliberately not exercised. Linking the reference's
//! `fp.c` (lvlx) for byte-equality against C is the documented next
//! increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_select, Fp, NWORDS_FIELD};

fn fp(bytes: &[u8]) -> Fp {
    let mut n = [0u64; NWORDS_FIELD];
    for (i, chunk) in bytes.chunks(8).take(NWORDS_FIELD).enumerate() {
        let mut w = [0u8; 8];
        w[..chunk.len()].copy_from_slice(chunk);
        n[i] = u64::from_le_bytes(w);
    }
    n
}

fn sign_extend_ctl(ctl: u32) -> u64 {
    (ctl as i32) as u64
}

fuzz_target!(|data: &[u8]| {
    // Three fp inputs: a0, a1, and the destination pre-fill.
    let limbs = NWORDS_FIELD * 8;
    if data.len() < 3 * limbs {
        return;
    }
    let a0 = fp(&data[..limbs]);
    let a1 = fp(&data[limbs..2 * limbs]);
    let prefill = fp(&data[2 * limbs..3 * limbs]);

    for ctl in [0x00000000u32, 0xFFFFFFFFu32] {
        let cw = sign_extend_ctl(ctl);
        let mut d: Fp = prefill;
        fp_select(&mut d, &a0, &a1, ctl);
        if ctl == 0 {
            assert_eq!(d, a0, "ctl == 0 must select a0 bit for bit");
        } else {
            assert_eq!(d, a1, "ctl == 0xFFFFFFFF must select a1 bit for bit");
        }
        for i in 0..NWORDS_FIELD {
            assert_eq!(
                d[i],
                (a0[i] & !cw) | (a1[i] & cw),
                "fp_select is not the per-bit blend at ctl=0x{ctl:08x}"
            );
        }
    }

    // Sign-extension subtlety: the two endpoint values must be exactly 0
    // and u64::MAX, with no intermediate truncation. Pinned outside the
    // call to keep the property visible.
    assert_eq!(sign_extend_ctl(0x00000000), 0u64);
    assert_eq!(sign_extend_ctl(0xFFFFFFFF), u64::MAX);
});
