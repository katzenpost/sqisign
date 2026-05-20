//! Invariant fuzz target for `sqisign_gf::fp_set_one`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_set_one` on a fuzzing host.
//!
//! `fp_set_one` wraps the reference's `modone`, which writes positional
//! `1` and then calls `nres` to convert it to its Montgomery
//! representative. The on-the-wire output is therefore the Montgomery
//! `ONE`, `[0x19, 0, 0, 0, 0x300000000000]`, matching `extern const ONE`
//! at `fp_p5248_64.c:526..530`. The fuzz input is used as the
//! destination pre-fill so a no-op or partial-write port would leave a
//! visible residue. Linking the reference's `fp_p5248_64.c` for
//! byte-equality against C is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_set_one, Fp, NWORDS_FIELD};

/// Montgomery representative of `1`; matches `extern const ONE` at
/// `fp_p5248_64.c:526..530`.
const MONTGOMERY_ONE: Fp = [
    0x0000_0000_0000_0019,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_0000_0000_0000,
    0x0000_3000_0000_0000,
];

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
    if data.len() < NWORDS_FIELD * 8 {
        return;
    }
    let prefill = fp(&data[..NWORDS_FIELD * 8]);

    let mut out: Fp = prefill;
    fp_set_one(&mut out);
    assert_eq!(
        out, MONTGOMERY_ONE,
        "fp_set_one did not leave the destination as the Montgomery ONE"
    );
});
