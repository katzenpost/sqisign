//! Invariant fuzz target for `sqisign_gf::fp_set_one`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_set_one` on a fuzzing host.
//!
//! `fp_set_one` is the reference's `modone`: a plain five-limb fill
//! writing `1` to limb `0` and `0` to limbs `1..5`. The asserted
//! invariant is the only one that transfers from the bit-exact form,
//! namely that the destination ends up as the canonical positional-one
//! limb vector `[1, 0, 0, 0, 0]` on arbitrary (possibly non-canonical)
//! pre-fills. The fuzz input is used as the destination pre-fill so a
//! no-op or partial-write port would leave a visible residue. Note that
//! this is the **positional** one (`1 = 1 * 2^0`), not the Montgomery
//! representative of one (`ONE = {0x19, 0, 0, 0, 0x300000000000}`).
//! Linking the reference's `fp_p5248_64.c` for byte-equality against C
//! is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_set_one, Fp, NWORDS_FIELD};

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
        out,
        [1u64, 0, 0, 0, 0],
        "fp_set_one did not leave the destination as the canonical positional-one limb vector"
    );
});
