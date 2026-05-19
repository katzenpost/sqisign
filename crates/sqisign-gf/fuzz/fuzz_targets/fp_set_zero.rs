//! Invariant fuzz target for `sqisign_gf::fp_set_zero`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run fp_set_zero` on a fuzzing host.
//!
//! `fp_set_zero` is the reference's `modzer`: a plain five-limb
//! zero-fill. The asserted invariant is therefore the only one that
//! transfers from the bit-exact form, namely that the destination ends
//! up as the canonical all-zero limb vector on arbitrary (possibly
//! non-canonical) pre-fills. The fuzz input is used as the destination
//! pre-fill so a no-op or partial-write port would leave a visible
//! non-zero residue. Linking the reference's `fp_p5248_64.c` for
//! byte-equality against C is the documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_gf::{fp_set_zero, Fp, NWORDS_FIELD};

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
    fp_set_zero(&mut out);
    assert_eq!(
        out,
        [0u64; NWORDS_FIELD],
        "fp_set_zero did not leave the destination as the canonical all-zero limb vector"
    );
});
