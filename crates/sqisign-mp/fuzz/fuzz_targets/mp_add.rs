//! Invariant fuzz target for `sqisign_mp::mp_add`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Ready for `cargo +nightly fuzz run mp_add` on a fuzzing host.
//!
//! Asserts the ring laws (commutativity, zero identity) and that a single
//! limb matches native wrapping addition, for an arbitrary limb count.
//! Linking the reference's `mp.c` for byte-equality against C is the
//! documented next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_add;

fn limbs(bytes: &[u8]) -> Vec<u64> {
    bytes
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect()
}

fuzz_target!(|data: &[u8]| {
    // Split the input into two equal-length limb vectors.
    let half = data.len() / 2;
    let a = limbs(&data[..half]);
    let mut b = limbs(&data[half..half + half]);
    if a.is_empty() {
        return;
    }
    b.resize(a.len(), 0);

    let mut c1 = vec![0u64; a.len()];
    let mut c2 = vec![0u64; a.len()];
    mp_add(&mut c1, &a, &b);
    mp_add(&mut c2, &b, &a);
    assert_eq!(c1, c2, "mp_add not commutative");

    let zero = vec![0u64; a.len()];
    let mut c3 = vec![0u64; a.len()];
    mp_add(&mut c3, &a, &zero);
    assert_eq!(c3, a, "zero is not the additive identity");

    if a.len() == 1 {
        assert_eq!(c1[0], a[0].wrapping_add(b[0]), "single limb != wrapping_add");
    }
});
