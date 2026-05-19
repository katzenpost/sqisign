//! Invariant fuzz target for `sqisign_mp::swap_ct`: the per-bit
//! conditional swap, with 0 / all-ones the no-op / full-swap endpoints,
//! and double-swap an involution, for an arbitrary width. Linking
//! `mp.c` for byte-equality vs C is the next increment. Runner deferred
//! (see ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::swap_ct;

fuzz_target!(|data: &[u8]| {
    if data.len() < 9 {
        return;
    }
    let option = {
        let mut w = [0u8; 8];
        w.copy_from_slice(&data[..8]);
        u64::from_le_bytes(w)
    };
    let rest = &data[8..];
    let half = rest.len() / 2;
    let to_limbs = |s: &[u8]| -> Vec<u64> {
        s.chunks(8)
            .map(|c| {
                let mut w = [0u8; 8];
                w[..c.len()].copy_from_slice(c);
                u64::from_le_bytes(w)
            })
            .collect()
    };
    let a0 = to_limbs(&rest[..half]);
    if a0.is_empty() {
        return;
    }
    let mut b0 = to_limbs(&rest[half..half + half]);
    b0.resize(a0.len(), 0);

    let mut a = a0.clone();
    let mut b = b0.clone();
    swap_ct(&mut a, &mut b, option);
    for i in 0..a0.len() {
        assert_eq!(a[i], (a0[i] & !option) | (b0[i] & option), "a not bit-swap");
        assert_eq!(b[i], (b0[i] & !option) | (a0[i] & option), "b not bit-swap");
    }
    // Involution.
    swap_ct(&mut a, &mut b, option);
    assert_eq!(a, a0, "double swap did not restore a");
    assert_eq!(b, b0, "double swap did not restore b");
});
