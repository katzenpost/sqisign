//! Invariant fuzz target for `sqisign_mp::mp_is_one`: equals the
//! canonical-one predicate, for an arbitrary non-empty width. Linking
//! `mp.c` for byte-equality vs C is the next increment. Runner deferred
//! (see ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_is_one;

fuzz_target!(|data: &[u8]| {
    let x: Vec<u64> = data
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    if x.is_empty() {
        return;
    }
    assert_eq!(
        mp_is_one(&x),
        x[0] == 1 && x[1..].iter().all(|&t| t == 0),
        "mp_is_one is not the canonical-one predicate"
    );
});
