//! Invariant fuzz target for `sqisign_mp::mp_copy`: it is the identity
//! into `b`, for an arbitrary width. Linking `mp.c` for byte-equality vs
//! C is the next increment. Runner deferred (see ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_copy;

fuzz_target!(|data: &[u8]| {
    let a: Vec<u64> = data
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    let mut b = vec![0xa5a5_a5a5_a5a5_a5a5u64; a.len()];
    mp_copy(&mut b, &a);
    assert_eq!(b, a, "mp_copy is not the identity");
});
