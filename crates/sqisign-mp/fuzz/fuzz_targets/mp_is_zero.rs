//! Invariant fuzz target for `sqisign_mp::mp_is_zero`: equals the
//! all-limbs-zero predicate, for an arbitrary width. Linking `mp.c` for
//! byte-equality vs C is the next increment. Runner deferred (see
//! ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::mp_is_zero;

fuzz_target!(|data: &[u8]| {
    let a: Vec<u64> = data
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    assert_eq!(
        mp_is_zero(&a),
        a.iter().all(|&x| x == 0),
        "mp_is_zero is not the all-zero predicate"
    );
});
