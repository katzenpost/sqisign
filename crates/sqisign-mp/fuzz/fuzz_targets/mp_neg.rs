//! Invariant fuzz target for `sqisign_mp::mp_neg`.
//!
//! Status: harness present, runner deferred (see ../../../FUZZING.md).
//! Asserts the faithful model (complement, +1 on limb 0, no carry) and
//! that it is true negation iff a[0] != 0, for an arbitrary width.
//! Linking `mp.c` for byte-equality vs C is the next increment.

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::{mp_add, mp_neg};

fuzz_target!(|data: &[u8]| {
    let v: Vec<u64> = data
        .chunks(8)
        .map(|c| {
            let mut w = [0u8; 8];
            w[..c.len()].copy_from_slice(c);
            u64::from_le_bytes(w)
        })
        .collect();
    if v.is_empty() {
        return;
    }

    let mut got = v.clone();
    mp_neg(&mut got);

    let mut model: Vec<u64> = v.iter().map(|x| !x).collect();
    model[0] = model[0].wrapping_add(1);
    assert_eq!(got, model, "mp_neg != faithful model");

    let mut sum = vec![0u64; v.len()];
    mp_add(&mut sum, &v, &got);
    let is_true_neg = sum.iter().all(|&w| w == 0);
    if v[0] != 0 {
        assert!(is_true_neg, "a[0]!=0 must give true negation");
    }
});
