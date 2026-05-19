//! Invariant fuzz target for `sqisign_mp::select_ct`: equals the per-bit
//! blend `(a & !mask) | (b & mask)`, and the 0 / all-ones endpoints
//! select a / b, for an arbitrary width. Linking `mp.c` for
//! byte-equality vs C is the next increment. Runner deferred (see
//! ../../../FUZZING.md).

#![no_main]

use libfuzzer_sys::fuzz_target;
use sqisign_mp::select_ct;

fuzz_target!(|data: &[u8]| {
    if data.len() < 9 {
        return;
    }
    let mask = {
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
    let a = to_limbs(&rest[..half]);
    if a.is_empty() {
        return;
    }
    let mut b = to_limbs(&rest[half..half + half]);
    b.resize(a.len(), 0);

    let mut c = vec![0u64; a.len()];
    select_ct(&mut c, &a, &b, mask);
    for i in 0..a.len() {
        assert_eq!(c[i], (a[i] & !mask) | (b[i] & mask), "not the bit blend");
    }

    let mut ca = vec![0u64; a.len()];
    select_ct(&mut ca, &a, &b, 0);
    assert_eq!(ca, a, "mask 0 must select a");
    let mut cb = vec![0u64; a.len()];
    select_ct(&mut cb, &a, &b, u64::MAX);
    assert_eq!(cb, b, "mask all-ones must select b");
});
