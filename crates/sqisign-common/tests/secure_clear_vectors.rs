//! Differential test of the ported `secure_clear` against the committed
//! C-derived vectors.
//!
//! Each vector records an input buffer, a `clear_len`, and the buffer as
//! the reference left it after `sqisign_secure_clear(buf, clear_len)`. The
//! replay clears exactly the leading `clear_len` bytes and bit-compares the
//! whole buffer: this proves both that the wiped span is zero and that the
//! tail beyond `clear_len` is untouched (the reference clears precisely
//! `size` bytes), the only observable contract the boundary has.

use sqisign_common::secure_clear;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/common/secure_clear.json"
);

fn le_u32(b: &[u8]) -> usize {
    assert_eq!(b.len(), 4, "clear_len must be a u32");
    u32::from_le_bytes([b[0], b[1], b[2], b[3]]) as usize
}

#[test]
fn secure_clear_matches_reference_vectors() {
    let file = load(VECTORS).expect("vector file must load, be canonical, and match the pin");
    assert_eq!(file.boundary, "sqisign_common::secure_clear");
    assert!(
        file.vectors.len() >= 1000,
        "expected the full battery, found {} vectors",
        file.vectors.len()
    );

    for v in &file.vectors {
        let mut buf = decode("input", &v.inputs["input"]).expect("input hex");
        let clear_len =
            le_u32(&decode("clear_len", &v.inputs["clear_len"]).expect("clear_len hex"));
        let expected = decode("output", &v.outputs["output"]).expect("output hex");

        assert!(
            clear_len <= buf.len(),
            "vector {}: clear_len > buffer",
            v.id
        );
        secure_clear(&mut buf[..clear_len]);
        assert_eq!(
            buf,
            expected,
            "vector {} diverged from the C reference (len={}, clear_len={})",
            v.id,
            expected.len(),
            clear_len
        );
    }
}
