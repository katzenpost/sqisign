//! Differential test of [`sqisign_hd::sample_random_index`] against the
//! C-derived vectors.
//!
//! Each record seeds a fresh KAT-only [`CtrDrbg`] with the recorded
//! 48-byte entropy block (no personalisation) and draws
//! `INDICES_PER_RECORD` consecutive indices. The bytes produced must
//! match the C reference's `sample_random_index` (file-local in
//! `the-sqisign/src/hd/ref/lvlx/theta_isogenies.c`) byte-for-byte.

use sqisign_common::CtrDrbg;
use sqisign_hd::sample_random_index;
use sqisign_vectors::{decode, load};

const VECTORS: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../vectors/hd/sample_random_index.json"
);

const INDICES_PER_RECORD: usize = 8;

fn fixed48(v: &[u8]) -> [u8; 48] {
    let mut a = [0u8; 48];
    assert_eq!(v.len(), 48, "entropy must be 48 bytes, got {}", v.len());
    a.copy_from_slice(v);
    a
}

#[test]
fn sample_random_index_matches_reference_vectors() {
    let f = load(VECTORS).expect("vector file load");
    assert_eq!(f.boundary, "sqisign_hd::sample_random_index");
    assert!(!f.vectors.is_empty(), "vector battery empty");

    for v in &f.vectors {
        let entropy = fixed48(&decode("entropy", &v.inputs["entropy"]).expect("entropy hex"));
        let indices_exp = decode("indices", &v.outputs["indices"]).expect("indices hex");
        assert_eq!(
            indices_exp.len(),
            INDICES_PER_RECORD,
            "vector {}: expected {} index bytes, got {}",
            v.id,
            INDICES_PER_RECORD,
            indices_exp.len()
        );

        let mut drbg = CtrDrbg::new(&entropy, None);
        let mut got = [0u8; INDICES_PER_RECORD];
        for slot in got.iter_mut() {
            *slot = sample_random_index(&mut drbg);
        }
        assert_eq!(
            &got[..],
            &indices_exp[..],
            "vector {}: index sequence diverged from C reference",
            v.id
        );
    }
}
