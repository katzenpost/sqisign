//! Randomness abstraction and the NIST AES-256 CTR-DRBG.
//!
//! The [`RngSource`] trait is the only RNG surface every other crate in
//! this workspace knows about. Every signing-path primitive that needs
//! randomness takes `&mut impl RngSource` and never sees a concrete
//! implementation. Production callers wire in whatever they like (a
//! future Rust port of the Katzenpost hpqc/rand library, an `OsRng`
//! adapter, a hardware RNG, ...).
//!
//! [`CtrDrbg`] is the NIST AES-256-CTR_DRBG (SP 800-90A, no derivation
//! function, no prediction resistance) as the SQIsign reference uses it
//! to drive the KAT test vectors. It is a **KAT-only** primitive in this
//! workspace: it exists exclusively so the differential tests against
//! `vendor/the-sqisign` can replay the upstream KAT seeds bit-for-bit.
//! The Katzenpost research team does not use NIST's DRBG in production
//! (its backdoor history is the reason); production builds inject a
//! different [`RngSource`] and never construct a [`CtrDrbg`].
//!
//! The block-cipher core, AES-256 in ECB on a single block, is the
//! audited RustCrypto `aes` crate rather than a bespoke reimplementation:
//! AES is a standardized primitive and the port's rule is wire-and-prove,
//! not reinvent. The DRBG construction *around* it is ported here and
//! proven bit-equal to the reference by the committed C-derived vectors
//! (`tests/ctr_drbg_vectors.rs`).

use aes::cipher::{Array, BlockCipherEncrypt, KeyInit};
use aes::Aes256;

/// The single abstraction over a byte source. Every RNG-driven primitive
/// in this workspace (`ibz_rand_*`, `quat_represent_integer`, the
/// randomized theta splitting, ...) takes `&mut impl RngSource`.
///
/// Mirrors the C reference's `randombytes(unsigned char *, size_t)` thunk
/// in shape; the concrete RNG state is held by the caller and threaded
/// explicitly, which is the difference that lets production swap NIST's
/// CTR-DRBG out for something the project actually trusts.
pub trait RngSource {
    /// Fill `out` with fresh random bytes. Implementations must succeed
    /// or panic; the SQIsign reference signature treats RNG failure as
    /// unrecoverable.
    fn fill(&mut self, out: &mut [u8]);
}

/// AES-256-ECB of one 16-byte block, matching the reference's
/// `AES256_ECB(key, ctr, buffer)` (a single-block encrypt).
fn aes256_ecb(key: &[u8; 32], block: &[u8; 16]) -> [u8; 16] {
    let cipher = Aes256::new(&Array(*key));
    let mut b = Array(*block);
    cipher.encrypt_block(&mut b);
    b.0
}

/// Increment the 128-bit big-endian counter `v` by one, exactly as the
/// reference's inlined loop (from the last byte, carrying on 0xff).
fn increment(v: &mut [u8; 16]) {
    for byte in v.iter_mut().rev() {
        if *byte == 0xff {
            *byte = 0x00;
        } else {
            *byte += 1;
            break;
        }
    }
}

/// An AES-256 CTR-DRBG instance. Construct with [`CtrDrbg::new`] (the
/// reference's `randombytes_init`), then draw bytes via the
/// [`RngSource::fill`] impl (its `randombytes`). The state evolves after
/// every draw.
///
/// **KAT-only.** Production builds wire a different [`RngSource`]; this
/// type exists in the workspace so the upstream KAT seeds can be replayed
/// in differential tests.
///
/// `Clone` is intentional: bench harnesses snapshot the post-keygen
/// state and replay it per iteration so the signer always sees the
/// byte stream it would see after a real keypair generation.
#[derive(Clone)]
pub struct CtrDrbg {
    key: [u8; 32],
    v: [u8; 16],
}

impl CtrDrbg {
    /// `AES256_CTR_DRBG_Update`: derive 48 fresh bytes by encrypting three
    /// successive counter values, optionally XOR the provided data, then
    /// resplit into the new key and counter.
    fn update(&mut self, provided: Option<&[u8; 48]>) {
        let mut temp = [0u8; 48];
        for chunk in temp.chunks_exact_mut(16) {
            increment(&mut self.v);
            chunk.copy_from_slice(&aes256_ecb(&self.key, &self.v));
        }
        if let Some(pd) = provided {
            for (t, p) in temp.iter_mut().zip(pd.iter()) {
                *t ^= *p;
            }
        }
        self.key.copy_from_slice(&temp[..32]);
        self.v.copy_from_slice(&temp[32..]);
    }

    /// `randombytes_init`: seed material is the entropy XOR the optional
    /// personalization string; key and counter start at zero, then one
    /// update folds the seed in.
    pub fn new(entropy: &[u8; 48], personalization: Option<&[u8; 48]>) -> Self {
        let mut seed = *entropy;
        if let Some(p) = personalization {
            for (s, p) in seed.iter_mut().zip(p.iter()) {
                *s ^= *p;
            }
        }
        let mut drbg = CtrDrbg {
            key: [0u8; 32],
            v: [0u8; 16],
        };
        drbg.update(Some(&seed));
        drbg
    }
}

impl RngSource for CtrDrbg {
    /// `randombytes`: fill `out` from the counter stream, then run an update
    /// with no provided data. A zero-length request still performs that
    /// update, advancing the state, exactly as the reference does.
    fn fill(&mut self, out: &mut [u8]) {
        let mut off = 0;
        while off < out.len() {
            increment(&mut self.v);
            let block = aes256_ecb(&self.key, &self.v);
            let take = core::cmp::min(16, out.len() - off);
            out[off..off + take].copy_from_slice(&block[..take]);
            off += take;
        }
        self.update(None);
    }
}

// Convenience: `&mut R` is itself an `RngSource` if `R` is. This lets
// callers thread an existing RNG through nested helper calls without
// reborrow ceremony.
impl<R: RngSource + ?Sized> RngSource for &mut R {
    fn fill(&mut self, out: &mut [u8]) {
        (**self).fill(out)
    }
}
