/*
 * sqisign.h: C ABI for level-1 SQIsign keypair generation, signing, and
 * signature verification.
 *
 * Mirrors the surface defined in crates/sqisign-ffi/src/lib.rs. Hand-written
 * (no cbindgen step) so the public ABI stays a small, stable, reviewable
 * contract. The crate produces both a cdylib (libsqisign_ffi.so on Linux,
 * libsqisign_ffi.dylib on macOS) and a staticlib (libsqisign_ffi.a); link
 * against either.
 *
 * Every entry point returns 1 on success and 0 on any failure. Failure
 * cases include: length mismatch on a fixed-size buffer, NULL pointer on
 * a non-empty buffer, the algorithm itself returning a non-success status,
 * and a Rust panic caught at the FFI boundary (panics cannot cross the
 * boundary; they are mapped to 0 via catch_unwind).
 *
 * License: GPL-3.0-or-later, matching the rest of the workspace.
 */

#ifndef SQISIGN_H
#define SQISIGN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Wire size of a serialized level-1 SQIsign public key, in bytes. */
#define SQISIGN_LVL1_PUBLIC_KEY_BYTES 65

/* Wire size of a serialized level-1 SQIsign secret key, in bytes. */
#define SQISIGN_LVL1_SECRET_KEY_BYTES 353

/* Wire size of a serialized level-1 SQIsign signature, in bytes. */
#define SQISIGN_LVL1_SIGNATURE_BYTES 148

/* Length of the entropy block consumed by keygen and sign. Matches the
 * NIST KAT format (48 bytes). The KAT-compatible CTR-DRBG is seeded from
 * this block. Production callers who do not want the CTR-DRBG should use
 * the Rust-level sqisign_sign::protocols_keygen / protocols_sign entry
 * points directly with their own RngSource implementation; this C ABI is
 * for KAT replay and for callers who supply hardware-derived entropy. */
#define SQISIGN_LVL1_ENTROPY_BYTES 48

/*
 * Generate a level-1 SQIsign keypair from a 48-byte entropy seed.
 *
 * pk, pk_len: output public key buffer; pk_len must equal
 *             SQISIGN_LVL1_PUBLIC_KEY_BYTES.
 * sk, sk_len: output secret key buffer; sk_len must equal
 *             SQISIGN_LVL1_SECRET_KEY_BYTES.
 * entropy, entropy_len: 48 bytes used to seed the KAT-compatible CTR-DRBG;
 *                       entropy_len must equal SQISIGN_LVL1_ENTROPY_BYTES.
 *
 * Returns 1 on success, 0 on any failure.
 *
 * The buffers are not retained past the call.
 */
int sqisign_lvl1_keygen(unsigned char *pk, size_t pk_len,
                         unsigned char *sk, size_t sk_len,
                         const unsigned char *entropy, size_t entropy_len);

/*
 * Sign a message with a level-1 SQIsign secret key, seeded by a 48-byte
 * entropy block. The output sig is the signature alone (length
 * SQISIGN_LVL1_SIGNATURE_BYTES); the NIST "sm = signature || msg"
 * concatenation is the caller's responsibility.
 *
 * sig, sig_len:  output signature buffer; sig_len must equal
 *                SQISIGN_LVL1_SIGNATURE_BYTES.
 * msg, msg_len:  input message. msg may be NULL iff msg_len == 0.
 * sk,  sk_len:   secret key; sk_len must equal SQISIGN_LVL1_SECRET_KEY_BYTES.
 * entropy, entropy_len: 48 bytes used to seed the signer's CTR-DRBG;
 *                       entropy_len must equal SQISIGN_LVL1_ENTROPY_BYTES.
 *
 * Returns 1 on success, 0 on any failure.
 */
int sqisign_lvl1_sign(unsigned char *sig, size_t sig_len,
                       const unsigned char *msg, size_t msg_len,
                       const unsigned char *sk, size_t sk_len,
                       const unsigned char *entropy, size_t entropy_len);

/*
 * Verify a level-1 SQIsign signature.
 *
 * sig, sig_len: signature bytes; sig_len must equal
 *               SQISIGN_LVL1_SIGNATURE_BYTES.
 * pk,  pk_len:  public key bytes; pk_len must equal
 *               SQISIGN_LVL1_PUBLIC_KEY_BYTES.
 * msg, msg_len: signed message. msg may be NULL iff msg_len == 0.
 *
 * Returns 1 on a valid signature, 0 on any failure (verification false,
 * length mismatch, NULL pointer where data was required, internal panic
 * caught at the FFI boundary).
 *
 * The buffers are not retained past the call.
 */
int sqisign_lvl1_verify(const unsigned char *sig, size_t sig_len,
                         const unsigned char *pk,  size_t pk_len,
                         const unsigned char *msg, size_t msg_len);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* SQISIGN_H */
