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
 * Keypair generation and signing come in two flavours; both produce the
 * same wire output given the same byte stream, the difference is which
 * RNG drives them:
 *
 *   - sqisign_lvl1_keygen / sqisign_lvl1_sign take a 48-byte entropy
 *     block, seed a NIST AES-256 CTR-DRBG with it, and drive the
 *     algorithm from there. These exist so the Rust differential KAT
 *     tests can replay the upstream NIST vectors bit-for-bit; the
 *     Katzenpost research team does not use NIST's DRBG in production
 *     and non-Rust callers should not use these entries.
 *
 *   - sqisign_lvl1_keygen_with_rng / sqisign_lvl1_sign_with_rng take a
 *     (callback, context) pair and pull every byte of randomness from
 *     the callback. No NIST DRBG, no hidden state. This is the
 *     production path for non-Rust callers (the cgo binding in
 *     bindings/go uses it to thread a Go io.Reader through; future
 *     Python bindings use it analogously).
 *
 * Every entry point returns 1 on success and 0 on any failure. Failure
 * cases include: length mismatch on a fixed-size buffer, NULL pointer
 * on a non-empty buffer or a NULL callback function pointer, the
 * algorithm itself returning a non-success status, and a Rust panic
 * caught at the FFI boundary (panics cannot cross the boundary; they
 * are mapped to 0 via catch_unwind).
 *
 * License: GPL-3.0-or-later, matching the rest of the workspace.
 */

#ifndef SQISIGN_H
#define SQISIGN_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Wire size of a serialized level-1 SQIsign public key, in bytes. */
#define SQISIGN_LVL1_PUBLIC_KEY_BYTES 65

/* Wire size of a serialized level-1 SQIsign secret key, in bytes. */
#define SQISIGN_LVL1_SECRET_KEY_BYTES 353

/* Wire size of a serialized level-1 SQIsign signature, in bytes. */
#define SQISIGN_LVL1_SIGNATURE_BYTES 148

/* Length of the entropy block consumed by sqisign_lvl1_keygen /
 * sqisign_lvl1_sign. Matches the NIST KAT format (48 bytes); the
 * KAT-compatible CTR-DRBG is seeded from this block. These entries
 * exist so the Rust differential KAT tests can replay the upstream
 * NIST vectors bit-for-bit; the Katzenpost research team does not use
 * NIST's DRBG in production. Non-Rust callers should use the
 * _with_rng entries below. */
#define SQISIGN_LVL1_ENTROPY_BYTES 48

/*
 * Caller-supplied-RNG callback signature. Implementations must, on
 * each call, fill exactly `len` bytes at `out` from the caller's
 * randomness source. `context` is an opaque uintptr_t the FFI threads
 * through unchanged; callers typically use it to carry a pointer to
 * RNG state through the C boundary (e.g. a Go go-pointer handle, a
 * Python self-pointer).
 *
 * The shape mirrors highctidh's ctidh_fillrandom callback, so the
 * same `mattn/go-pointer`-style context-smuggling works here.
 */
typedef void (*sqisign_fill_random_fn)(unsigned char *out, size_t len,
                                       uintptr_t context);

/*
 * KAT-replay keypair entry: drive the algorithm from a CTR-DRBG seeded
 * by the supplied 48-byte entropy block.
 *
 * Non-Rust callers should use sqisign_lvl1_keygen_with_rng instead.
 *
 * Returns 1 on success, 0 on any failure. Buffers are not retained
 * past the call.
 */
int sqisign_lvl1_keygen(unsigned char *pk, size_t pk_len,
                         unsigned char *sk, size_t sk_len,
                         const unsigned char *entropy, size_t entropy_len);

/*
 * KAT-replay signing entry: drive the algorithm from a CTR-DRBG seeded
 * by the supplied 48-byte entropy block. The output sig is the
 * signature alone (length SQISIGN_LVL1_SIGNATURE_BYTES); the NIST
 * "sm = signature || msg" concatenation is the caller's
 * responsibility.
 *
 * Non-Rust callers should use sqisign_lvl1_sign_with_rng instead.
 *
 * Returns 1 on success, 0 on any failure.
 */
int sqisign_lvl1_sign(unsigned char *sig, size_t sig_len,
                       const unsigned char *msg, size_t msg_len,
                       const unsigned char *sk, size_t sk_len,
                       const unsigned char *entropy, size_t entropy_len);

/*
 * Production keypair entry: drive the algorithm from a caller-supplied
 * RNG callback.
 *
 * pk, pk_len: output public key buffer; pk_len must equal
 *             SQISIGN_LVL1_PUBLIC_KEY_BYTES.
 * sk, sk_len: output secret key buffer; sk_len must equal
 *             SQISIGN_LVL1_SECRET_KEY_BYTES.
 * fill_random: function pointer the algorithm calls every time it
 *              needs randomness. Must be non-NULL.
 * rng_context: opaque value passed through to every fill_random call.
 *              Not interpreted by this library.
 *
 * Returns 1 on success, 0 on any failure (including a NULL
 * fill_random).
 */
int sqisign_lvl1_keygen_with_rng(unsigned char *pk, size_t pk_len,
                                  unsigned char *sk, size_t sk_len,
                                  sqisign_fill_random_fn fill_random,
                                  uintptr_t rng_context);

/*
 * Production signing entry: drive the algorithm from a caller-supplied
 * RNG callback. The output sig is the signature alone (length
 * SQISIGN_LVL1_SIGNATURE_BYTES).
 *
 * Returns 1 on success, 0 on any failure.
 */
int sqisign_lvl1_sign_with_rng(unsigned char *sig, size_t sig_len,
                                const unsigned char *msg, size_t msg_len,
                                const unsigned char *sk, size_t sk_len,
                                sqisign_fill_random_fn fill_random,
                                uintptr_t rng_context);

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
