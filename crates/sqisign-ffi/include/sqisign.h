/*
 * sqisign.h: C ABI for level-1 SQIsign signature verification.
 *
 * This header mirrors the surface defined in crates/sqisign-ffi/src/lib.rs.
 * It is hand-written (no cbindgen step) to keep the public ABI a small,
 * stable, and reviewable contract: two size constants and one verify
 * function. The crate produces both a cdylib (libsqisign_ffi.so) and a
 * staticlib (libsqisign_ffi.a); link against either.
 *
 * The verifier is verify-only by design. Signing will be exposed once the
 * LLL / dpe quaternion paths land in the underlying Rust port.
 *
 * License: LGPL-3.0-or-later, matching the rest of the workspace.
 */

#ifndef SQISIGN_H
#define SQISIGN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Wire size of a serialized level-1 SQIsign public key, in bytes. */
#define SQISIGN_LVL1_PUBLIC_KEY_BYTES 65

/* Wire size of a serialized level-1 SQIsign signature, in bytes. */
#define SQISIGN_LVL1_SIGNATURE_BYTES 148

/*
 * Verify a level-1 SQIsign signature.
 *
 * Parameters:
 *   sig, sig_len: signature bytes; sig_len must equal
 *                 SQISIGN_LVL1_SIGNATURE_BYTES.
 *   pk,  pk_len:  public key bytes; pk_len must equal
 *                 SQISIGN_LVL1_PUBLIC_KEY_BYTES.
 *   msg, msg_len: signed message. msg may be NULL iff msg_len == 0.
 *
 * Return value:
 *   1 on a valid signature.
 *   0 on any failure (verification false, length mismatch, NULL pointer
 *     where data was required, or an internal panic that was caught at
 *     the FFI boundary). No other return value is possible.
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
