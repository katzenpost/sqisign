"""ctypes loader and wrapper functions for libsqisign_ffi.

This module isolates everything that touches ctypes so the public
:mod:`sqisign` surface is small and the loading policy is documented in
one place. The loader resolution order is:

1. ``SQISIGN_FFI_LIB`` environment variable, if set, used verbatim.
2. ``../../target/release/libsqisign_ffi.<ext>`` relative to this file,
   where ``<ext>`` is ``so`` on Linux/BSD and ``dylib`` on macOS. This
   is the layout produced by ``cargo build --release -p sqisign-ffi``
   in the workspace at the repository root.
3. Bare ``libsqisign_ffi.so`` / ``libsqisign_ffi.dylib`` resolved via
   the operating system's default library search path (``LD_LIBRARY_PATH``
   on Linux, ``DYLD_LIBRARY_PATH`` on macOS).
"""

from __future__ import annotations

import ctypes
import os
import sys
from pathlib import Path

# Public byte sizes, mirroring the SQISIGN_LVL1_*_BYTES C macros in
# crates/sqisign-ffi/include/sqisign.h. These are kept in sync with
# the C ABI by inspection; if the FFI ever changes a size, this file
# is the first place the mismatch will surface (the round-trip test
# would fail with a length-mismatch error from the Rust side).
PUBLIC_KEY_BYTES = 65
SECRET_KEY_BYTES = 353
SIGNATURE_BYTES = 148
ENTROPY_BYTES = 48


class SqisignError(RuntimeError):
    """Raised on any failure from the FFI layer.

    The C ABI maps every failure mode to a 0 return value; this
    exception carries the originating call site as its message so the
    Python user can distinguish keygen failure from sign failure from
    a caller-side buffer-size bug.
    """


def _library_name() -> str:
    if sys.platform == "darwin":
        return "libsqisign_ffi.dylib"
    return "libsqisign_ffi.so"


def _candidate_paths() -> list[str]:
    env = os.environ.get("SQISIGN_FFI_LIB")
    if env:
        return [env]
    here = Path(__file__).resolve().parent
    libname = _library_name()
    # Repository layout: this file is at bindings/python/sqisign/_lib.py,
    # and the workspace target directory is at the repo root.
    workspace_target = here.parent.parent.parent / "target" / "release" / libname
    return [str(workspace_target), libname]


def _load() -> ctypes.CDLL:
    last_err: OSError | None = None
    for path in _candidate_paths():
        try:
            return ctypes.CDLL(path)
        except OSError as err:
            last_err = err
    raise SqisignError(
        "could not load libsqisign_ffi; "
        "run `cargo build --release -p sqisign-ffi` first, or set "
        "SQISIGN_FFI_LIB to the absolute path of the shared library. "
        f"Last error: {last_err}"
    )


_lib = _load()

# int sqisign_lvl1_keygen(unsigned char *pk, size_t pk_len,
#                        unsigned char *sk, size_t sk_len,
#                        const unsigned char *entropy, size_t entropy_len);
_lib.sqisign_lvl1_keygen.argtypes = [
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
]
_lib.sqisign_lvl1_keygen.restype = ctypes.c_int

# int sqisign_lvl1_sign(unsigned char *sig, size_t sig_len,
#                      const unsigned char *msg, size_t msg_len,
#                      const unsigned char *sk, size_t sk_len,
#                      const unsigned char *entropy, size_t entropy_len);
_lib.sqisign_lvl1_sign.argtypes = [
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
]
_lib.sqisign_lvl1_sign.restype = ctypes.c_int

# int sqisign_lvl1_verify(const unsigned char *sig, size_t sig_len,
#                        const unsigned char *pk,  size_t pk_len,
#                        const unsigned char *msg, size_t msg_len);
_lib.sqisign_lvl1_verify.argtypes = [
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
    ctypes.c_char_p,
    ctypes.c_size_t,
]
_lib.sqisign_lvl1_verify.restype = ctypes.c_int


def keygen(entropy: bytes) -> tuple[bytes, bytes]:
    """Generate a fresh SQIsign level-1 keypair.

    :param entropy: exactly :data:`ENTROPY_BYTES` (48) bytes used to
        seed the KAT-compatible NIST CTR-DRBG inside the FFI.
    :returns: ``(public_key, secret_key)`` of lengths
        :data:`PUBLIC_KEY_BYTES` and :data:`SECRET_KEY_BYTES`.
    :raises SqisignError: if the entropy is the wrong length or
        the FFI's keygen routine returned a non-success status.
    """
    if not isinstance(entropy, (bytes, bytearray)) or len(entropy) != ENTROPY_BYTES:
        raise SqisignError(f"entropy must be exactly {ENTROPY_BYTES} bytes")
    pk = ctypes.create_string_buffer(PUBLIC_KEY_BYTES)
    sk = ctypes.create_string_buffer(SECRET_KEY_BYTES)
    rc = _lib.sqisign_lvl1_keygen(
        pk, PUBLIC_KEY_BYTES,
        sk, SECRET_KEY_BYTES,
        bytes(entropy), ENTROPY_BYTES,
    )
    if rc != 1:
        raise SqisignError("sqisign_lvl1_keygen returned failure")
    return pk.raw[:PUBLIC_KEY_BYTES], sk.raw[:SECRET_KEY_BYTES]


def sign(secret_key: bytes, message: bytes, entropy: bytes) -> bytes:
    """Produce a SQIsign signature over ``message``.

    :param secret_key: exactly :data:`SECRET_KEY_BYTES` (353) bytes,
        as returned by :func:`keygen`.
    :param message: arbitrary-length payload; may be empty.
    :param entropy: exactly :data:`ENTROPY_BYTES` (48) bytes used to
        seed the signer's CTR-DRBG.
    :returns: signature of length :data:`SIGNATURE_BYTES` (148).
    :raises SqisignError: on any wrong-sized buffer or an FFI failure.
    """
    if not isinstance(secret_key, (bytes, bytearray)) or len(secret_key) != SECRET_KEY_BYTES:
        raise SqisignError(f"secret_key must be exactly {SECRET_KEY_BYTES} bytes")
    if not isinstance(entropy, (bytes, bytearray)) or len(entropy) != ENTROPY_BYTES:
        raise SqisignError(f"entropy must be exactly {ENTROPY_BYTES} bytes")
    if message is None:
        message = b""
    if not isinstance(message, (bytes, bytearray)):
        raise SqisignError("message must be bytes-like or None")
    sig = ctypes.create_string_buffer(SIGNATURE_BYTES)
    msg_buf = bytes(message)
    rc = _lib.sqisign_lvl1_sign(
        sig, SIGNATURE_BYTES,
        msg_buf if msg_buf else None, len(msg_buf),
        bytes(secret_key), SECRET_KEY_BYTES,
        bytes(entropy), ENTROPY_BYTES,
    )
    if rc != 1:
        raise SqisignError("sqisign_lvl1_sign returned failure")
    return sig.raw[:SIGNATURE_BYTES]


def verify(signature: bytes, public_key: bytes, message: bytes) -> bool:
    """Verify a SQIsign signature.

    :param signature: exactly :data:`SIGNATURE_BYTES` (148) bytes.
    :param public_key: exactly :data:`PUBLIC_KEY_BYTES` (65) bytes.
    :param message: arbitrary-length payload that was signed.
    :returns: ``True`` iff the signature is valid for the public key
        and message; ``False`` otherwise.
    :raises SqisignError: if any of the buffers has the wrong size.

    Note that a ``False`` return is not an error in the usual sense:
    the FFI verified the inputs and rejected the signature. Callers
    distinguishing "rejected" from "malformed input" should catch
    :class:`SqisignError` separately.
    """
    if not isinstance(signature, (bytes, bytearray)) or len(signature) != SIGNATURE_BYTES:
        raise SqisignError(f"signature must be exactly {SIGNATURE_BYTES} bytes")
    if not isinstance(public_key, (bytes, bytearray)) or len(public_key) != PUBLIC_KEY_BYTES:
        raise SqisignError(f"public_key must be exactly {PUBLIC_KEY_BYTES} bytes")
    if message is None:
        message = b""
    if not isinstance(message, (bytes, bytearray)):
        raise SqisignError("message must be bytes-like or None")
    msg_buf = bytes(message)
    rc = _lib.sqisign_lvl1_verify(
        bytes(signature), SIGNATURE_BYTES,
        bytes(public_key), PUBLIC_KEY_BYTES,
        msg_buf if msg_buf else None, len(msg_buf),
    )
    return rc == 1
