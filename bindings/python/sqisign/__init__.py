"""Python binding for the SQIsign level-1 signature scheme.

Wraps the cdylib emitted by

    cargo build --release -p sqisign-ffi

at ``../../target/release/libsqisign_ffi.so`` (Linux) or
``libsqisign_ffi.dylib`` (macOS), resolved relative to this package.
Callers who install the shared library elsewhere can override the
search path by setting the ``SQISIGN_FFI_LIB`` environment variable to
the absolute path of the library before importing this module.

The binding is a thin ``ctypes`` wrapper: no compiler is needed at
install time, only the shared library at run time. Every entry point
in the FFI returns 1 on success and 0 on any failure (the Rust side
catches panics so they cannot cross the boundary); the Python wrappers
translate that into ``(bytes, ...)`` returns or raise
:class:`SqisignError` for caller bugs and algorithmic failures.

SQIsign is a NIST Round 2 candidate. This library has not been
audited. See SECURITY.md at the repository root before any production
use.
"""

from ._lib import (
    ENTROPY_BYTES,
    PUBLIC_KEY_BYTES,
    SECRET_KEY_BYTES,
    SIGNATURE_BYTES,
    SqisignError,
    keygen,
    sign,
    verify,
)

__all__ = [
    "ENTROPY_BYTES",
    "PUBLIC_KEY_BYTES",
    "SECRET_KEY_BYTES",
    "SIGNATURE_BYTES",
    "SqisignError",
    "keygen",
    "sign",
    "verify",
]

__version__ = "0.1.0a0"
