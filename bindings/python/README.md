# sqisign Python binding

A `ctypes` wrapper around the `sqisign-ffi` C ABI. Loads the cdylib at
import time; no compiler is needed at install time, only the shared
library at run time.

## Build

From the repository root:

```sh
cargo build --release -p sqisign-ffi
cd bindings/python && python3 -m pytest tests/
```

The first command produces `target/release/libsqisign_ffi.so` (Linux)
or `libsqisign_ffi.dylib` (macOS); the loader in
`sqisign/_lib.py` picks it up by relative path. If you install the
shared library elsewhere, set `SQISIGN_FFI_LIB=/path/to/the/library`
in the environment before importing the package.

## API

```python
import sqisign

sqisign.PUBLIC_KEY_BYTES  # 65
sqisign.SECRET_KEY_BYTES  # 353
sqisign.SIGNATURE_BYTES   # 148
sqisign.ENTROPY_BYTES     # 48

public_key, secret_key = sqisign.keygen(entropy: bytes)
signature              = sqisign.sign(secret_key, message, entropy)
ok                     = sqisign.verify(signature, public_key, message)
```

`verify` returns ``True`` or ``False``; wrong-sized buffers raise
``sqisign.SqisignError`` rather than returning ``False``, so callers
distinguishing "rejected" from "malformed input" can catch the
exception separately.

## Example

```python
import secrets, sqisign

entropy = secrets.token_bytes(sqisign.ENTROPY_BYTES)
public_key, secret_key = sqisign.keygen(entropy)

message = b"hello sqisign"
signature = sqisign.sign(secret_key, message, entropy)

assert sqisign.verify(signature, public_key, message)
```

## Status and warning

This binding inherits every limitation of the Rust port. See
`SECURITY.md` at the repository root. SQIsign is a NIST Round 2
candidate and this implementation has not been audited; do not use
it in production systems where a forgery would cause harm.
