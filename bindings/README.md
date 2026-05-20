# Language bindings

Bindings for the SQIsign level-1 signature scheme. The Rust crate
`sqisign-ffi` exports a small C ABI (three functions: keygen, sign,
verify), and these directories wrap that ABI in idiomatic Go and
Python.

| Language | Binding type                | Tests                           | README                  |
|----------|-----------------------------|---------------------------------|-------------------------|
| Go       | `cgo`, static link          | `bindings/go/sqisign/`          | [go/README.md](go/README.md) |
| Python   | `ctypes`, dynamic load      | `bindings/python/tests/`        | [python/README.md](python/README.md) |

Both bindings expect the workspace's FFI artefacts to be present at
`target/release/libsqisign_ffi.{a,so,dylib}`. Build them once with

```sh
cargo build --release -p sqisign-ffi
```

from the repository root and both bindings will find the library.

Neither binding adds API surface beyond what the C ABI already
exposes: the SQIsign level-1 byte sizes, plus `keygen(entropy)`,
`sign(secret_key, message, entropy)`, and `verify(signature,
public_key, message)`. Both translate the FFI's 1/0 success indicator
into idiomatic error reporting for their language.

The binding implementations are intentionally thin. They are not the
place to add hybrid constructions (Ed25519+SQIsign or similar), key
serialisation policies, or higher-level signature objects; those
concerns belong in the consumers (for Katzenpost: `hpqc` and the
dirauth signing path).

See `SECURITY.md` at the repository root before any production use.
The bindings inherit every threat-model caveat the Rust port has.
