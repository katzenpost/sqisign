# sqisign Go binding

A `cgo` wrapper around the `sqisign-ffi` C ABI. Statically links
against `target/release/libsqisign_ffi.a`, which means the resulting
Go binary has no run-time dependency on the shared library: the
SQIsign code is folded straight into the executable.

## Build

From the repository root:

```sh
cargo build --release -p sqisign-ffi
cd bindings/go && go test ./sqisign/...
```

The first command produces `target/release/libsqisign_ffi.a`; the
`#cgo LDFLAGS:` directive in `sqisign/sqisign.go` picks it up by a
path relative to the binding's source directory.

The cgo header (`sqisign.h`) is vendored alongside the binding so
the include resolves cleanly whether the binding is consumed in-tree
or from Go's module cache. The companion staticlib path, however, is
still a `${SRCDIR}/../../../target/release/...` reference that only
exists for in-tree builds; consumers who pull this module via
`go get` and want to actually link against it must either copy the
prebuilt `libsqisign_ffi.a` into the expected location, override the
binding's `CGO_LDFLAGS` directive, or vendor the binding into their
own tree where the path resolves. Solving the staticlib-distribution
question more thoroughly (env-driven `-L`, pkg-config, or per-arch
release artifacts) is left as future work.

## Randomness

SQIsign is a randomised signature scheme; every keypair and every
signature consumes a stream of random bytes. The binding deliberately
offers no NIST CTR-DRBG path: keypair and signing entries take an
`io.Reader`, and every byte the algorithm demands is pulled from that
reader through a cgo callback (the highctidh `fillrandom` pattern).

Production callers pass `hpqc/rand.Reader`. Deterministic tests pass
something like `hpqc/rand.DeterministicRandReader` (a ChaCha20
keystream); the binding's own tests reproduce that construction
inline with `golang.org/x/crypto/chacha20` to keep its dependency
graph trivial.

The underlying C ABI also exposes entropy-block entries
(`sqisign_lvl1_keygen`, `sqisign_lvl1_sign`) that seed a CTR-DRBG.
They exist only so the Rust differential KAT tests can replay the
upstream NIST vectors bit-for-bit; they are intentionally not exported
through this binding.

## API

```go
import "github.com/katzenpost/sqisign/bindings/go/sqisign"

const (
    sqisign.PublicKeyBytes = 65
    sqisign.SecretKeyBytes = 353
    sqisign.SignatureBytes = 148
)

func KeyGen(rng io.Reader) (publicKey, secretKey []byte, err error)
func Sign(rng io.Reader, secretKey, msg []byte) (signature []byte, err error)
func Verify(signature, publicKey, msg []byte) (bool, error)
```

`Verify` returns `(false, nil)` when the signature is structurally
well-formed but does not validate; it returns `(false, err)` only
when a buffer has the wrong size (a caller-side bug).

`KeyGen` and `Sign` return `ErrNilRNG` for a nil reader and panic
through the cgo callback if the reader fails mid-call: the reference
algorithm has no recovery story for partial reads.

## Example

```go
package main

import (
    "fmt"
    "log"

    hpqcrand "github.com/katzenpost/hpqc/rand"

    "github.com/katzenpost/sqisign/bindings/go/sqisign"
)

func main() {
    pk, sk, err := sqisign.KeyGen(hpqcrand.Reader)
    if err != nil { log.Fatal(err) }

    msg := []byte("hello sqisign")
    sig, err := sqisign.Sign(hpqcrand.Reader, sk, msg)
    if err != nil { log.Fatal(err) }

    ok, err := sqisign.Verify(sig, pk, msg)
    if err != nil { log.Fatal(err) }
    fmt.Println("verified:", ok)
}
```

## Status and warning

This binding inherits every limitation of the Rust port. See
`SECURITY.md` at the repository root. SQIsign is a NIST Round 2
candidate and this implementation has not been audited; do not use
it in production systems where a forgery would cause harm.
