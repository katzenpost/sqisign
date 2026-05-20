# sqisign Go binding

A `cgo` wrapper around the `sqisign-ffi` C ABI. Statically links
against `target/release/libsqisign_ffi.a`, which means the resulting
Go binary has no run-time dependency on the shared library: the
SQIsign code is folded straight into the executable.

## Build

The binding is self-contained from cgo's point of view: both the C
header (`sqisign.h`) and a prebuilt `libsqisign_ffi.a` are vendored
alongside the .go files, the latter under `lib/<GOOS>_<GOARCH>/`.
Consumers pulling the module via `go get` link cleanly against the
vendored staticlib with no separate cargo step.

```sh
cd bindings/go && go test ./sqisign/...
```

### Supported platforms

The binding currently ships a vendored staticlib for `linux/amd64`
only. On every other GOOS/GOARCH the package compiles via a
pure-Go stub (`sqisign_unsupported.go`), `KeyGen`/`Sign`/`Verify`
return `ErrUnsupported`, and the constants stay valid for buffer
sizing and hybrid composition.

### Refreshing the vendored staticlib

After a Rust-side change to `sqisign-ffi`, regenerate the linux/amd64
staticlib from the repository root:

```sh
cargo build --release -p sqisign-ffi
cp target/release/libsqisign_ffi.a \
   bindings/go/sqisign/lib/linux_amd64/libsqisign_ffi.a
```

### Adding a new target

Build `libsqisign_ffi.a` for the target platform and drop it under
`bindings/go/sqisign/lib/<GOOS>_<GOARCH>/libsqisign_ffi.a`. Add the
matching `#cgo <GOOS>,<GOARCH> LDFLAGS: ...` directive in
`sqisign.go` and widen its `//go:build` constraint accordingly. The
stub's negated constraint then narrows automatically.

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
