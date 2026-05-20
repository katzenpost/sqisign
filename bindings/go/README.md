# sqisign Go binding

A `cgo` wrapper around the `sqisign-ffi` C ABI. Statically links
against `target/release/libsqisign_ffi.a`, which means the resulting
Go binary has no run-time dependency on the shared library — the
SQIsign code is folded straight into the executable.

## Build

From the repository root:

```sh
cargo build --release -p sqisign-ffi
cd bindings/go && go test ./sqisign/...
```

The first command produces `target/release/libsqisign_ffi.a`; the
`#cgo LDFLAGS:` directive in `sqisign/sqisign.go` picks it up by
relative path. If you move this binding outside the workspace, edit
that directive or set `CGO_LDFLAGS` to point at wherever you placed
the staticlib.

## API

```go
import "github.com/katzenpost/sqisign/bindings/go/sqisign"

const (
    sqisign.PublicKeyBytes = 65
    sqisign.SecretKeyBytes = 353
    sqisign.SignatureBytes = 148
    sqisign.EntropyBytes   = 48
)

func KeyGen(entropy []byte) (publicKey, secretKey []byte, err error)
func Sign(secretKey, msg, entropy []byte) (signature []byte, err error)
func Verify(signature, publicKey, msg []byte) (bool, error)
```

`Verify` returns `(false, nil)` when the signature is structurally
well-formed but does not validate; it returns `(false, err)` only
when a buffer has the wrong size (a caller-side bug).

## Example

```go
entropy := make([]byte, sqisign.EntropyBytes)
if _, err := rand.Read(entropy); err != nil {
    log.Fatal(err)
}
pk, sk, err := sqisign.KeyGen(entropy)
if err != nil { log.Fatal(err) }

msg := []byte("hello sqisign")
sig, err := sqisign.Sign(sk, msg, entropy)
if err != nil { log.Fatal(err) }

ok, err := sqisign.Verify(sig, pk, msg)
if err != nil { log.Fatal(err) }
fmt.Println("verified:", ok)
```

## Status and warning

This binding inherits every limitation of the Rust port. See
`SECURITY.md` at the repository root. SQIsign is a NIST Round 2
candidate and this implementation has not been audited; do not use
it in production systems where a forgery would cause harm.
