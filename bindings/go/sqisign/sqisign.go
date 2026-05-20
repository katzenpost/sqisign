// Package sqisign is a Go binding for the SQIsign level-1 signature
// scheme, wrapping the sqisign-ffi crate at the repository root.
//
// The binding links against the staticlib emitted by
//
//	cargo build --release -p sqisign-ffi
//
// at ../../../target/release/libsqisign_ffi.a, resolved relative to
// this source file. Users who keep this binding inside the sqisign
// workspace get the link automatically; users who vendor the binding
// elsewhere need to either copy the staticlib in or override CGO_LDFLAGS.
//
// The C ABI guarantees that Rust panics never cross the boundary; every
// entry point either returns 1 for success or 0 for any failure. The
// Go wrappers translate that into idiomatic (result, error) pairs.
//
// Randomness model. SQIsign is a randomised signature scheme; every
// keypair and every signature consumes a stream of random bytes. The
// Go binding deliberately offers no NIST CTR-DRBG path: keypair and
// signing entries take an io.Reader, and every byte the algorithm
// demands is pulled from that reader through a cgo callback. Production
// callers pass hpqc/rand's Reader; deterministic tests pass a
// chacha20-keystream reader (see hpqc/rand.DeterministicRandReader).
// The KAT-replay entries in the underlying C ABI exist for Rust-side
// differential testing and are intentionally not exported here.
//
// SQIsign is a NIST Round 2 candidate. This library has not been
// audited. See SECURITY.md at the repository root before any
// production use.
package sqisign

// #cgo LDFLAGS: ${SRCDIR}/../../../target/release/libsqisign_ffi.a -lm -ldl -lpthread
// #include "sqisign.h"
//
// /* The Go-exported symbol the Rust FFI calls back into when it needs
//    randomness. Declared here so we can pass its address as a
//    sqisign_fill_random_fn function pointer. The implementation lives
//    in Go (see fillrandom.go), where it can restore the io.Reader the
//    caller saved with gopointer.Save(). */
// extern void sqisign_go_fillrandom(unsigned char *out, size_t len, uintptr_t context);
import "C"

import (
	"errors"
	"io"
	"unsafe"

	gopointer "github.com/mattn/go-pointer"
)

// Wire sizes of the SQIsign level-1 byte buffers. These match the C
// ABI constants in sqisign.h byte-for-byte and the
// SQISIGN_LVL1_*_BYTES Rust constants in crates/sqisign-ffi/src/lib.rs.
const (
	PublicKeyBytes = C.SQISIGN_LVL1_PUBLIC_KEY_BYTES
	SecretKeyBytes = C.SQISIGN_LVL1_SECRET_KEY_BYTES
	SignatureBytes = C.SQISIGN_LVL1_SIGNATURE_BYTES
)

// Errors returned by this package. The C ABI maps every failure mode to
// a single 0 return; we surface a small handful of Go errors that
// distinguish "caller passed a wrong-sized buffer" (a programming bug)
// from "the algorithm itself failed" (rare, but possible because both
// keygen and sign internally loop until success in the reference) and
// from "no rng was supplied".
var (
	ErrNilRNG           = errors.New("sqisign: rng must not be nil")
	ErrInvalidPublicKey = errors.New("sqisign: public key must be exactly PublicKeyBytes long")
	ErrInvalidSecretKey = errors.New("sqisign: secret key must be exactly SecretKeyBytes long")
	ErrInvalidSignature = errors.New("sqisign: signature must be exactly SignatureBytes long")
	ErrKeygenFailed     = errors.New("sqisign: keygen returned failure")
	ErrSignFailed       = errors.New("sqisign: sign returned failure")
)

// KeyGen produces a fresh SQIsign level-1 keypair, drawing every byte
// of randomness from rng. The returned public key has length
// PublicKeyBytes; the secret key has length SecretKeyBytes.
//
// rng is held only for the duration of the call; the algorithm may
// invoke rng.Read many times. A read returning fewer bytes than
// requested, or a non-nil error, panics inside the cgo callback (the
// reference algorithm treats RNG failure as unrecoverable). Callers
// who want a non-panicking surface should wrap their reader with one
// that retries internally.
func KeyGen(rng io.Reader) (publicKey, secretKey []byte, err error) {
	if rng == nil {
		return nil, nil, ErrNilRNG
	}
	pk := make([]byte, PublicKeyBytes)
	sk := make([]byte, SecretKeyBytes)
	ctx := gopointer.Save(rng)
	defer gopointer.Unref(ctx)

	rc := C.sqisign_lvl1_keygen_with_rng(
		(*C.uchar)(unsafe.Pointer(&pk[0])), C.size_t(len(pk)),
		(*C.uchar)(unsafe.Pointer(&sk[0])), C.size_t(len(sk)),
		C.sqisign_fill_random_fn(C.sqisign_go_fillrandom),
		C.uintptr_t(uintptr(ctx)),
	)
	if rc != 1 {
		return nil, nil, ErrKeygenFailed
	}
	return pk, sk, nil
}

// Sign produces a SQIsign signature over msg using secretKey, drawing
// every byte of randomness from rng. The returned signature has length
// SignatureBytes. msg may be nil if its length is also 0.
//
// The NIST "sm = signature || msg" concatenation is the caller's
// responsibility; this entry point returns the signature alone.
func Sign(rng io.Reader, secretKey, msg []byte) ([]byte, error) {
	if rng == nil {
		return nil, ErrNilRNG
	}
	if len(secretKey) != SecretKeyBytes {
		return nil, ErrInvalidSecretKey
	}
	sig := make([]byte, SignatureBytes)
	var msgPtr *C.uchar
	if len(msg) != 0 {
		msgPtr = (*C.uchar)(unsafe.Pointer(&msg[0]))
	}
	ctx := gopointer.Save(rng)
	defer gopointer.Unref(ctx)

	rc := C.sqisign_lvl1_sign_with_rng(
		(*C.uchar)(unsafe.Pointer(&sig[0])), C.size_t(len(sig)),
		msgPtr, C.size_t(len(msg)),
		(*C.uchar)(unsafe.Pointer(&secretKey[0])), C.size_t(len(secretKey)),
		C.sqisign_fill_random_fn(C.sqisign_go_fillrandom),
		C.uintptr_t(uintptr(ctx)),
	)
	if rc != 1 {
		return nil, ErrSignFailed
	}
	return sig, nil
}

// Verify returns true iff the signature is valid for msg under the
// given public key. A non-nil error indicates a caller-side bug (a
// wrong-sized buffer); a false return with no error indicates the
// signature did not verify.
func Verify(signature, publicKey, msg []byte) (bool, error) {
	if len(signature) != SignatureBytes {
		return false, ErrInvalidSignature
	}
	if len(publicKey) != PublicKeyBytes {
		return false, ErrInvalidPublicKey
	}
	var msgPtr *C.uchar
	if len(msg) != 0 {
		msgPtr = (*C.uchar)(unsafe.Pointer(&msg[0]))
	}
	rc := C.sqisign_lvl1_verify(
		(*C.uchar)(unsafe.Pointer(&signature[0])), C.size_t(len(signature)),
		(*C.uchar)(unsafe.Pointer(&publicKey[0])), C.size_t(len(publicKey)),
		msgPtr, C.size_t(len(msg)),
	)
	return rc == 1, nil
}
