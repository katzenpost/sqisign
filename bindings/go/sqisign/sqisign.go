// Package sqisign is a Go binding for the SQIsign level-1 signature
// scheme, wrapping the sqisign-ffi crate at the repository root.
//
// The binding links against the staticlib emitted by
//
//	cargo build --release -p sqisign-ffi
//
// at ../../../target/release/libsqisign_ffi.a, resolved relative to
// this source file. Users who keep this binding inside the sqisign-rs
// workspace get the link automatically; users who vendor the binding
// elsewhere need to either copy the staticlib in or override CGO_LDFLAGS.
//
// The C ABI exposed by sqisign-ffi guarantees that Rust panics never
// cross the boundary; every entry point either returns 1 for success or
// 0 for any failure. The Go wrappers translate that into idiomatic
// (result, error) pairs.
//
// SQIsign is a NIST Round 2 candidate. This library has not been
// audited. See SECURITY.md at the repository root before any
// production use.
package sqisign

// #cgo LDFLAGS: ${SRCDIR}/../../../target/release/libsqisign_ffi.a -lm -ldl -lpthread
// #include "../../../crates/sqisign-ffi/include/sqisign.h"
import "C"

import (
	"errors"
	"unsafe"
)

// Wire sizes of the SQIsign level-1 byte buffers. These match the C
// ABI constants in sqisign.h byte-for-byte and the
// SQISIGN_LVL1_*_BYTES Rust constants in crates/sqisign-ffi/src/lib.rs.
const (
	PublicKeyBytes = C.SQISIGN_LVL1_PUBLIC_KEY_BYTES
	SecretKeyBytes = C.SQISIGN_LVL1_SECRET_KEY_BYTES
	SignatureBytes = C.SQISIGN_LVL1_SIGNATURE_BYTES
	EntropyBytes   = C.SQISIGN_LVL1_ENTROPY_BYTES
)

// Errors returned by this package. The C ABI maps every failure mode to
// a single 0 return; we surface a small handful of Go errors that
// distinguish "caller passed a wrong-sized buffer" (a programming bug)
// from "the algorithm itself failed" (rare, but possible because both
// keygen and sign internally loop until success in the C reference).
var (
	ErrInvalidEntropy   = errors.New("sqisign: entropy must be exactly EntropyBytes long")
	ErrInvalidPublicKey = errors.New("sqisign: public key must be exactly PublicKeyBytes long")
	ErrInvalidSecretKey = errors.New("sqisign: secret key must be exactly SecretKeyBytes long")
	ErrInvalidSignature = errors.New("sqisign: signature must be exactly SignatureBytes long")
	ErrKeygenFailed     = errors.New("sqisign: keygen returned failure")
	ErrSignFailed       = errors.New("sqisign: sign returned failure")
)

// KeyGen produces a fresh SQIsign level-1 keypair seeded by the given
// 48-byte entropy block. The returned public key has length
// PublicKeyBytes; the secret key has length SecretKeyBytes.
//
// The entropy block is consumed by a KAT-compatible NIST CTR-DRBG
// inside the FFI. Production callers who want to avoid that DRBG
// should call the Rust-level sqisign_sign::protocols_keygen entry
// point directly with their own RngSource; this Go binding is for the
// case where the caller provides hardware-derived entropy and is happy
// with the KAT-compatible DRBG.
func KeyGen(entropy []byte) (publicKey, secretKey []byte, err error) {
	if len(entropy) != EntropyBytes {
		return nil, nil, ErrInvalidEntropy
	}
	pk := make([]byte, PublicKeyBytes)
	sk := make([]byte, SecretKeyBytes)
	rc := C.sqisign_lvl1_keygen(
		(*C.uchar)(unsafe.Pointer(&pk[0])), C.size_t(len(pk)),
		(*C.uchar)(unsafe.Pointer(&sk[0])), C.size_t(len(sk)),
		(*C.uchar)(unsafe.Pointer(&entropy[0])), C.size_t(len(entropy)),
	)
	if rc != 1 {
		return nil, nil, ErrKeygenFailed
	}
	return pk, sk, nil
}

// Sign produces a SQIsign signature over msg using secretKey, seeded by
// the given 48-byte entropy block. The returned signature has length
// SignatureBytes. msg may be nil if its length is also 0.
//
// The NIST "sm = signature || msg" concatenation is the caller's
// responsibility; this entry point returns the signature alone.
func Sign(secretKey, msg, entropy []byte) ([]byte, error) {
	if len(secretKey) != SecretKeyBytes {
		return nil, ErrInvalidSecretKey
	}
	if len(entropy) != EntropyBytes {
		return nil, ErrInvalidEntropy
	}
	sig := make([]byte, SignatureBytes)
	var msgPtr *C.uchar
	if len(msg) != 0 {
		msgPtr = (*C.uchar)(unsafe.Pointer(&msg[0]))
	}
	rc := C.sqisign_lvl1_sign(
		(*C.uchar)(unsafe.Pointer(&sig[0])), C.size_t(len(sig)),
		msgPtr, C.size_t(len(msg)),
		(*C.uchar)(unsafe.Pointer(&secretKey[0])), C.size_t(len(secretKey)),
		(*C.uchar)(unsafe.Pointer(&entropy[0])), C.size_t(len(entropy)),
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
