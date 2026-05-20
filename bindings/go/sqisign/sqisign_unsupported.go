// SPDX-FileCopyrightText: (c) 2026 David Stainton
// SPDX-License-Identifier: GPL-3.0-or-later

//go:build !linux || !amd64

// Pure-Go stub for platforms with no vendored sqisign-ffi staticlib.
//
// The cgo implementation in sqisign.go is gated on linux/amd64 because
// that is the only platform for which a prebuilt staticlib is currently
// shipped under lib/. On every other GOOS/GOARCH this file compiles in
// its place, satisfies the public API surface, and reports
// ErrUnsupported at runtime. Constants stay valid everywhere so callers
// can size buffers and compose the binding into hybrid schemes without
// platform-specific imports.

package sqisign

import (
	"errors"
	"io"
)

// Wire sizes of the SQIsign level-1 byte buffers, mirroring the values
// the cgo implementation gets from the C ABI. Kept in sync with the
// SQISIGN_LVL1_*_BYTES constants in crates/sqisign-ffi/include/sqisign.h.
const (
	PublicKeyBytes = 65
	SecretKeyBytes = 353
	SignatureBytes = 148
)

// The non-Unsupported errors are declared here too so downstream code
// that checks for them under errors.Is compiles on every platform
// without conditional imports. The stub itself never returns them.
var (
	ErrNilRNG           = errors.New("sqisign: rng must not be nil")
	ErrInvalidPublicKey = errors.New("sqisign: public key must be exactly PublicKeyBytes long")
	ErrInvalidSecretKey = errors.New("sqisign: secret key must be exactly SecretKeyBytes long")
	ErrInvalidSignature = errors.New("sqisign: signature must be exactly SignatureBytes long")
	ErrKeygenFailed     = errors.New("sqisign: keygen returned failure")
	ErrSignFailed       = errors.New("sqisign: sign returned failure")
)

// KeyGen always returns ErrUnsupported on platforms without a vendored
// staticlib.
func KeyGen(_ io.Reader) (publicKey, secretKey []byte, err error) {
	return nil, nil, ErrUnsupported
}

// Sign always returns ErrUnsupported on platforms without a vendored
// staticlib.
func Sign(_ io.Reader, _, _ []byte) ([]byte, error) {
	return nil, ErrUnsupported
}

// Verify always returns (false, ErrUnsupported) on platforms without a
// vendored staticlib.
func Verify(_, _, _ []byte) (bool, error) {
	return false, ErrUnsupported
}
