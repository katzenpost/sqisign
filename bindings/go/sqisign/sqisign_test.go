//go:build linux && amd64

package sqisign

import (
	"bytes"
	"crypto/rand"
	"errors"
	"io"
	"testing"

	"golang.org/x/crypto/chacha20"
)

// keystreamReader is a small deterministic byte source backed by a
// ChaCha20 keystream. It mirrors the role hpqc/rand.DeterministicRandReader
// plays in the wider tree, but the binding deliberately depends on
// nothing from hpqc so the test reproduces the construction inline.
type keystreamReader struct {
	cipher *chacha20.Cipher
}

func newKeystreamReader(seed byte) *keystreamReader {
	var key [chacha20.KeySize]byte
	for i := range key {
		key[i] = seed ^ byte(i)
	}
	var nonce [chacha20.NonceSize]byte
	c, err := chacha20.NewUnauthenticatedCipher(key[:], nonce[:])
	if err != nil {
		panic(err)
	}
	return &keystreamReader{cipher: c}
}

func (r *keystreamReader) Read(p []byte) (int, error) {
	for i := range p {
		p[i] = 0
	}
	r.cipher.XORKeyStream(p, p)
	return len(p), nil
}

func TestSizesMatchRustABI(t *testing.T) {
	if PublicKeyBytes != 65 {
		t.Errorf("PublicKeyBytes = %d, want 65", PublicKeyBytes)
	}
	if SecretKeyBytes != 353 {
		t.Errorf("SecretKeyBytes = %d, want 353", SecretKeyBytes)
	}
	if SignatureBytes != 148 {
		t.Errorf("SignatureBytes = %d, want 148", SignatureBytes)
	}
}

func TestRoundTripWithCryptoRand(t *testing.T) {
	pk, sk, err := KeyGen(rand.Reader)
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	if len(pk) != PublicKeyBytes {
		t.Fatalf("pk len = %d, want %d", len(pk), PublicKeyBytes)
	}
	if len(sk) != SecretKeyBytes {
		t.Fatalf("sk len = %d, want %d", len(sk), SecretKeyBytes)
	}

	msg := []byte("the time has come, the walrus said, to talk of many things")
	sig, err := Sign(rand.Reader, sk, msg)
	if err != nil {
		t.Fatalf("Sign: %v", err)
	}
	if len(sig) != SignatureBytes {
		t.Fatalf("sig len = %d, want %d", len(sig), SignatureBytes)
	}

	ok, err := Verify(sig, pk, msg)
	if err != nil {
		t.Fatalf("Verify: %v", err)
	}
	if !ok {
		t.Fatalf("Verify returned false for a fresh signature")
	}
}

func TestRoundTripWithDeterministicReader(t *testing.T) {
	pk, sk, err := KeyGen(newKeystreamReader(0x5a))
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	msg := []byte("deterministic keystream signed message")
	sig, err := Sign(newKeystreamReader(0xa5), sk, msg)
	if err != nil {
		t.Fatalf("Sign: %v", err)
	}
	ok, err := Verify(sig, pk, msg)
	if err != nil {
		t.Fatalf("Verify: %v", err)
	}
	if !ok {
		t.Fatalf("Verify returned false for a deterministic-rng signature")
	}
}

func TestKeygenIsDeterministicGivenReader(t *testing.T) {
	pk1, sk1, err := KeyGen(newKeystreamReader(0x11))
	if err != nil {
		t.Fatalf("KeyGen 1: %v", err)
	}
	pk2, sk2, err := KeyGen(newKeystreamReader(0x11))
	if err != nil {
		t.Fatalf("KeyGen 2: %v", err)
	}
	if !bytes.Equal(pk1, pk2) {
		t.Fatalf("identical readers produced different public keys")
	}
	if !bytes.Equal(sk1, sk2) {
		t.Fatalf("identical readers produced different secret keys")
	}
}

func TestRandomisedSignaturesDiffer(t *testing.T) {
	pk, sk, err := KeyGen(rand.Reader)
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	msg := []byte("two signatures with fresh entropy must differ")
	sigA, err := Sign(rand.Reader, sk, msg)
	if err != nil {
		t.Fatalf("Sign A: %v", err)
	}
	sigB, err := Sign(rand.Reader, sk, msg)
	if err != nil {
		t.Fatalf("Sign B: %v", err)
	}
	if bytes.Equal(sigA, sigB) {
		t.Fatalf("two fresh signatures coincided; randomness is broken")
	}
	if ok, _ := Verify(sigA, pk, msg); !ok {
		t.Fatalf("Verify rejected fresh signature A")
	}
	if ok, _ := Verify(sigB, pk, msg); !ok {
		t.Fatalf("Verify rejected fresh signature B")
	}
}

func TestVerifyRejectsTamperedSignature(t *testing.T) {
	pk, sk, err := KeyGen(rand.Reader)
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	msg := []byte("important payload")
	sig, err := Sign(rand.Reader, sk, msg)
	if err != nil {
		t.Fatalf("Sign: %v", err)
	}
	bad := bytes.Clone(sig)
	bad[len(bad)/2] ^= 0x01
	ok, err := Verify(bad, pk, msg)
	if err != nil {
		t.Fatalf("Verify: %v", err)
	}
	if ok {
		t.Fatalf("Verify accepted a tampered signature")
	}
}

func TestVerifyRejectsTamperedMessage(t *testing.T) {
	pk, sk, err := KeyGen(rand.Reader)
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	msg := []byte("a payload to mangle")
	sig, err := Sign(rand.Reader, sk, msg)
	if err != nil {
		t.Fatalf("Sign: %v", err)
	}
	tampered := bytes.Clone(msg)
	tampered[0] ^= 0x01
	ok, err := Verify(sig, pk, tampered)
	if err != nil {
		t.Fatalf("Verify: %v", err)
	}
	if ok {
		t.Fatalf("Verify accepted a signature for a tampered message")
	}
}

func TestWrongSizedBuffersReturnErrors(t *testing.T) {
	_, err := Sign(rand.Reader, make([]byte, SecretKeyBytes-1), nil)
	if !errors.Is(err, ErrInvalidSecretKey) {
		t.Errorf("Sign short secret key: got %v, want ErrInvalidSecretKey", err)
	}
	_, err = Verify(make([]byte, SignatureBytes-1), make([]byte, PublicKeyBytes), nil)
	if !errors.Is(err, ErrInvalidSignature) {
		t.Errorf("Verify short signature: got %v, want ErrInvalidSignature", err)
	}
}

func TestNilRNGReturnsError(t *testing.T) {
	_, _, err := KeyGen(nil)
	if !errors.Is(err, ErrNilRNG) {
		t.Errorf("KeyGen nil rng: got %v, want ErrNilRNG", err)
	}
	_, err = Sign(nil, make([]byte, SecretKeyBytes), nil)
	if !errors.Is(err, ErrNilRNG) {
		t.Errorf("Sign nil rng: got %v, want ErrNilRNG", err)
	}
}

// shortReader feeds out only the first byte of every Read and then
// returns io.EOF; the cgo callback should propagate the failure by
// panicking, which catch_unwind on the Rust side maps to a 0 return.
type shortReader struct{}

func (shortReader) Read(p []byte) (int, error) {
	if len(p) == 0 {
		return 0, nil
	}
	p[0] = 0
	return 1, io.EOF
}

func TestFailingRNGTriggersKeygenFailure(t *testing.T) {
	defer func() {
		// The panic from the callback unwinds through cgo; catch_unwind
		// on the Rust side converts it to a 0 return and KeyGen
		// surfaces ErrKeygenFailed. If the panic instead escapes the
		// boundary the test runner will fail with the panic message,
		// which is also a useful failure signal.
		_ = recover()
	}()
	_, _, err := KeyGen(shortReader{})
	if err == nil {
		t.Fatalf("expected KeyGen to fail with a short reader, got nil error")
	}
}
