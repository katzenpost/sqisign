package sqisign

import (
	"bytes"
	"testing"
)

// fixedEntropy returns a deterministic 48-byte block so tests are
// reproducible across runs.
func fixedEntropy(seed byte) []byte {
	out := make([]byte, EntropyBytes)
	for i := range out {
		out[i] = seed ^ byte(i)
	}
	return out
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
	if EntropyBytes != 48 {
		t.Errorf("EntropyBytes = %d, want 48", EntropyBytes)
	}
}

func TestRoundTrip(t *testing.T) {
	entropy := fixedEntropy(0x5a)
	pk, sk, err := KeyGen(entropy)
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
	sig, err := Sign(sk, msg, fixedEntropy(0xa5))
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

func TestVerifyRejectsTamperedSignature(t *testing.T) {
	entropy := fixedEntropy(0x42)
	pk, sk, err := KeyGen(entropy)
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	msg := []byte("important payload")
	sig, err := Sign(sk, msg, fixedEntropy(0x24))
	if err != nil {
		t.Fatalf("Sign: %v", err)
	}
	// Flip a bit deep inside the signature.
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
	entropy := fixedEntropy(0x11)
	pk, sk, err := KeyGen(entropy)
	if err != nil {
		t.Fatalf("KeyGen: %v", err)
	}
	msg := []byte("a payload to mangle")
	sig, err := Sign(sk, msg, fixedEntropy(0x22))
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
	_, _, err := KeyGen(make([]byte, EntropyBytes-1))
	if err != ErrInvalidEntropy {
		t.Errorf("KeyGen short entropy: got %v, want ErrInvalidEntropy", err)
	}
	_, err = Sign(make([]byte, SecretKeyBytes-1), nil, fixedEntropy(0))
	if err != ErrInvalidSecretKey {
		t.Errorf("Sign short secret key: got %v, want ErrInvalidSecretKey", err)
	}
	_, err = Verify(make([]byte, SignatureBytes-1), make([]byte, PublicKeyBytes), nil) //nolint:errcheck
	if err != ErrInvalidSignature {
		t.Errorf("Verify short signature: got %v, want ErrInvalidSignature", err)
	}
}
