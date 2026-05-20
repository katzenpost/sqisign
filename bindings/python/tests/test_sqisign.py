"""Smoke tests for the Python binding.

Exercises round-trip keygen + sign + verify, plus the two negative
cases (tampered signature, tampered message) and the wrong-sized
buffer guards. Run with ``python -m pytest tests/`` from the
``bindings/python`` directory.
"""

from __future__ import annotations

import pytest

import sqisign


def fixed_entropy(seed: int) -> bytes:
    """Deterministic 48-byte block so the tests are reproducible."""
    return bytes((seed ^ i) & 0xFF for i in range(sqisign.ENTROPY_BYTES))


def test_sizes_match_rust_abi():
    assert sqisign.PUBLIC_KEY_BYTES == 65
    assert sqisign.SECRET_KEY_BYTES == 353
    assert sqisign.SIGNATURE_BYTES == 148
    assert sqisign.ENTROPY_BYTES == 48


def test_roundtrip():
    pk, sk = sqisign.keygen(fixed_entropy(0x5A))
    assert len(pk) == sqisign.PUBLIC_KEY_BYTES
    assert len(sk) == sqisign.SECRET_KEY_BYTES

    msg = b"the time has come, the walrus said, to talk of many things"
    sig = sqisign.sign(sk, msg, fixed_entropy(0xA5))
    assert len(sig) == sqisign.SIGNATURE_BYTES

    assert sqisign.verify(sig, pk, msg) is True


def test_verify_rejects_tampered_signature():
    pk, sk = sqisign.keygen(fixed_entropy(0x42))
    msg = b"important payload"
    sig = sqisign.sign(sk, msg, fixed_entropy(0x24))
    bad = bytearray(sig)
    bad[len(bad) // 2] ^= 0x01
    assert sqisign.verify(bytes(bad), pk, msg) is False


def test_verify_rejects_tampered_message():
    pk, sk = sqisign.keygen(fixed_entropy(0x11))
    msg = b"a payload to mangle"
    sig = sqisign.sign(sk, msg, fixed_entropy(0x22))
    tampered = bytearray(msg)
    tampered[0] ^= 0x01
    assert sqisign.verify(sig, pk, bytes(tampered)) is False


def test_empty_message_round_trip():
    pk, sk = sqisign.keygen(fixed_entropy(0x33))
    sig = sqisign.sign(sk, b"", fixed_entropy(0x44))
    assert sqisign.verify(sig, pk, b"") is True
    # Also accept None as a convenience for empty payloads.
    assert sqisign.verify(sig, pk, None) is True


def test_wrong_sized_entropy_raises():
    with pytest.raises(sqisign.SqisignError):
        sqisign.keygen(b"\x00" * (sqisign.ENTROPY_BYTES - 1))


def test_wrong_sized_secret_key_raises():
    with pytest.raises(sqisign.SqisignError):
        sqisign.sign(b"\x00" * (sqisign.SECRET_KEY_BYTES - 1), b"msg", fixed_entropy(0))


def test_wrong_sized_signature_raises():
    pk, _ = sqisign.keygen(fixed_entropy(0x77))
    with pytest.raises(sqisign.SqisignError):
        sqisign.verify(b"\x00" * (sqisign.SIGNATURE_BYTES - 1), pk, b"msg")


def test_wrong_sized_public_key_raises():
    pk, sk = sqisign.keygen(fixed_entropy(0x88))
    sig = sqisign.sign(sk, b"msg", fixed_entropy(0x99))
    with pytest.raises(sqisign.SqisignError):
        sqisign.verify(sig, b"\x00" * (sqisign.PUBLIC_KEY_BYTES - 1), b"msg")
