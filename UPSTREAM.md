# Pinned upstream reference

The C reference at <https://github.com/SQISign/the-sqisign> is the
canonical oracle for this port. It is **not vendored** in this
repository.

Specification: SQIsign Round 2.

| Field            | Value                                      |
|------------------|--------------------------------------------|
| Repository       | https://github.com/SQISign/the-sqisign     |
| Pinned commit    | `dd133d7aca576c361a270c8e6434832535b42ecc` |
| Pinned at (date) | 2026-05-17                                 |
| Subject          | make RNG state thread_local                |

## What this pin commits us to

The lvl1 KAT response file at `kat/PQCsignKAT_353_SQIsign_lvl1.rsp` is
a byte copy of the upstream's `KAT/PQCsignKAT_353_SQIsign_lvl1.rsp`
at the pinned commit. The end-to-end KAT round-trip test in
`crates/sqisign-sign/tests/kat_sign.rs` (and the verify-only counterpart
in `crates/sqisign-verify/tests/kat_verify.rs`) replays each recorded
entry against the Rust port and bit-compares the produced public key,
secret key, and signed message against the upstream's.

The `vectors/precomp/` directory carries canonical-bytes JSON copies
of the lvl1 precomputed cryptographic constants (`EXTREMAL_ORDERS`,
`TORSION_PLUS_2POWER`, etc.) at the same pinned commit. These are
embedded into the `sqisign-precomp` crate at compile time via
`include_str!`; they are runtime data, not test fixtures.

## Reproducing the audit from the pin

```sh
git clone https://github.com/SQISign/the-sqisign /tmp/the-sqisign
git -C /tmp/the-sqisign checkout dd133d7aca576c361a270c8e6434832535b42ecc

# Verify our carried KAT matches the upstream's
diff /tmp/the-sqisign/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp \
     kat/PQCsignKAT_353_SQIsign_lvl1.rsp     # must be empty

# Run the round-trip
cargo test -p sqisign-sign --release -- --ignored kat_lvl1_full
```

## Updating the pin

Bumping this commit is never automatic. It requires:

1. Updating the commit and date in this file, in the same commit that
   copies the new `kat/PQCsignKAT_353_SQIsign_lvl1.rsp` from the
   upstream tree.
2. Re-regenerating the `vectors/precomp/` JSON files if the upstream
   constants moved, and reviewing the diff.
3. A human reviewer on the sync PR.
