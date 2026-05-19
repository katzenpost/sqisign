# Differential fuzzing

The plan requires a `cargo-fuzz` target per crate that interprets libFuzzer
bytes as inputs, runs the Rust implementation (and, ultimately, the linked
C reference), and asserts bit-equal behaviour. Any divergence is a P0 bug.

## Status: harnesses ship, runner deferred

`cargo-fuzz` and libFuzzer are not installed in the current environment,
exactly as `valgrind`/`iai-callgrind` are not (see BENCHMARKS.md). This is
recorded, not silently dropped: the fuzz **harnesses are committed and
build**; only the *running* is deferred to a fuzzing host with the
toolchain.

| Crate | Target | Asserts today | Next increment |
|---|---|---|---|
| `sqisign-common` | `shake256` | determinism, XOF prefix stability | link `fips202.c` for byte-equality vs C |
| `sqisign-common` | `shake128` | determinism, XOF prefix stability, distinct from SHAKE256 | link `fips202.c` for byte-equality vs C |
| `sqisign-common` | `shake_inc` | incremental == one-shot for arbitrary absorb/squeeze chunking, both rates | link `fips202.c` for byte-equality vs C |
| `sqisign-common` | `sha3` | determinism, intrinsic lengths, mutual distinctness of the three widths | link `fips202.c` for byte-equality vs C |
| `sqisign-common` | `ctr_drbg` | determinism, draw reproducible from seed, state evolves identically across calls | link `randombytes_ctrdrbg.c` for byte-equality vs C |
| `sqisign-common` | `secure_clear` | chosen prefix zeroed, tail untouched, for arbitrary buffer/split | link `mem.c` for byte-equality vs C |
| `sqisign-mp` | `mp_add` | commutativity, zero identity, single-limb == wrapping_add | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_sub` | (a+b)-b == a, a-a == 0, single-limb == wrapping_sub | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_mul` | nwords>=2 commutes & == u128 low half; nwords==1 == reproduced 2*(a*b) defect | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_mul2` | reproduced partial product c == a*b - (a1*b0)*2^64 | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_mod_2exp` | per-limb low-bit mask, idempotent, over-width no-op | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_shiftl` | x<<1 == x+x, low shift bits clear, arbitrary width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_shiftr` | returned bit == entry parity, top zero-filled, single-limb == native >> | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `multiple_mp_shiftl` | == mp_shiftl in 1..=63, over-width == 0, arbitrary width | link `mp.c` for byte-equality vs C |

## Running (on a host with the toolchain)

```sh
cargo install cargo-fuzz
cd crates/sqisign-common/fuzz
cargo +nightly fuzz run shake256
```

The fuzz crate is its own Cargo workspace by design, so it does not perturb
the main workspace toolchain.

## Triage

Per the plan: shrink the input, identify the offending function, commit a
fixed-input regression test under the crate's `tests/`, then fix. A
divergence the property tests should have caught means the property tests
are themselves wrong (an escalation trigger).
