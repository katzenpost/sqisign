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
