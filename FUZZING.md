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
| `sqisign-mp` | `mp_neg` | faithful model (no carry past limb 0), == -a iff a[0]!=0 | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_copy` | identity into b, arbitrary width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_compare` | reflexive, antisymmetric, consistent with top differing limb | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_is_zero` | equals the all-limbs-zero predicate, arbitrary width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_is_one` | equals the canonical-one predicate, arbitrary width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `select_ct` | per-bit blend `(a&!m)|(b&m)`; 0/~0 select a/b | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `swap_ct` | per-bit conditional swap; 0/~0 no-op/swap; double-swap involution | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_inv_2e` | a*b == 1 mod 2^e for odd a, e within width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_shiftl` | x<<1 == x+x, low shift bits clear, arbitrary width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_shiftr` | returned bit == entry parity, top zero-filled, single-limb == native >> | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `multiple_mp_shiftl` | == mp_shiftl in 1..=63, over-width == 0, arbitrary width | link `mp.c` for byte-equality vs C |
| `sqisign-mp` | `mp_invert_matrix` | odd-det input never panics, main diagonal of M*Minv == 1 mod 2^e (off-diagonal inherits mp_neg defect) | link `mp.c` for byte-equality vs C |
| `sqisign-gf` | `fp_add` | bit-exact commutativity, structural carry invariant (limbs 0..=3 < 2^51; limb 4 unmasked by design) | link `fp_p5248_64.c` for byte-equality vs C |
| `sqisign-gf` | `fp_sub` | `fp_sub(a,a)` is the canonical all-zero limb vector, structural carry invariant (limbs 0..=3 < 2^51; limb 4 unmasked by design) | link `fp_p5248_64.c` for byte-equality vs C |
| `sqisign-gf` | `fp_neg` | `fp_neg(0)` is the canonical all-zero limb vector, structural carry invariant (limbs 0..=3 < 2^51; limb 4 unmasked by design) | link `fp_p5248_64.c` for byte-equality vs C |
| `sqisign-gf` | `fp_copy` | bit-exact identity into `out`, for arbitrary five-limb inputs | link `fp_p5248_64.c` for byte-equality vs C |
| `sqisign-gf` | `fp_set_zero` | bit-exact canonical all-zero limb vector into `out`, for arbitrary five-limb destination pre-fill | link `fp_p5248_64.c` for byte-equality vs C |
| `sqisign-gf` | `fp_set_one` | bit-exact positional-one limb vector `[1, 0, 0, 0, 0]` into `out`, for arbitrary five-limb destination pre-fill | link `fp_p5248_64.c` for byte-equality vs C |
| `sqisign-gf` | `fp_select` | per-limb bit blend `d[i] = a0[i] ^ (cw & (a0[i] ^ a1[i]))` at the two declared `ctl` endpoints, with `cw` the sign-extended widening of `ctl`; `ctl == 0` selects `a0`, `ctl == 0xFFFFFFFF` selects `a1` | link `fp.c` (lvlx) for byte-equality vs C |

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
