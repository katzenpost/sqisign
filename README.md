# sqisign-rs

A pure-Rust port of [SQIsign](https://sqisign.org/) Round 2, level 1.

SQIsign is a post-quantum digital signature scheme based on isogenies of
supersingular elliptic curves and the Deuring correspondence between
ideals and isogenies. It is one of the most compact post-quantum
signatures (the level-1 public key is 65 bytes and the signature is 148
bytes), at the cost of comparatively expensive key generation and
signing. This crate is a from-scratch transcription of the
[C reference implementation](https://github.com/SQISign/the-sqisign)
into safe Rust, intended as an audit-friendly alternative for projects
that want to ship SQIsign without taking on a C dependency.

The full C reference is vendored as a git submodule at
`vendor/the-sqisign`; every boundary in the port is witnessed against
the same upstream commit (see `UPSTREAM.md`) by differential test
vectors plus the NIST KAT round-trip.

## Status

| Level   | Keygen | Sign | Verify | KAT |
|---------|:------:|:----:|:------:|:---:|
| lvl1    |   ✓    |  ✓   |   ✓    |  ✓  |
| lvl3    |   -    |  -   |   -    |  -  |
| lvl5    |   -    |  -   |   -    |  -  |

Only level 1 (NIST category 1, 128-bit classical / 64-bit quantum
security) is ported. The lvl3 and lvl5 parameter sets share the same
algorithm structure but use different precomputed constants and field
moduli; the work to port them is well-defined but has not been done.

All workspace tests pass under `cargo test --workspace`, including the
KAT round-trip witness against the first vector of the upstream-recorded
`PQCsignKAT_353_SQIsign_lvl1.rsp`. The full 100-vector battery is
gated behind `#[ignore]` for speed and is invoked explicitly:

```sh
cargo test -p sqisign-sign --release -- --ignored kat_lvl1_full
```

## Performance

Measured on AMD EPYC Turin at 2.4 GHz, release build, single-threaded,
no AVX assembly. The right column is the upstream C reference (`ref`
build, mini-gmp, same host, also no AVX) for comparison.

| Operation | sqisign-rs   | upstream C `ref` | ratio (rs / C) |
|-----------|-------------:|-----------------:|---------------:|
| verify    |    3.64 ms   |          3.20 ms |        1.14×   |
| keygen    |   36.68 ms   |         28.40 ms |        1.29×   |
| sign      |  179.67 ms   |         74.38 ms |        2.42×   |

The verify ratio is small because the verifier spends most of its time
in field arithmetic, which is a portable radix-2^51 limb representation
in safe Rust on both sides. The sign ratio is the largest because the
signing path runs the lattice reduction inner loop (`quat_lll_core`)
many times, which is where GMP's hand-tuned assembler beats
`num-bigint` most cleanly. Upstream's `opt` and `broadwell` builds
add AVX2-tuned p5248 arithmetic on top, widening the gap further.

If a smaller gap matters for your deployment, see
`SPIKE-MALACHITE.md` for a recorded microbench of one alternative
backend (`malachite`) plus a discussion of the principled paths to
parity (`rug` / `gmp-mpfr-sys`, or a bespoke fixed-width `Ibz`).

To reproduce the numbers above:

```sh
cargo bench -p sqisign-verify --bench verify
cargo bench -p sqisign-sign   --bench keygen
cargo bench -p sqisign-sign   --bench sign
```

## Building

The port is plain Rust with no system dependencies, no `unsafe`
outside the FFI crate, and no build scripts. A stock Rust toolchain
suffices:

```sh
cargo build --workspace --release
cargo test  --workspace
```

The C reference is only needed if you want to regenerate test vectors
(see "Verifying correctness" below). It is provided as a git submodule;
clone with `--recurse-submodules` or run
`git submodule update --init --recursive` after cloning.

## Using the library

### From Rust

The high-level API mirrors the NIST SQIsign API shape but expects the
caller to supply an `RngSource` implementation rather than relying on a
process-wide DRBG. The `sqisign-common` crate ships a KAT-compatible
NIST CTR-DRBG; production callers are expected to plug in their own
RNG (a hashed system-entropy source, a hardware RNG, or any source
that satisfies the `RngSource` trait).

```rust
use sqisign_common::CtrDrbg;
use sqisign_sign::{protocols_keygen, protocols_sign, SecretKey};
use sqisign_verify::{protocols_verify, PublicKey, Signature};

// 1. Generate a keypair.
let entropy = [0u8; 48];   // replace with a real entropy source
let mut rng = CtrDrbg::new(&entropy, None);
let mut pk = PublicKey::zero();
let mut sk = SecretKey::new();
assert_eq!(protocols_keygen(&mut rng, &mut pk, &mut sk), 1);

// 2. Sign a message.
let mut sig = Signature::zero();
let message = b"hello sqisign";
assert_eq!(protocols_sign(&mut rng, &mut sig, &pk, &mut sk, message), 1);

// 3. Verify.
assert!(protocols_verify(&sig, &pk, message));
```

The byte-level serializers (`public_key_to_bytes`, `secret_key_to_bytes`,
`signature_to_bytes`, and their inverses) match the wire formats of the
NIST KAT submission. See the test in `crates/sqisign-sign/tests/kat_sign.rs`
for a complete keygen+sign+verify round-trip against the recorded KAT.

### From C

The `sqisign-ffi` crate produces both a `cdylib` and a `staticlib` with
a small C ABI:

```c
#include <stdint.h>

int sqisign_lvl1_keygen(uint8_t *pk, size_t pk_len,
                        uint8_t *sk, size_t sk_len,
                        const uint8_t *entropy, size_t entropy_len);

int sqisign_lvl1_sign(uint8_t *sig, size_t sig_len,
                      const uint8_t *msg, size_t msg_len,
                      const uint8_t *sk, size_t sk_len,
                      const uint8_t *entropy, size_t entropy_len);

int sqisign_lvl1_verify(const uint8_t *sig, size_t sig_len,
                        const uint8_t *pk, size_t pk_len,
                        const uint8_t *msg, size_t msg_len);
```

Buffer sizes are fixed: 65 bytes for `pk`, 353 bytes for `sk`, 148 bytes
for `sig`, 48 bytes for `entropy`. Each function returns `1` on
success, `0` on any failure (length mismatch, null pointer, panic
caught at the FFI boundary, or the algorithm returning a non-success
status). The Rust panics that should never reach a production caller
are turned into a `0` return via `catch_unwind` at every FFI boundary.

A complete C example is at `crates/sqisign-ffi/examples/`.

## Architecture

The port mirrors the C reference's module structure:

```
crates/
  sqisign-common/      common helpers: SHAKE, SHA3, CTR-DRBG, RngSource trait
  sqisign-mp/          multi-precision integer helpers (fixed-width digits)
  sqisign-gf/          finite-field arithmetic: GF(p), GF(p^2), radix-2^51
  sqisign-ec/          Montgomery curve arithmetic, isogeny evaluation, pairings
  sqisign-hd/          (2,2)-isogenies on abelian surfaces (theta model)
  sqisign-precomp/     precomputed constants for the lvl1 parameter set
  sqisign-quaternion/  quaternion orders, ideals, LLL, normeq, lat-ball sampler
  sqisign-id2iso/      ideal-to-isogeny translation (Clapotis style)
  sqisign-verify/      protocols_verify, signature/public-key (de)serialization
  sqisign-sign/        protocols_keygen, protocols_sign, secret-key serialization
  sqisign-ffi/         C ABI cdylib + staticlib (verify + sign + keygen)
  sqisign/             top-level meta-crate
  sqisign-vectors/     differential-vector loader (used only by tests)

vendor/
  the-sqisign/         pinned C reference; the canonical oracle

tools/
  cdump/               C program that dumps deterministic test vectors
  vector-gen/          Rust program that converts cdump output to JSON
  gen-vectors.sh       regenerate every vector from the C reference

vectors/               committed reference vectors, organised by module
```

The `sqisign-verify` crate is deliberately small: a Katzenpost mix node
embedding it pulls in `sqisign-common`, `sqisign-mp`, `sqisign-gf`,
`sqisign-ec`, `sqisign-hd`, `sqisign-precomp`, `sqisign-id2iso` (the
deterministic subset only), and `sqisign-verify` itself. The signing
path (`sqisign-sign`) adds the LLL, normeq, and dim-2 paths inside
`sqisign-quaternion` and `sqisign-id2iso`; it is never pulled into a
verify-only deployment.

## Verifying correctness

This is the most important section of this document. SQIsign is a
recent cryptographic scheme and a port like this carries a substantial
risk of subtle correctness errors that may not show up in normal use
but make signatures forgeable or leak information. The port is not a
substitute for a security audit (see the warning at the end), but the
witnesses below are what a careful human reviewer can use to
independently convince themselves that the port computes the same
function as the upstream C reference.

### The upstream pin

The C reference at <https://github.com/SQISign/the-sqisign> is pinned
to a single commit in `UPSTREAM.md`. Every test vector and the KAT
witness records that same commit hash in its `upstream_commit` field;
CI fails if `UPSTREAM.md`, the vendored submodule, and the vectors
ever disagree on the hash. To inspect the pin:

```sh
cat UPSTREAM.md
git -C vendor/the-sqisign rev-parse HEAD
```

### Differential test vectors

Every primitive (~262 boundaries at the time of writing) is witnessed
by a JSON file of input/output pairs that were dumped from the pinned
C reference. The convention is:

- The C harness in `tools/cdump/` exposes each upstream function as
  a deterministic battery: a fixed-seed driver builds inputs, calls
  the upstream function, and emits a binary record containing
  canonical bytes of the inputs and outputs.
- The Rust program `tools/vector-gen/` converts that binary into a
  JSON file under `vectors/<module>/<boundary>.json`. Each JSON
  records the upstream commit hash, a `generated_at` timestamp (fixed
  to the pin date, not wall-clock), and a list of vectors with
  canonical bytes for every named field.
- Each Rust crate has a per-boundary test file (e.g.
  `crates/sqisign-quaternion/tests/quat_lll_core_vectors.rs`) that
  loads the JSON, replays the inputs through the Rust port, and
  asserts byte-equality with the recorded output.

To reproduce the vectors from a clean checkout:

```sh
git submodule update --init --recursive
./tools/gen-vectors.sh
git status vectors/                  # should be clean
```

The `./tools/gen-vectors.sh` script regenerates every committed JSON
file. Because the harness uses fixed seeds and `generated_at` is taken
from `UPSTREAM.md` rather than the current clock, a clean rerun is
byte-identical to the committed JSON. A non-empty `git status
vectors/` after running it means the pinned commit, the harness, or
the vector format drifted; CI runs the same comparison.

To run the full per-boundary differential suite:

```sh
cargo test --workspace                  # exercises every committed vector
```

### KAT round-trip

The strongest end-to-end witness is the NIST KAT response file at
`vendor/the-sqisign/KAT/PQCsignKAT_353_SQIsign_lvl1.rsp`. It records,
for each of 100 fixed-seed tests, the input entropy and the resulting
public key, secret key, and signed message that the upstream
implementation produces.

`crates/sqisign-sign/tests/kat_sign.rs` parses that file, seeds the
Rust CTR-DRBG with the recorded entropy, runs `protocols_keygen`
followed by `protocols_sign`, and bit-compares the serialised public
key, secret key, and signed message against the recorded bytes for
every entry. A single divergence anywhere in the keygen or sign chain
is detected at the byte level.

- The single-vector witness `kat_lvl1_count_0` runs by default in
  `cargo test --workspace` (~120 ms in release).
- The full 100-vector battery is gated behind `#[ignore]` for default
  test speed; invoke explicitly:

  ```sh
  cargo test -p sqisign-sign --release -- --ignored kat_lvl1_full
  ```

  Expect roughly 12 seconds on a modest desktop. Every vector must
  match byte-for-byte.

### A review checklist a human can follow

1. **Pin parity.** `cat UPSTREAM.md` and verify the date and commit are
   the ones quoted in any paper, advisory, or release note you trust.
   Then `git -C vendor/the-sqisign rev-parse HEAD` and confirm it
   matches.
2. **Vector reproducibility.** Run `./tools/gen-vectors.sh` from a
   clean checkout; `git status vectors/` must be clean afterwards. A
   diff there means the C harness, the vendored upstream, or the
   vector format silently changed.
3. **Per-boundary parity.** Run `cargo test --workspace`. Every
   committed vector must pass; a failure means the corresponding Rust
   port has diverged from the C reference at that primitive.
4. **End-to-end parity.** Run `cargo test -p sqisign-sign --release
   -- --ignored kat_lvl1_full`. Every one of the 100 KAT round-trips
   must match byte-for-byte.
5. **Spot audit.** Pick one or two boundaries (the LLL reducer in
   `crates/sqisign-quaternion/src/lll.rs` and the rejection sampler
   in `crates/sqisign-quaternion/src/lat_ball.rs` are good choices)
   and read the Rust source side-by-side with
   `vendor/the-sqisign/src/quaternion/ref/generic/lll/l2.c` and
   `vendor/the-sqisign/src/quaternion/ref/generic/lat_ball.c`. The
   port mirrors the C reference's control flow line-for-line where
   feasible, which is what makes that kind of diff readable.

None of these by themselves rules out a subtle bug; together they
narrow the search.

## Provenance

Some directories under `vendor/the-sqisign` are governed by the
upstream's own license (Apache 2.0, with portions under MIT or BSD
from the underlying NIST submission packaging). This crate as a whole
is licensed under the GNU General Public License v3 or later; see
`LICENSE` for the full text.

## Acknowledgements

The mathematics is the work of the SQIsign team (see
<https://sqisign.org>). The C reference implementation that this port
follows line-for-line is at <https://github.com/SQISign/the-sqisign>,
under Apache 2.0. Any correctness this port has, it owes to that
reference. Any bug it has is the port's own.

---

> **⚠️ This library has not been audited.**
>
> sqisign-rs is a from-scratch transcription of a research-stage
> post-quantum signature scheme by a single author (with model
> assistance). It is byte-compatible with the published C reference
> on the lvl1 KAT, which is a strong correctness witness, but
> byte-compatibility is not the same as security. The port has not
> been reviewed by a cryptography auditor; it has not been hardened
> against side channels; it relies on `num-bigint`, which is not
> constant-time; and SQIsign itself is a Round 2 NIST candidate, not
> a standardised algorithm. Do not use this crate in production
> systems where a forgery or key-extraction would cause harm. If you
> need a deployable post-quantum signature today, use one of the
> NIST-standardised algorithms (ML-DSA or SLH-DSA) from a vetted
> implementation.
