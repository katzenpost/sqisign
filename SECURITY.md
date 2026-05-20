# Security policy and threat model

> **This library has not been audited.** It is a from-scratch Rust port
> of a NIST Round 2 candidate signature scheme. It has not been reviewed
> by a cryptography auditor; it has not been hardened against every
> side channel; the scheme it implements (SQIsign) is itself research
> output, not a standardised algorithm. Do not deploy this code in any
> setting where a forgery or a key-extraction would cause harm without
> at least the steps documented at the bottom of this file.

This document is the security backstop for the one-paragraph "not
audited" warning at the end of the README. It enumerates what the
library claims, what it does not claim, the known gaps, and how to
report a vulnerability.

## What this library claims

1. **KAT equivalence with the C reference.** For the level-1 parameter
   set, the byte output of `protocols_keygen` followed by
   `protocols_sign` (and the byte verdict of `protocols_verify`) is
   bit-identical to the upstream C reference at the commit pinned in
   `UPSTREAM.md`, across the full 100-entry NIST KAT response file at
   `kat/PQCsignKAT_353_SQIsign_lvl1.rsp`. This is witnessed by
   `crates/sqisign-sign/tests/kat_sign.rs::kat_lvl1_full`. The
   single-vector smoke `kat_lvl1_count_0` runs by default in
   `cargo test --workspace`.

2. **Memory safety, except at the FFI boundary.** Every crate except
   `sqisign-ffi` is `#![forbid(unsafe_code)]`. `sqisign-ffi` has
   `unsafe` blocks at the C-callable surface (to turn raw pointers
   into Rust slices) and wraps each entry point in `catch_unwind` to
   prevent a Rust panic from becoming undefined behaviour on the C
   side.

## What this library does not claim

The list below is non-exhaustive. A real audit would extend it.

1. **No constant-time guarantees.** The underlying `num-bigint`
   library branches on operand values during multiplication,
   division, and comparison. Any code path that handles
   secret-derived integers (`SecretKey.secret_ideal`,
   `SecretKey.mat_BAcan_to_BA0_two`, `resp_quat` and the lattice
   intersections inside `protocols_sign`, the rejection sampler in
   `quat_lattice_sample_from_ball`) leaks timing information correlated
   with the secret. The upstream C reference has the same property;
   we inherit it. Closing this gap requires switching to a
   constant-time integer backend or a bespoke fixed-width `Ibz`.

2. **No memory-disclosure mitigation for secret material.** The
   `Ibz` newtype wraps `num_bigint::BigInt`, which heap-allocates its
   limbs. When an `Ibz`-bearing value is dropped, Rust deallocates the
   limb storage but does not overwrite it. A memory disclosure
   between deallocation and the next allocation reusing the same
   bytes (core dump, swap-out, debugger, page-table side channel,
   use-after-free) can recover the secret. This affects `SecretKey`,
   every intermediate value during `protocols_sign`, and the
   intermediates during `protocols_keygen`. Closing this gap
   requires switching to `num-bigint-dig` (zeroize-aware) or a
   bespoke fixed-width integer type, and then wiring `ZeroizeOnDrop`
   on the secret-bearing structs.

3. **No randomness-supplier guarantees.** The signing-path primitives
   take `&mut impl RngSource`. The trait says nothing about the
   quality of the byte source; a biased or predictable `RngSource`
   implementation produces forgeable signatures. Production callers
   are responsible for providing a cryptographically-secure RNG
   (`OsRng`, a vetted hardware RNG, or an equivalent). The
   `sqisign_common::CtrDrbg` implementation is KAT-only and **must
   not be used in production**: it is the NIST AES-256-CTR_DRBG, which
   the Katzenpost research team has documented reservations about.

4. **No algorithmic guarantees beyond the scheme itself.** SQIsign is
   a NIST Round 2 candidate, not a standardised signature scheme. The
   research community is still examining its concrete-security
   properties. A faithful port inherits whatever weaknesses are
   discovered in the scheme; the port itself cannot defend against
   them.

5. **No supply-chain guarantees beyond what `cargo audit` catches.**
   We pin a small set of third-party crates (`num-bigint`, `aes`,
   `sha3`, `zeroize`, `proptest`, `criterion`, `serde`, `serde_json`,
   `hex`). The CI workflow `.github/workflows/audit.yml` runs
   `cargo audit` on every push and daily, and the build fails on any
   advisory match. That catches *known* advisories; it does not
   detect novel supply-chain attacks on those crates.

## Threat model the library protects against (with the caveats above)

| Adversary capability                           | Protected? |
|------------------------------------------------|:----------:|
| Adversary submits chosen messages and observes outputs | yes, modulo (1) |
| Adversary submits malformed signatures for verification | yes (verifier rejects, never UB) |
| Adversary submits malformed public keys for verification | yes (decoder rejects, never UB) |
| Adversary observes wall-clock timing of `protocols_sign` | **no**, see (1) |
| Adversary observes wall-clock timing of `protocols_verify` | partially (the lattice arithmetic is on public data; the hash and pairing layers branch on public values only, but `num-bigint` arithmetic on public values still has data-dependent timing) |
| Adversary inspects the process's heap after a signing operation | **no**, see (2) |
| Adversary controls the `RngSource` | **no**, see (3) |
| Adversary breaks SQIsign as a scheme | **no**, see (4) |
| Adversary substitutes a backdoored version of a third-party crate | **no**, see (5) |

## Pre-deployment checklist

Before using this library in any setting beyond research replay of the
KAT, do at least:

1. **Have one SQIsign-fluent reviewer read the Rust source side-by-side
   with the upstream C reference at the pinned commit.** The README's
   "Verifying correctness" section has a five-step checklist; that is
   the qualitative backstop the KAT cannot provide. Two to three
   focused days of reviewer time is the right order of magnitude.

2. **Run the fuzz harnesses overnight.** `crates/sqisign-common/fuzz/`,
   `crates/sqisign-gf/fuzz/`, `crates/sqisign-mp/fuzz/`. Each target
   tests determinism and basic invariants on adversary-controlled
   bytes; cargo-fuzz on a fuzzing host with libFuzzer installed,
   eight to twelve hours per target, triage anything that pops.

3. **Confirm the upstream pin matches the version of SQIsign you are
   committing to.** Compare `UPSTREAM.md` against any paper, advisory,
   or release note you trust. The pin can be updated; see the
   instructions at the bottom of `UPSTREAM.md`.

4. **Decide whether the threat model in the table above matches the
   deployment.** If "Adversary inspects the process's heap" is in
   scope (a multi-tenant host, a memory-disclosure CVE class), this
   library is not ready until the BigInt heap-leak gap is closed.

5. **Confirm the RNG you wire matches the quality bar the scheme
   needs.** A real `OsRng` adapter is the conservative default. A
   hardware RNG is acceptable if it passes the platform's health
   tests. The KAT-only `CtrDrbg` must not be wired into production.

## Reporting a vulnerability

This repository ships on GitHub. If you find a vulnerability, please
**do not file a public issue**. Use GitHub's "Report a vulnerability"
feature on the repository's Security tab; that opens a private
advisory visible only to the maintainers.

If GitHub's private advisory mechanism is not available to you (for
example, you found the vulnerability in a downstream fork), reach the
Katzenpost maintainers via the project's documented contact channels
at <https://katzenpost.network>.

We will acknowledge the report within a reasonable window, work with
you on a disclosure timeline, and credit you in the advisory if you
wish.

## Provenance of this document

This `SECURITY.md` was written by the same author who wrote the port
(David Stainton with model assistance). It is best-effort honest about
what the port does and does not protect against, but it is not itself
an audit deliverable; an audit would produce a longer and more careful
document. The right way to read this file is "the maintainer's
self-assessment, intended to surface gaps a reviewer should check."
