# Pinned upstream reference

The C reference at <https://github.com/SQISign/the-sqisign> is the
canonical oracle for this port. It is **not vendored** in this repository.
Every committed test vector records the pinned commit below in its
`upstream_commit` field; the CI vector-regeneration workflow clones the
upstream at this exact commit, regenerates every vector, and fails if
the regenerated set differs from the committed one.

Specification: SQIsign Round 2.

| Field            | Value                                      |
|------------------|--------------------------------------------|
| Repository       | https://github.com/SQISign/the-sqisign     |
| Pinned commit    | `dd133d7aca576c361a270c8e6434832535b42ecc` |
| Pinned at (date) | 2026-05-17                                 |
| Subject          | make RNG state thread_local                |

## Reproducing the audit from the pin

```sh
git clone https://github.com/SQISign/the-sqisign /tmp/the-sqisign
git -C /tmp/the-sqisign checkout dd133d7aca576c361a270c8e6434832535b42ecc
UPSTREAM_REF=/tmp/the-sqisign ./tools/gen-vectors.sh
git status vectors/                                  # must be clean
```

The KAT response file used by the round-trip test is carried in this
repository at `kat/PQCsignKAT_353_SQIsign_lvl1.rsp`; it is a byte copy
of the upstream's own `KAT/PQCsignKAT_353_SQIsign_lvl1.rsp` at the
pinned commit. Its license is the upstream's (Apache 2.0).

## Updating the pin

Bumping this commit is never automatic. It requires:

1. Updating the commit and date in this file in the same commit that
   updates the `kat/` file (if the upstream's KAT changed) and the
   vectors (which the script regenerates).
2. Reviewing the vector and KAT diffs as part of the upstream-sync PR; a
   non-empty diff at a stable boundary is itself a signal worth
   understanding.
3. A human reviewer on the sync PR.
