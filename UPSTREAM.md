# Pinned upstream reference

The C reference at <https://github.com/SQISign/the-sqisign> is the canonical
oracle for this port. It is vendored as a git submodule at
`vendor/the-sqisign`, pinned to the commit below. Every committed test vector
records this same commit in its `upstream_commit` field; CI fails if the
submodule, this file, and the vectors disagree.

Specification: SQIsign Round 2.

| Field            | Value                                      |
|------------------|--------------------------------------------|
| Repository       | https://github.com/SQISign/the-sqisign     |
| Pinned commit    | `dd133d7aca576c361a270c8e6434832535b42ecc` |
| Pinned at (date) | 2026-05-17                                 |
| Subject          | make RNG state thread_local                |

## Updating the pin

Bumping this commit is never automatic. It requires:

1. Advancing the submodule and this file in the same commit.
2. Regenerating every vector file under `vectors/`.
3. Reviewing the vector diff as part of the upstream-sync PR (a non-empty
   diff at a stable boundary is itself a signal worth understanding).
4. A human reviewer on the sync PR. See "Escalation Triggers" in the plan.
