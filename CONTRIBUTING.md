# Contributing

## Running CI locally before pushing

The `test` job in `.github/workflows/ci.yml` runs five steps against a
pinned Rust toolchain (currently 1.94). Newer or older toolchains
emit a different lint set; running plain `cargo clippy` against a
host-default toolchain therefore leaves a real risk of pushing a
change that passes locally and fails in CI on a lint your toolchain
does not yet (or no longer) flag.

### One-time setup

Install the pinned toolchain and its components alongside whatever
default you keep:

```sh
rustup toolchain install 1.94 --component rustfmt clippy
```

This adds 1.94 as an additional toolchain. Your default is
untouched; `cargo` by itself continues to use whatever you had
configured. The CI-equivalent invocations all go through `cargo
+1.94 ...`.

### The convenience script

```sh
./scripts/ci-local.sh
```

runs all five `test`-job steps in order against the pinned
toolchain, stopping at the first failure. Individual steps:

```sh
./scripts/ci-local.sh fmt
./scripts/ci-local.sh clippy
./scripts/ci-local.sh build
./scripts/ci-local.sh test
./scripts/ci-local.sh verify
```

`verify` is the small `cargo tree` check that asserts the
verify-only umbrella crate does not pull in `sqisign-sign`. The
other four names map directly to the corresponding workflow step.

### When the workflow changes

`scripts/ci-local.sh` mirrors `ci.yml` by hand. If you add, remove,
or reorder a step in the workflow, update the script in the same
commit so a local pre-push run continues to reflect what CI will
actually do.

## Style

See `CLAUDE.md` (in the parent directory of this repo if you have
the Katzenpost monorepo nearby, or the abbreviated notes in the
sqisign repo itself) for the project conventions:

- No em dashes; use comma, colon, semicolon, parentheses, or split
  the sentence.
- No banner-style section comments. If code needs section markers
  to be readable, extract helper functions instead.
- Keep functions and methods short, with logical sections extracted
  into well-named helpers.

## Reporting issues

Security-sensitive reports: see `SECURITY.md`.

Other bugs and feature requests: file an issue against
[katzenpost/sqisign](https://github.com/katzenpost/sqisign).
