#!/usr/bin/env bash
# Reproduce the .github/workflows/ci.yml `test` job locally.
#
# The steps below mirror the workflow line by line, against the same
# toolchain version it pins. Run this before pushing to catch the
# failures CI would otherwise catch on your behalf.
#
# One-time prerequisite:
#
#   rustup toolchain install 1.94 --component rustfmt clippy
#
# Usage:
#
#   ./scripts/ci-local.sh           # all steps, stops at first failure
#   ./scripts/ci-local.sh fmt       # just the fmt step
#   ./scripts/ci-local.sh clippy    # just the clippy step
#   ./scripts/ci-local.sh build     # just the build step
#   ./scripts/ci-local.sh test      # just the test step
#   ./scripts/ci-local.sh verify    # just the verify-only-tree check
#
# Keep this script in step with .github/workflows/ci.yml. If you
# change one, change the other.

set -euo pipefail

# The toolchain pin must match .github/workflows/ci.yml.
TOOLCHAIN="${SQISIGN_CI_TOOLCHAIN:-1.94}"

cd "$(dirname "$0")/.."

step_fmt() {
    echo ">> cargo +${TOOLCHAIN} fmt --all --check"
    cargo "+${TOOLCHAIN}" fmt --all --check
}

step_clippy() {
    echo ">> cargo +${TOOLCHAIN} clippy --workspace --all-targets -- -D warnings"
    cargo "+${TOOLCHAIN}" clippy --workspace --all-targets -- -D warnings
}

step_build() {
    echo ">> cargo +${TOOLCHAIN} build --workspace --all-targets"
    cargo "+${TOOLCHAIN}" build --workspace --all-targets
}

step_test() {
    echo ">> cargo +${TOOLCHAIN} test --workspace --all-targets"
    cargo "+${TOOLCHAIN}" test --workspace --all-targets
}

step_verify() {
    # Katzenpost mix nodes pull only the verification path; building the
    # umbrella crate with default features must not drag in sqisign-sign.
    local tree
    tree="$(mktemp)"
    trap 'rm -f "${tree}"' RETURN
    echo ">> cargo +${TOOLCHAIN} tree -p sqisign --edges normal"
    cargo "+${TOOLCHAIN}" tree -p sqisign --edges normal | tee "${tree}"
    if grep -q '^.*sqisign-sign' "${tree}"; then
        echo "error: sqisign-sign present in default (verify-only) tree" >&2
        return 1
    fi
}

main() {
    local target="${1:-all}"
    case "${target}" in
        fmt)     step_fmt ;;
        clippy)  step_clippy ;;
        build)   step_build ;;
        test)    step_test ;;
        verify)  step_verify ;;
        all)
            step_fmt
            step_clippy
            step_build
            step_test
            step_verify
            ;;
        *)
            echo "usage: $0 [fmt|clippy|build|test|verify|all]" >&2
            exit 2
            ;;
    esac
    echo "ok"
}

main "$@"
