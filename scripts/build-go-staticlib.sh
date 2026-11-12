#!/bin/sh
# Rebuild the Go binding's vendored linux/amd64 staticlib.
#
# musl target: the archive then references POSIX libc names that both
# musl (Alpine) and glibc resolve; see bindings/go/README.md.
set -eu

root=$(cd "$(dirname "$0")/.." && pwd)
target=x86_64-unknown-linux-musl
built=$root/target/$target/release/libsqisign_ffi.a
out=$root/bindings/go/sqisign/lib/linux_amd64

rustup target add "$target"
cargo build --manifest-path "$root/Cargo.toml" --release --target "$target" -p sqisign-ffi

# Rust std's prebuilt rlibs embed LLVM bitcode (.llvmbc/.llvmcmd) whatever
# our build flags say. Go >= 1.26.11 rejects cgo archives carrying it
# ("disallowed LTO section"), so strip it.
objcopy --remove-section=.llvmbc --remove-section=.llvmcmd "$built" "$out/libsqisign_ffi.a"
if readelf -S "$out/libsqisign_ffi.a" | grep -q '\.llvm\(bc\|cmd\)'; then
	echo "error: LLVM bitcode left in $out/libsqisign_ffi.a" >&2
	exit 1
fi

{ rustc -V; cargo -V; objcopy --version | head -n 1; } > "$out/BUILDINFO"
cat "$out/BUILDINFO"
