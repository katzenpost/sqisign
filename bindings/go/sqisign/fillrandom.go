// SPDX-FileCopyrightText: (c) 2026 David Stainton
// SPDX-License-Identifier: GPL-3.0-or-later

package sqisign

// The Go-exported callback the Rust FFI invokes whenever it needs
// randomness. The C side passes the uintptr_t context we handed it
// via gopointer.Save; we restore it back to the io.Reader the caller
// originally supplied and drain it for `len` bytes.
//
// This is the only Go code on the hot path for randomness, and it
// must satisfy the FillRandom callback contract: fill exactly `len`
// bytes at `out`, or panic if the reader cannot deliver them. The
// reference algorithm has no recovery story for partial reads, so
// panicking out of the cgo callback (which catch_unwind on the Rust
// side maps to a 0 return) is the safest available behaviour.
//
// Lives in its own file so the `//export` directive does not collide
// with the `import "C"` preamble in sqisign.go (cgo forbids both in
// the same file).

import (
	"io"
	"unsafe"

	gopointer "github.com/mattn/go-pointer"
)

// #include <stddef.h>
// #include <stdint.h>
import "C"

//export sqisign_go_fillrandom
func sqisign_go_fillrandom(out *C.uchar, length C.size_t, context C.uintptr_t) {
	if length == 0 {
		return
	}
	rng, ok := gopointer.Restore(unsafe.Pointer(uintptr(context))).(io.Reader)
	if !ok {
		panic("sqisign: callback context did not restore to an io.Reader")
	}
	buf := unsafe.Slice((*byte)(unsafe.Pointer(out)), int(length))
	if _, err := io.ReadFull(rng, buf); err != nil {
		panic("sqisign: rng read failed: " + err.Error())
	}
}
