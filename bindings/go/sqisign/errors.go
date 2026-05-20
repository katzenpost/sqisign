// SPDX-FileCopyrightText: (c) 2026 David Stainton
// SPDX-License-Identifier: GPL-3.0-or-later

// Cross-platform sentinel errors. Declared in this constraint-free
// file so callers can reference them under errors.Is on every GOOS
// and GOARCH, whether or not the cgo implementation is linked in on
// this target.

package sqisign

import "errors"

// ErrUnsupported is returned by KeyGen, Sign, and Verify on any
// GOOS/GOARCH for which no prebuilt sqisign-ffi staticlib is vendored
// in this binding. Adding a new target reduces to building the
// staticlib for that platform and dropping it into
// lib/<GOOS>_<GOARCH>/libsqisign_ffi.a alongside an updated cgo
// LDFLAGS directive in sqisign.go.
var ErrUnsupported = errors.New("sqisign: no vendored libsqisign_ffi for this GOOS/GOARCH")
