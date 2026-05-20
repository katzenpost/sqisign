/*
 * c_demo.c: self-contained C program that verifies one level-1 SQIsign
 *           KAT entry (count = 0 in PQCsignKAT_353_SQIsign_lvl1.rsp) via
 *           the sqisign-ffi C ABI. Useful as a hand-compile smoke check
 *           of the staticlib + header pair.
 *
 * Build (from the workspace root):
 *
 *     cargo build -p sqisign-ffi --release
 *
 *     # Link against the staticlib (no runtime ld dependency):
 *     cc -O2 -Wall -Wextra \
 *        -I crates/sqisign-ffi/include \
 *        crates/sqisign-ffi/examples/c_demo.c \
 *        target/release/libsqisign_ffi.a \
 *        -lpthread -ldl -lm \
 *        -o /tmp/sqisign_c_demo
 *
 *     # Or link against the cdylib:
 *     cc -O2 -Wall -Wextra \
 *        -I crates/sqisign-ffi/include \
 *        crates/sqisign-ffi/examples/c_demo.c \
 *        -L target/release -lsqisign_ffi \
 *        -o /tmp/sqisign_c_demo_dyn
 *     LD_LIBRARY_PATH=target/release /tmp/sqisign_c_demo_dyn
 *
 * Expected output:
 *
 *     sqisign_lvl1_verify -> 1 (valid)
 *
 * This file is not built by cargo; cargo's `examples/` directory under a
 * cdylib crate is treated as a *Rust* example. The KAT vector embedded
 * below is the first record of PQCsignKAT_353_SQIsign_lvl1.rsp; the
 * canonical KAT file lives at kat/PQCsignKAT_353_SQIsign_lvl1.rsp at
 * the repository root.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sqisign.h"

/* Hex decoder. Strict: lowercase or uppercase; returns 0 on success. */
static int hex_decode(const char *hex, unsigned char *out, size_t out_len) {
    size_t expected = out_len * 2;
    if (strlen(hex) != expected) {
        return -1;
    }
    for (size_t i = 0; i < out_len; i++) {
        unsigned int byte;
        if (sscanf(hex + (i * 2), "%2x", &byte) != 1) {
            return -1;
        }
        out[i] = (unsigned char)byte;
    }
    return 0;
}

/* KAT count=0, lifted from PQCsignKAT_353_SQIsign_lvl1.rsp. */
static const char PK_HEX[] =
    "07CCD21425136F6E865E497D2D4D208F0054AD81372066E817480787AAF7B202"
    "9550C89E892D618CE3230F23510BFBE68FCCDDAEA51DB1436B462ADFAF008A01"
    "0B";

/* sm = signature (148 bytes) || message (33 bytes). */
static const char SM_HEX[] =
    "84228651F271B0F39F2F19F2E8718F31ED3365AC9E5CB303AFE663D0CFC11F04"
    "55D891B0CA6C7E653F9BA2667730BB77BEFE1B1A31828404284AF8FD7BAACC01"
    "0001D974B5CA671FF65708D8B462A5A84A1443EE9B5FED7218767C9D85CEED04"
    "DB0A69A2F6EC3BE835B3B2624B9A0DF68837AD00BCACC27D1EC806A44840267471"
    "D86EFF3447018ADB0A6551EE8322AB30010202"
    "D81C4D8D734FCBFBEADE3D3F8A039FAA2A2C9957E835AD55B22E75BF57BB556AC8";

#define MSG_LEN 33
#define SM_LEN  (SQISIGN_LVL1_SIGNATURE_BYTES + MSG_LEN)

int main(void) {
    unsigned char pk[SQISIGN_LVL1_PUBLIC_KEY_BYTES];
    unsigned char sm[SM_LEN];

    if (hex_decode(PK_HEX, pk, sizeof pk) != 0) {
        fprintf(stderr, "pk hex decode failed\n");
        return EXIT_FAILURE;
    }
    if (hex_decode(SM_HEX, sm, sizeof sm) != 0) {
        fprintf(stderr, "sm hex decode failed\n");
        return EXIT_FAILURE;
    }

    const unsigned char *sig = sm;
    const unsigned char *msg = sm + SQISIGN_LVL1_SIGNATURE_BYTES;

    int r = sqisign_lvl1_verify(
        sig, SQISIGN_LVL1_SIGNATURE_BYTES,
        pk,  SQISIGN_LVL1_PUBLIC_KEY_BYTES,
        msg, MSG_LEN);

    printf("sqisign_lvl1_verify -> %d (%s)\n",
           r, r == 1 ? "valid" : "invalid");

    return (r == 1) ? EXIT_SUCCESS : EXIT_FAILURE;
}
