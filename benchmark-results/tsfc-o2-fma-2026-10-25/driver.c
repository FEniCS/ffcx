/* Kernel-only timings: no Python, assembly, or zeroing within timed regions. */
#define _POSIX_C_SOURCE 200809L
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "inputs.h"

#if TSFC
#if ACTION
extern void KERNEL(double *, const double *, const double *);
typedef void (*kernel_t)(double *, const double *, const double *);
#define CALL fn(A, coords, coefficients)
#else
extern void KERNEL(double *, const double *);
typedef void (*kernel_t)(double *, const double *);
#define CALL fn(A, coords)
#endif
#else
extern void KERNEL(double *, const double *, const double *, const double *, const int *, const uint8_t *, void *);
typedef void (*kernel_t)(double *, const double *, const double *, const double *, const int *, const uint8_t *, void *);
#define CALL fn(A, coefficients, constants, coords, NULL, NULL, NULL)
#endif

static double now(void) {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec + 1e-9 * t.tv_nsec;
}

int main(int argc, char **argv) {
    double *A = calloc(TENSOR_SIZE, sizeof(double));
    if (!A) return 2;
    volatile kernel_t fn = KERNEL;
    if (argc == 2) {
        CALL;
        FILE *file = fopen(argv[1], "wb");
        if (!file) return 3;
        size_t written = fwrite(A, sizeof(double), TENSOR_SIZE, file);
        fclose(file);
        free(A);
        return written != TENSOR_SIZE;
    }
    size_t iterations = 1;
    for (;;) {
        memset(A, 0, TENSOR_SIZE * sizeof(double));
        double start = now();
        for (size_t i = 0; i < iterations; ++i) { CALL; }
        if (now() - start >= 0.01) break;
        iterations *= 2;
    }
    printf("{\"iterations\":%zu,\"ns\":[", iterations);
    for (int sample = 0; sample < 7; ++sample) {
        memset(A, 0, TENSOR_SIZE * sizeof(double));
        double start = now();
        for (size_t i = 0; i < iterations; ++i) { CALL; }
        double elapsed = now() - start;
        printf("%s%.9g", sample ? "," : "", elapsed * 1e9 / iterations);
    }
    printf("],\"checksum\":%.17g}\n", A[0]);
    free(A);
    return 0;
}
