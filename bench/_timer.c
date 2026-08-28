/* Copyright (C) 2026 Garth N. Wells
 *
 * This file is part of FFCx. (https://www.fenicsproject.org)
 *
 * SPDX-License-Identifier:    LGPL-3.0-or-later
 *
 * Tiny C timing loop for bench/harness.py. Calls a tabulate_tensor_float64
 * kernel (passed by address, since it lives in a separately-compiled shared
 * library built by the FFCx JIT) `nrepeat` times back to back, entirely in
 * C, so the reported per-call cost isn't swamped by Python/FFI call
 * overhead -- which would otherwise dominate for the cheapest kernels (a
 * plain P1 mass matrix is tens of nanoseconds).
 */
#include <stdint.h>
#include <string.h>
#include <time.h>

typedef void (*tabulate_tensor_float64)(double *restrict A, const double *restrict w,
                                         const double *restrict c,
                                         const double *restrict coordinate_dofs,
                                         const int *restrict entity_local_index,
                                         const unsigned char *restrict quadrature_permutation,
                                         void *custom_data);

/* Returns nanoseconds per call, averaged over nrepeat calls. A is zeroed
 * before each call since every kernel accumulates into it (+=). */
double bench_run(uintptr_t fn_addr, double *A, long A_size, const double *w, const double *c,
                  const double *coordinate_dofs, long nrepeat) {
  tabulate_tensor_float64 fn = (tabulate_tensor_float64)fn_addr;

  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  for (long i = 0; i < nrepeat; ++i) {
    memset(A, 0, sizeof(double) * (size_t)A_size);
    fn(A, w, c, coordinate_dofs, 0, 0, 0);
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);

  double ns = (double)(t1.tv_sec - t0.tv_sec) * 1e9 + (double)(t1.tv_nsec - t0.tv_nsec);
  return ns / (double)nrepeat;
}
