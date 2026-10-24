"""Time generated simplex kernels against a memory-traffic-only bound.

Each cell owns a 64-byte aligned record holding its coordinates, coefficients
and output tensor, so a sweep over many records streams through DRAM in the
same way as a batched assembly loop. The ``null`` kernel reads every input and
updates every output entry without arithmetic beyond a sum; its time is the
per-cell DRAM bound on the measuring core. ``hot`` times one record repeatedly
so that it stays in L1 and shows the arithmetic cost alone.

Generate kernels for each revision with ``generate.py`` into separate
directories, then compare them, for example::

    python benchmarks/streaming/roofline.py --variant head=$W/baseline \
        --variant new=$W/new --output roofline.json
"""

import argparse
import ctypes
import json
import os
import re
import subprocess
from pathlib import Path

import numpy as np

from benchmarks.elasticity import inputs

FLAGS = ["-O3", "-march=native", "-ffp-contract=fast", "-ffast-math", "-fPIC"]
NULL = """#include <stdint.h>
void null_kernel(double *A, const double *w, const double *c, const double *x,
                 const int *e, const uint8_t *p, void *d)
{
  double s = 0;
  for (int i = 0; i < NX; ++i) s += x[i];
  for (int i = 0; i < NW; ++i) s += w[i];
  for (int i = 0; i < SIZE; ++i) A[i] += s;
}
"""


def cases(root):
    """Yield simplex case names, kernel file names, inputs and checks."""
    for suite in ["poisson", "elasticity"]:
        manifest = json.loads((root / suite / "manifest-ffcx.json").read_text())["cases"]
        for case, m in manifest.items():
            if m["cell"] not in ("triangle", "tetrahedron"):
                continue
            data = np.load(root / suite / f"{case}.npz")
            yield (
                f"{suite}-{case}",
                Path(suite) / f"ffcx_{case}.c",
                np.ascontiguousarray(data["ffcx_coordinates"]),
                np.zeros((32, 0)),
                np.ascontiguousarray(data["expected"].reshape(32, -1)),
            )
    for degree in [1, 2, 3]:
        nodes = np.load(root / "hyper" / f"nodes{degree}.npy")
        ws, xs = inputs(nodes)
        for kind in ["residual", "tangent"]:
            n = 3 * len(nodes)
            size = n if kind == "residual" else n * n
            yield (f"hyper-{kind}{degree}", Path("hyper") / f"{kind}{degree}.c", xs, ws, size)


def driver(function, nx, nw, size):
    """Return C source for preparing, evaluating and timing one kernel."""
    wo = max(128, (nx * 8 + 63) // 64 * 64)
    ao = wo + (nw * 8 + 63) // 64 * 64
    stride = ao + (size * 8 + 63) // 64 * 64
    return f"""#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
void {function}(double*,const double*,const double*,const double*,const int*,const uint8_t*,void*);
static unsigned char *records;
static long cells;
double checksum;
void prepare(const double *xs, const double *ws, long bytes)
{{
  cells = bytes / {stride};
  if (posix_memalign((void **)&records, 64, cells * {stride})) abort();
  memset(records, 0, cells * {stride});
  for (long k = 0; k < cells; ++k) {{
    memcpy(records + k * {stride}, xs + (k % 32) * {nx}, {nx} * 8);
    memcpy(records + k * {stride} + {wo}, ws + (k % 32) * {nw}, {nw} * 8);
  }}
}}
void release(void) {{ free(records); records = 0; }}
void evaluate(double *A, const double *x, const double *w)
{{
  {function}(A, w, NULL, x, NULL, NULL, NULL);
}}
long cell_count(void) {{ return cells; }}
double record_error(const double *ref)
{{
  double error = 0;
  for (long k = 0; k < cells; ++k) {{
    const double *A = (const double *)(records + k * {stride} + {ao});
    for (long j = 0; j < {size}; ++j) {{
      double e = A[j] - ref[(k % 32) * {size} + j];
      e = e < 0 ? -e : e;
      if (!(e <= error)) error = e;
    }}
  }}
  return error;
}}
static double now(void)
{{
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec + 1e-9 * t.tv_nsec;
}}
double measure(void)
{{
  double check = 0, t0 = now();
  for (long k = 0; k < cells; ++k) {{
    unsigned char *r = records + k * {stride};
    double *A = (double *)(r + {ao});
    memset(A, 0, {size} * 8);
    {function}(A, (const double *)(r + {wo}), NULL, (const double *)r, NULL, NULL, NULL);
    check += A[k % {size}];
  }}
  double t = now() - t0;
  checksum = check;
  return t;
}}
double hot(long count)
{{
  unsigned char *r = records;
  double *A = (double *)(r + {ao}), check = 0, t0 = now();
  for (long k = 0; k < count; ++k) {{
    memset(A, 0, {size} * 8);
    {function}(A, (const double *)(r + {wo}), NULL, (const double *)r, NULL, NULL, NULL);
    check += A[k % {size}];
  }}
  double t = now() - t0;
  checksum = check;
  return t;
}}
"""


def build(work, name, source, function, nx, nw, size, cc, defines=()):
    """Compile a kernel and its driver as separate objects into a library."""
    out = work / cc
    out.mkdir(parents=True, exist_ok=True)
    driver_c = out / f"{name}-driver.c"
    driver_c.write_text(driver(function, nx, nw, size))
    objects = []
    for src, extra in [(source, [f"-I{source.parent}", *defines]), (driver_c, [])]:
        obj = out / f"{name}-{src.stem}.o"
        subprocess.run([cc, *FLAGS, *extra, "-c", str(src), "-o", str(obj)], check=True)
        objects.append(str(obj))
    lib = out / f"{name}.so"
    subprocess.run([cc, "-shared", *objects, "-lm", "-o", str(lib)], check=True)
    mod = ctypes.CDLL(str(lib))
    ptr = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")
    mod.prepare.argtypes = [ptr, ptr, ctypes.c_long]
    mod.evaluate.argtypes = [ptr, ptr, ptr]
    mod.record_error.argtypes = [ptr]
    mod.record_error.restype = ctypes.c_double
    mod.measure.restype = ctypes.c_double
    mod.hot.argtypes = [ctypes.c_long]
    mod.hot.restype = ctypes.c_double
    mod.cell_count.restype = ctypes.c_long
    return mod


def evaluate(mod, xs, ws, size, init=None):
    """Evaluate all validation cells, returning the increment to the output."""
    out = []
    for k in range(len(xs)):
        a0 = np.zeros(size) if init is None else np.sin(np.arange(size) + k) * 0.23
        a = a0.copy()
        mod.evaluate(a, np.ascontiguousarray(xs[k]), np.ascontiguousarray(ws[k]))
        out.append(a - a0)
    return np.array(out)


def main():
    """Validate all variants, then time them interleaved with the null kernel."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variant", action="append", required=True, help="name=directory")
    parser.add_argument("--work", type=Path, default=Path(os.environ.get("TMPDIR", "/tmp")))
    parser.add_argument("--cc", nargs="+", default=["gcc-16", "clang-22"])
    parser.add_argument("--filter", default="")
    parser.add_argument("--mib", type=int, default=512)
    parser.add_argument("--reps", type=int, default=5)
    parser.add_argument("--cpu", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    os.sched_setaffinity(0, {args.cpu})
    variants = dict(v.split("=", 1) for v in args.variant)
    variants = {k: Path(v).resolve() for k, v in variants.items()}
    first = next(iter(variants.values()))
    work = (args.work / "roofline-build").resolve()
    report = json.loads(args.output.read_text()) if args.output and args.output.exists() else {}
    for name, rel, xs, ws, expected in cases(first):
        if not re.search(args.filter, name):
            continue
        size = expected if isinstance(expected, int) else expected.shape[1]
        nx, nw = xs.shape[1], ws.shape[1]
        rows = {}
        for cc in args.cc:
            mods = {}
            reference = None if isinstance(expected, int) else expected
            null_c = work / f"{name}-null.c"
            null_c.parent.mkdir(parents=True, exist_ok=True)
            null_c.write_text(NULL)
            defines = [f"-DNX={nx}", f"-DNW={nw}", f"-DSIZE={size}"]
            mods["null"] = build(
                work, f"{name}-null", null_c, "null_kernel", nx, nw, size, cc, defines
            )
            for label, root in variants.items():
                source = root / rel
                function = re.search(r"void (tabulate_tensor_\w+)\(", source.read_text())[1]
                mod = build(work, f"{name}-{label}", source, function, nx, nw, size, cc)
                actual = evaluate(mod, xs, ws, size)
                if reference is None:
                    reference = actual
                np.testing.assert_allclose(
                    actual, reference, rtol=1e-10, atol=1e-11, err_msg=f"{name} {label} {cc}"
                )
                np.testing.assert_allclose(
                    evaluate(mod, xs, ws, size, True),
                    reference,
                    rtol=1e-10,
                    atol=1e-11,
                    err_msg=f"{name} {label} {cc} nonzero",
                )
                mods[label] = mod
            for mod in mods.values():
                mod.prepare(xs, np.ascontiguousarray(ws), args.mib * 1024**2)
                mod.measure()
            cells = mods["null"].cell_count()
            dram = {k: [] for k in mods}
            hot = {k: [] for k in mods}
            labels = list(mods)
            for rep in range(args.reps):
                for label in labels[rep % len(labels) :] + labels[: rep % len(labels)]:
                    dram[label].append(mods[label].measure() / cells * 1e9)
                    hot[label].append(mods[label].hot(200000) / 200000 * 1e9)
            for label, mod in mods.items():
                if label != "null":
                    err = mod.record_error(np.ascontiguousarray(reference))
                    assert err < 1e-9, (name, label, cc, err)
                mod.release()
                rows[f"{cc}-{label}"] = {
                    "dram_ns": min(dram[label]),
                    "dram_median_ns": float(np.median(dram[label])),
                    "hot_ns": min(hot[label]),
                    "dram_samples_ns": dram[label],
                    "hot_samples_ns": hot[label],
                    "cells": int(cells),
                    "max_absolute_tensor_error": err if label != "null" else None,
                }
        report[name] = rows
        line = [f"{name:24s}"]
        for cc in args.cc:
            bound = rows[f"{cc}-null"]["dram_ns"]
            line.append(f"{cc} null {bound:7.1f}")
            for label in variants:
                r = rows[f"{cc}-{label}"]
                line.append(
                    f"{label} {r['dram_ns']:8.1f} ({bound / r['dram_ns']:4.0%}, hot "
                    f"{r['hot_ns']:7.1f})"
                )
        print(" | ".join(line), flush=True)
        if args.output:
            args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
