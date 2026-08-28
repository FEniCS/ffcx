# Copyright (C) 2026 Garth N. Wells
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Compile a benchmark case and time its generated cell kernel.

Deliberately compiles at -O3 -march=native rather than the correctness
tests' -O1 (see test/conftest.py) -- the vectoriser behaviour these
benchmarks are meant to catch is invisible at -O1. Timing itself happens
in a tight C loop (bench/_timer.c), not a Python loop, so per-call FFI
overhead doesn't swamp the cheapest kernels.
"""

import ctypes
import re
import subprocess
import sys
import sysconfig
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

import ffcx.codegeneration
import ffcx.codegeneration.jit
from bench.kernels import BenchCase

_UFCX_INCLUDE_DIR = Path(ffcx.codegeneration.__file__).parent

_HERE = Path(__file__).parent
_CACHE_DIR = _HERE / ".cache"
_COMPILE_ARGS = ["-O3", "-march=native"]
_NREPEAT = 20_000


@dataclass
class BenchResult:
    """Timing and static stack-usage result for one benchmark case."""

    name: str
    ns_per_call: float
    stack_bytes: int | None
    tensor_size: int


def _build_timer_library() -> ctypes.CDLL:
    """Compile (once, cached) and load bench/_timer.c's timing loop."""
    _CACHE_DIR.mkdir(exist_ok=True)
    src = _HERE / "_timer.c"
    lib_path = _CACHE_DIR / ("_timer" + (".dylib" if sys.platform == "darwin" else ".so"))
    if not lib_path.exists() or lib_path.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(
            ["cc", *_COMPILE_ARGS, "-shared", "-fPIC", str(src), "-o", str(lib_path)],
            check=True,
        )
    lib = ctypes.CDLL(str(lib_path))
    lib.bench_run.restype = ctypes.c_double
    lib.bench_run.argtypes = [
        ctypes.c_size_t,  # fn_addr
        ctypes.POINTER(ctypes.c_double),  # A
        ctypes.c_long,  # A_size
        ctypes.POINTER(ctypes.c_double),  # w
        ctypes.POINTER(ctypes.c_double),  # c
        ctypes.POINTER(ctypes.c_double),  # coordinate_dofs
        ctypes.c_long,  # nrepeat
    ]
    return lib


def _tensor_size(case: BenchCase) -> int:
    """Local element tensor size: product of each argument's element dimension."""
    size = 1
    for arg in case.form.arguments():
        size *= arg.ufl_function_space().ufl_element().dim
    return size


def _stack_usage_bytes(source_path: Path, cellname: str) -> int | None:
    """Peak static stack usage (bytes) of the kernel's tabulate_tensor function.

    Recompiles the exact generated source with -fstack-usage into a scratch
    object file purely to read the .su sidecar; the object itself is discarded.
    Returns None if the compiler doesn't support -fstack-usage or nothing
    matches (best-effort, not required for timing to work).
    """
    with tempfile.TemporaryDirectory() as tmp:
        obj = Path(tmp) / "kernel.o"
        # The generated .c file is cffi's API-mode module source, which embeds
        # the FFCx kernel code alongside Python.h-based extension glue -- so
        # this recompile needs Python's include dir even though nothing here
        # touches the Python API.
        include = sysconfig.get_path("include")
        result = subprocess.run(
            [
                "cc",
                *_COMPILE_ARGS,
                f"-I{include}",
                f"-I{_UFCX_INCLUDE_DIR}",
                "-fstack-usage",
                "-c",
                str(source_path),
                "-o",
                str(obj),
            ],
            cwd=tmp,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0:
            return None
        su_path = obj.with_suffix(".su")
        if not su_path.exists():
            return None
        best = None
        for line in su_path.read_text().splitlines():
            m = re.match(r".*:(tabulate_tensor_\S*" + re.escape(cellname) + r")\t(\d+)\t", line)
            if m:
                size = int(m.group(2))
                best = size if best is None else max(best, size)
        return best


def run_case(case: BenchCase, timer: ctypes.CDLL, nrepeat: int = _NREPEAT) -> BenchResult:
    """Compile, JIT-build, and time one benchmark case."""
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        [case.form],
        options={"scalar_type": "float64"},
        cache_dir=_CACHE_DIR,
        cffi_extra_compile_args=_COMPILE_ARGS,
    )
    compiled_form = compiled_forms[0]
    ffi = module.ffi
    cell = module.lib.cell

    offsets = compiled_form.form_integral_offsets
    integral = compiled_form.form_integrals[offsets[cell]]
    fn_addr = int(ffi.cast("uintptr_t", integral.tabulate_tensor_float64))

    a_size = _tensor_size(case)
    A = np.zeros(a_size, dtype=np.float64)
    w = np.ascontiguousarray(case.w, dtype=np.float64)
    c = np.ascontiguousarray(case.c, dtype=np.float64)
    coords = np.ascontiguousarray(case.coords, dtype=np.float64)

    def ptr(arr: np.ndarray) -> ctypes._Pointer:
        return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    ns_per_call = timer.bench_run(fn_addr, ptr(A), a_size, ptr(w), ptr(c), ptr(coords), nrepeat)

    module_name = module.__name__.rsplit(".", 1)[-1]
    source_path = _CACHE_DIR / f"{module_name}.c"
    stack_bytes = _stack_usage_bytes(source_path, case.cellname) if source_path.exists() else None

    return BenchResult(case.name, ns_per_call, stack_bytes, tensor_size=a_size)
