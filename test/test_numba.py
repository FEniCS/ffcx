# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import ctypes
import importlib
from pathlib import Path

import numpy as np
import numpy.typing as npt
import pytest

import ffcx.main
from ffcx.codegeneration.utils import dtype_to_scalar_dtype, numba_ufcx_kernel_signature

numba = pytest.importorskip("numba")

def wrap_kernel(scalar_type, real_type):
    c_signature = numba_ufcx_kernel_signature(scalar_type, real_type)
    return numba.cfunc(c_signature, nopython=True)


def as_C_array(np_array: npt.NDArray):
    dtype_C = np.ctypeslib.as_ctypes_type(np_array.dtype)
    return np_array.ctypes.data_as(ctypes.POINTER(dtype_C))


@pytest.mark.parametrize("scalar_type", ["float32", "float64"])  # TODO: complex limited by ctypes
def test_integral(scalar_type: str) -> None:
    opts = f"--language numba --scalar_type {scalar_type}"
    dir = Path(__file__).parent
    assert ffcx.main.main([str(dir / "poisson.py"), *opts.split(" ")]) == 0

    poisson = importlib.import_module("poisson_numba")

    dtype = np.dtype(scalar_type).type
    dtype_r = dtype_to_scalar_dtype(dtype)

    kernel_a = wrap_kernel(dtype, dtype_r)(poisson.form_poisson_a.form_integrals[0].tabulate_tensor)

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype_r)
    empty = np.empty((0,), dtype=dtype_r)

    kernel_a(
        as_C_array(A),
        as_C_array(w),
        as_C_array(c),
        as_C_array(coords),
        as_C_array(empty),
        as_C_array(empty),
        0,
    )
    A_expected = np.array(
        [[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=scalar_type
    )

    assert np.allclose(A, np.trace(kappa_value) * A_expected)

    kernel_L = wrap_kernel(dtype, dtype_r)(poisson.form_poisson_L.form_integrals[0].tabulate_tensor)

    b = np.zeros((3,), dtype=dtype)
    w = np.full((3,), 0.5, dtype=dtype)
    c = np.empty((0,), dtype=dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype_r)
    empty = np.empty((0,), dtype=dtype_r)

    kernel_L(
        as_C_array(b),
        as_C_array(w),
        as_C_array(c),
        as_C_array(coords),
        as_C_array(empty),
        as_C_array(empty),
        0,
    )

    b_expected = np.full((3,), 1 / 6, dtype=np.float64)
    assert np.allclose(b, 0.5 * b_expected)


@pytest.mark.parametrize("scalar_type", ["float32", "float64"])  # TODO: complex limited by ctypes
def test_expression(scalar_type: str) -> None:
    opts = f"--language numba --scalar_type {scalar_type}"
    dir = Path(__file__).parent
    assert ffcx.main.main([str(dir / "poisson.py"), *opts.split(" ")]) == 0

    poisson = importlib.import_module("poisson_numba")

    dtype = np.dtype(scalar_type).type
    dtype_r = dtype_to_scalar_dtype(dtype)

    kernel_expr = wrap_kernel(dtype, dtype_r)(poisson.expression_poisson_0.tabulate_tensor)

    e = np.zeros((6 * 3,), dtype=dtype)
    w = np.array(
        [1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11], dtype=dtype
    )
    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=dtype_r)
    empty = np.empty((0,), dtype=dtype_r)

    kernel_expr(
        as_C_array(e),
        as_C_array(w),
        as_C_array(c),
        as_C_array(coords),
        as_C_array(empty),
        as_C_array(empty),
        0,
    )
    e_expected = np.array(
        [5, 7, 8, 11, 15, 18, 14, 16, 17, 32, 36, 39, 23, 25, 26, 53, 57, 60], dtype=dtype
    )
    assert np.allclose(e, e_expected)
