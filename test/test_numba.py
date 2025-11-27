# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import ctypes
import importlib
import subprocess
from pathlib import Path

import numba
import numpy as np

from ffcx.codegeneration.utils import dtype_to_scalar_dtype, numba_ufcx_kernel_signature


def test_poisson():
    subprocess.run(
        ["ffcx", Path(__file__).parent / "laplace.py", "--language", "numba"], check=True
    )

    laplace = importlib.import_module("laplace_numba")

    print(dir(laplace))

    def wrapper(scalar_type, real_type):
        c_signature = numba_ufcx_kernel_signature(scalar_type, real_type)
        return numba.cfunc(c_signature, nopython=True)

    dtype = np.float64
    realtype = np.dtype(dtype(0).real)
    kernel_a = wrapper(dtype, realtype)(laplace.form_laplace_a.form_integrals[0].tabulate_tensor)

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    empty = np.empty((0,), dtype=realtype)

    # c_type = "double"
    # c_xtype = "double"
    kernel_a(
        A.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        empty.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        empty.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        0,
    )

    A_expected = np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)

    assert np.allclose(A, np.trace(kappa_value) * A_expected)

    kernel_L = wrapper(dtype, realtype)(laplace.form_laplace_L.form_integrals[0].tabulate_tensor)

    b = np.zeros((3,), dtype=dtype)
    w = np.full((3,), 0.5, dtype=dtype)
    c = np.empty((0,), dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    empty = np.empty((0,), dtype=realtype)

    # c_type = "double"
    # c_xtype = "double"
    kernel_L(
        b.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        w.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        coords.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        empty.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        empty.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        0,
    )

    b_expected = np.full((3,), 1 / 6, dtype=np.float64)

    assert np.allclose(b, 0.5 * b_expected)
