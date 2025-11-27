# Copyright (C) 2025 Paul T. Kühner
#
# This file is part of FFCx.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
import ctypes
import importlib
import subprocess
from pathlib import Path

import cffi
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
    a_kernel = wrapper(dtype, realtype)(laplace.form_laplace_a.form_integrals[0].tabulate_tensor)

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    ffi = cffi.FFI()
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    empty = np.empty((0,), dtype=realtype)

    c_type = "double"
    # c_xtype = "double"
    print(a_kernel)
    print(type(ffi.cast(f"{c_type} *", A.ctypes.data)))
    print(dir(ffi.cast(f"{c_type} *", A.ctypes.data)))
    a_kernel(
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
