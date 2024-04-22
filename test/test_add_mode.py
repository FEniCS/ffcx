# Copyright (C) 2019 Chris Richardson
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


@pytest.mark.parametrize(
    "dtype",
    [
        "float32",
        "float64",
        "complex64",
        "complex128",
    ],
)
def test_additive_facet_integral(dtype, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(u, v) * ufl.ds
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    integral_offsets = form0.form_integral_offsets
    ex = module.lib.exterior_facet
    assert integral_offsets[ex + 1] - integral_offsets[ex] == 1
    integral_id = form0.form_integral_ids[integral_offsets[ex]]
    assert integral_id == -1

    default_integral = form0.form_integrals[integral_offsets[ex]]

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)
    facets = np.array([0], dtype=np.int32)
    perm = np.array([0], dtype=np.uint8)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array(
        [0.0, 2.0, 0.0, np.sqrt(3.0), -1.0, 0.0, -np.sqrt(3.0), -1.0, 0.0], dtype=xdtype
    )

    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    for i in range(3):
        facets[0] = i
        kernel(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.cast(f"{c_type} *", w.ctypes.data),
            ffi.cast(f"{c_type} *", c.ctypes.data),
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.cast("int *", facets.ctypes.data),
            ffi.cast("uint8_t *", perm.ctypes.data),
        )
        assert np.isclose(A.sum(), np.sqrt(12) * (i + 1))


@pytest.mark.parametrize(
    "dtype",
    [
        "float32",
        "float64",
        "complex64",
        "complex128",
    ],
)
def test_additive_cell_integral(dtype, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    cell = module.lib.cell
    offsets = form0.form_integral_offsets
    num_integrals = offsets[cell + 1] - offsets[cell]
    assert num_integrals == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1

    default_integral = form0.form_integrals[offsets[cell]]

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array(
        [0.0, 2.0, 0.0, np.sqrt(3.0), -1.0, 0.0, -np.sqrt(3.0), -1.0, 0.0], dtype=xdtype
    )

    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
    )

    A0 = np.array(A)
    for i in range(3):
        kernel(
            ffi.cast(f"{c_type} *", A.ctypes.data),
            ffi.cast(f"{c_type} *", w.ctypes.data),
            ffi.cast(f"{c_type} *", c.ctypes.data),
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
        )

        assert np.all(np.isclose(A, (i + 2) * A0))
