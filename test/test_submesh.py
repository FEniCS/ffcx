# Copyright (C) 2024 Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from __future__ import annotations

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


def compute_tensor(forms: list[ufl.form.Form], dtype: str, compile_args: list[str]):
    """Helper-function to compute matrix for a P1-Lagrange problem"""
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]
    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1

    default_integral = form0.form_integrals[offsets[cell]]

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[1.0, 2.0, 0.0], [1.5, 2.3, 0.0], [6.0, 1.8, 0.0]], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
    )
    return A


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        "complex128",
    ],
)
def test_multiple_mesh_codim0(dtype, compile_args):
    # Define coordinate element and element used in parent and sub-mesh
    element = basix.ufl.element("Lagrange", "triangle", 1)
    coordinate_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))

    domain = ufl.Mesh(coordinate_element)
    space = ufl.FunctionSpace(domain, element)
    u_parent = ufl.TrialFunction(space)

    # Create submesh and functionspace on submesh
    sub_domain = ufl.Mesh(coordinate_element)
    subspace = ufl.FunctionSpace(sub_domain, element)
    v_sub = ufl.TestFunction(subspace)

    #
    a = ufl.inner(u_parent.dx(0), v_sub.dx(0)) * ufl.dx(domain=domain)

    A = compute_tensor([a], dtype, compile_args)

    # Compute reference solution on with test and trial function from same mesh
    v_parent = ufl.TestFunction(space)
    a_org = ufl.inner(u_parent.dx(0), v_parent.dx(0)) * ufl.dx(domain=domain)
    A_org = compute_tensor([a_org], dtype, compile_args)

    np.testing.assert_allclose(A, A_org)
