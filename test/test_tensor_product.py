# Copyright (C) 2023 Igor A. Baratta
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


def cell_to_gdim(cell_type):
    """Return geometric dimension of cell."""
    if cell_type == basix.CellType.quadrilateral:
        return 2
    elif cell_type == basix.CellType.hexahedron:
        return 3
    else:
        raise NotImplementedError


def create_tensor_product_element(cell_type, degree, variant, shape=None):
    """Create tensor product element."""
    family = basix.ElementFamily.P
    element = basix.create_tp_element(family, cell_type, degree, variant)
    uflelement = basix.ufl.wrap_element(element)
    if shape is None:
        return uflelement
    else:
        return basix.ufl.blocked_element(uflelement, shape=shape)


def generate_kernel(forms, dtype, options):
    """Generate kernel for given forms."""
    # use a different cache directory for each option
    sf = options.get("sum_factorization", False)
    cache_dir = f"./ffcx-cache-{sf}"

    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, cache_dir=cache_dir, options=options
    )
    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0]

    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1
    default_integral = form0.form_integrals[offsets[cell]]
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    return kernel, code, module


@pytest.mark.parametrize("dtype", ["float32", "float64"])
@pytest.mark.parametrize("P", [1, 2, 3])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_bilinear_form(dtype, P, cell_type):
    gdim = cell_to_gdim(cell_type)
    element = create_tensor_product_element(cell_type, P, basix.LagrangeVariant.gll_warped)
    coords = create_tensor_product_element(
        cell_type, 1, basix.LagrangeVariant.gll_warped, shape=(gdim,)
    )
    mesh = ufl.Mesh(coords)
    V = ufl.FunctionSpace(mesh, element)

    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    ndofs = element.dim

    A = np.zeros((ndofs, ndofs), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    if cell_type == basix.CellType.quadrilateral:
        coords = np.array(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]], dtype=xdtype
        )
    elif cell_type == basix.CellType.hexahedron:
        coords = np.array(
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [1.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [1.0, 0.0, 1.0],
                [0.0, 1.0, 1.0],
                [1.0, 1.0, 1.0],
            ],
            dtype=xdtype,
        )

    c_type = dtype_to_c_type(dtype)
    c_xtype = dtype_to_c_type(xdtype)
    kernel, code, module = generate_kernel([a], dtype, options={"scalar_type": dtype})
    ffi = module.ffi
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
    )

    # Use sum factorization
    A1 = np.zeros((ndofs, ndofs), dtype=dtype)
    kernel, code, module = generate_kernel(
        [a], dtype, options={"scalar_type": dtype, "sum_factorization": True}
    )
    ffi = module.ffi
    kernel(
        ffi.cast(f"{c_type} *", A1.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
    )

    np.testing.assert_allclose(A, A1, rtol=1e-6, atol=1e-6)
