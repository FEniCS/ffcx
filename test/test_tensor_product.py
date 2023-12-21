# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest

import basix.ufl
import ffcx.codegeneration.jit
import ufl
from ffcx.codegeneration.utils import cdtype_to_numpy, scalar_to_value_type


def create_tensor_product_element(cell_type, degree, variant, shape=None):
    """Create tensor product element."""
    if cell_type == basix.CellType.quadrilateral:
        gdim = 2
    elif cell_type == basix.CellType.hexahedron:
        gdim = 3
    else:
        raise NotImplementedError("Only quadrilateral and hexahedron supported")

    family = basix.ElementFamily.P
    ref = basix.create_element(family, cell_type, degree, variant)
    factors = ref.get_tensor_product_representation()[0]
    perm = factors[1]
    dof_ordering = np.argsort(perm)
    element = basix.create_element(family, cell_type, degree, variant,
                                   dof_ordering=dof_ordering)
    uflelement = basix.ufl._BasixElement(element, gdim=gdim)
    if shape is None:
        return uflelement
    else:
        return basix.ufl.blocked_element(uflelement, shape=shape, gdim=gdim)


def generate_kernel(forms, mode, options):

    # string to a different cache folder for sum factorization
    sf = options.get("sum_factorization", False)
    cache_dir = f"./ffcx-cache-{sf}"

    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(forms, cache_dir=cache_dir, options=options)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0]

    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1

    default_integral = form0.form_integrals[offsets[cell]]

    np_type = cdtype_to_numpy(mode)
    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")

    return kernel, code, module


@pytest.mark.parametrize("mode, P", [("double", 1), ("double", 2), ("double", 3), ("double", 4),
                                     ("float", 1), ("float", 2), ("float", 3), ("float", 4)])
def test_bilinear_form(mode, P):
    cell_type = basix.CellType.quadrilateral
    element = create_tensor_product_element(cell_type, P, basix.LagrangeVariant.gll_warped)
    coords = create_tensor_product_element(cell_type, 1, basix.LagrangeVariant.gll_warped, shape=(2, ))
    mesh = ufl.Mesh(coords)
    V = ufl.FunctionSpace(mesh, element)

    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    np_type = cdtype_to_numpy(mode)
    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    ndofs = element.dim

    A = np.zeros((ndofs, ndofs), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0],
                       [1.0, 1.0, 0.0]], dtype=np_gtype)

    kernel, code, module = generate_kernel([a], mode, options={"scalar_type": mode})
    ffi = module.ffi
    kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    # Use sum factorization
    A1 = np.zeros((ndofs, ndofs), dtype=np_type)
    kernel, code, module = generate_kernel([a], mode, options={"scalar_type": mode, "sum_factorization": True})
    ffi = module.ffi
    kernel(ffi.cast('{type} *'.format(type=mode), A1.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    np.testing.assert_allclose(A, A1, rtol=1e-6, atol=1e-6)
