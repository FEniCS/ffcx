# Copyright (C) 2018 Chris Richardson
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest
import cffi

import ffc.codegeneration.jit
import ufl


def float_to_type(name):
    """Map a string name to C and NumPy types"""
    if name == "double":
        return "double", np.float64
    elif name == "double complex":
        return "double _Complex", np.complex128
    else:
        raise RuntimeError("Unknown C type for: {}".format(name))


@pytest.fixture(scope="module")
def lagrange_element():
    """Compile list of Lagrange elements"""
    cell = ufl.tetrahedron
    elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
    compiled_elements, module = ffc.codegeneration.jit.compile_elements(elements)
    return elements, compiled_elements, module


@pytest.fixture(scope="module")
def hexahedral_element():
    """Compile list of Lagrange elements"""
    cell = ufl.hexahedron
    elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
    compiled_elements, module = ffc.codegeneration.jit.compile_elements(elements)
    return elements, compiled_elements, module


def test_dim_degree(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        assert compiled_e[0].geometric_dimension == 3
        assert compiled_e[0].topological_dimension == 3
        assert e.degree() == compiled_e[0].degree


def test_tabulate_reference_dof_coordinates(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        tdim = compiled_e[0].topological_dimension
        space_dim = compiled_e[0].space_dimension
        X = np.zeros([space_dim, tdim])
        X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
        compiled_e[0].tabulate_reference_dof_coordinates(X_ptr)


def test_evaluate_reference_basis_hex(hexahedral_element):
    ufl_elements, compiled_elements, module = hexahedral_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        space_dim = compiled_e[0].space_dimension
        X = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5], [1.0, 0.0, 0.0],
                      [0.0, 1.0, 0.0], [1.0, 1.0, 1.0]])
        npoint = X.shape[0]
        X_ptr = module.ffi.cast("const double *", module.ffi.from_buffer(X))
        vals = np.zeros([npoint, space_dim])
        vals_ptr = module.ffi.cast("double *", module.ffi.from_buffer(vals))
        compiled_e[0].evaluate_reference_basis(vals_ptr, npoint, X_ptr)
        assert np.isclose(np.sum(vals), npoint)


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[0.5, -1 / 6, -1 / 6, -1 / 6],
                         [-1 / 6, 1 / 6, 0.0, 0.0],
                         [-1 / 6, 0.0, 1 / 6, 0.0],
                         [-1 / 6, 0.0, 0.0, 1 / 6]], dtype=np.float64)),
    ("double complex",
     np.array(
         [[0.5 + 0j, -1 / 6 + 0j, -1 / 6 + 0j, -1 / 6 + 0j],
          [-1 / 6 + 0j, 1 / 6 + 0j, 0.0 + 0j, 0.0 + 0j],
          [-1 / 6 + 0j, 0.0 + 0j, 1 / 6 + 0j, 0.0 + 0j],
          [-1 / 6 + 0j, 0.0 + 0j, 0.0 + 0j, 1 / 6 + 0j]],
         dtype=np.complex128)),
])
def test_laplace_bilinear_form_3d(mode, expected_result):
    cell = ufl.tetrahedron
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(forms, parameters={'scalar_type': mode})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)

    c_type, np_type = float_to_type(mode)
    A = np.zeros((4, 4), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    assert np.allclose(A, expected_result)
