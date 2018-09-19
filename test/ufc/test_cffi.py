# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFC (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest
import cffi

import ffc.backends.ufc.jit
import ufl


@pytest.fixture(scope="module")
def lagrange_element():
    """Compile list of Lagrange elements"""
    cell = ufl.triangle
    elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
    compiled_elements, module = ffc.backends.ufc.jit.compile_elements(elements)
    return elements, compiled_elements, module


def test_dim_degree(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        assert compiled_e.geometric_dimension == 2
        assert compiled_e.topological_dimension == 2
        assert e.degree() == compiled_e.degree


def test_tabulate_reference_dof_coordinates(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        # test = module.ffi.string(compiled_e.family)
        tdim = compiled_e.topological_dimension
        space_dim = compiled_e.space_dimension
        X = np.zeros([space_dim, tdim])
        X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
        compiled_e.tabulate_reference_dof_coordinates(X_ptr)
        # print(X)


def test_form():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TestFunction(element), ufl.TrialFunction(element)
    a = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.backends.ufc.jit.compile_forms(forms)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_default_cell_integral()
    A = np.zeros((3, 3), dtype=np.float64)
    w0 = np.array([], dtype=np.float64)
    w1 = np.array([w0.ctypes.data], dtype=np.uint64)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0,
                       1.0, 0.0,
                       0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(ffi.cast('double  *', A.ctypes.data),
                          ffi.cast('double  * *', w1.ctypes.data),
                          ffi.cast('double  *', coords.ctypes.data), 0)
    A_analytic = np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)
    A_diff = (A - A_analytic)
    assert A_diff.max() == 0
    assert A_diff.min() == 0

def test_form_coefficient():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TestFunction(element), ufl.TrialFunction(element)
    g = ufl.Coefficient(element)
    a = g * ufl.inner(u, v) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.backends.ufc.jit.compile_forms(forms)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_default_cell_integral()
    A = np.zeros((3, 3), dtype=np.float64)
    w0 = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    w1 = np.array([w0.ctypes.data], dtype=np.uint64)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0,
                       1.0, 0.0,
                       0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(ffi.cast('double  *', A.ctypes.data),
                          ffi.cast('double  * *', w1.ctypes.data),
                          ffi.cast('double  *', coords.ctypes.data), 0)

    A_analytic = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64)/24.0
    A_diff = (A - A_analytic)
    assert np.isclose(A_diff.max(), 0.0)
    assert np.isclose(A_diff.min(), 0.0)

def test_complex():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.backends.ufc.jit.compile_forms(forms, parameters={'scalar_type': 'double complex'})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_default_cell_integral()
    A = np.zeros((3, 3), dtype=np.complex128)
    w0 = np.array([], dtype=np.complex128)
    w1 = np.array([w0.ctypes.data], dtype=np.uint64)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0,
                       1.0, 0.0,
                       0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(ffi.cast('double _Complex *', A.ctypes.data),
                          ffi.cast('double _Complex * *', w1.ctypes.data),
                          ffi.cast('double *', coords.ctypes.data), 0)
    print(A)
    A_analytic = np.array([[1.0+0j, -0.5+0j, -0.5+0j],
                           [-0.5+0j, 0.5+0j, 0.0+0j],
                           [-0.5+0j, 0.0+0j, 0.5+0j]], dtype=np.complex128)
    A_diff = (A - A_analytic)
    assert A_diff.max() == 0
    assert A_diff.min() == 0



# cell = ufl.triangle
# elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
# compiled_elements, module = ffc.backends.ufc.jit.compile_elements(elements)

# for e, compiled_e in zip(elements, compiled_elements):
#     assert compiled_e.geometric_dimension == 2
#     assert compiled_e.topological_dimension == 2
#     assert e.degree() == compiled_e.degree
