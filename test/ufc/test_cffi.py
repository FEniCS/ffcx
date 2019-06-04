# -*- coding: utf-8 -*-
# Copyright (C) 2018 Garth N. Wells
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
    elif name == "float":
        return "float", np.float32
    elif name == "float complex":
        return "float _Complex", np.complex64
    elif name == "long double":
        return "long double", np.longdouble
    else:
        raise RuntimeError("Unknown C type for: {}".format(name))


@pytest.fixture(scope="module")
def lagrange_element():
    """Compile list of Lagrange elements"""
    cell = ufl.triangle
    elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
    compiled_elements, module = ffc.codegeneration.jit.compile_elements(elements)
    return elements, compiled_elements, module


@pytest.fixture(scope="module")
def quadrilateral_element():
    """Compile list of Quadrilateral elements"""
    cell = ufl.quadrilateral
    elements = [ufl.FiniteElement("Lagrange", cell, p) for p in range(1, 5)]
    compiled_elements, module = ffc.codegeneration.jit.compile_elements(elements)
    return elements, compiled_elements, module


def test_dim_degree(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        assert compiled_e[0].geometric_dimension == 2
        assert compiled_e[0].topological_dimension == 2
        assert e.degree() == compiled_e[0].degree


def test_tabulate_reference_dof_coordinates(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        # test = module.ffi.string(compiled_e.family)
        tdim = compiled_e[0].topological_dimension
        space_dim = compiled_e[0].space_dimension
        X = np.zeros([space_dim, tdim])
        X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
        compiled_e[0].tabulate_reference_dof_coordinates(X_ptr)
        # print(X)


def test_evaluate_reference_basis(lagrange_element):
    ufl_elements, compiled_elements, module = lagrange_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        space_dim = compiled_e[0].space_dimension
        X = np.array([[0.0, 0.0], [0.5, 0.5], [0.25, 0.25],
                      [1 / 3, 1 / 3], [1.0, 0.0], [0.0, 1.0]])
        npoint = X.shape[0]
        X_ptr = module.ffi.cast("const double *", module.ffi.from_buffer(X))
        vals = np.zeros([npoint, space_dim])
        vals_ptr = module.ffi.cast("double *", module.ffi.from_buffer(vals))
        compiled_e[0].evaluate_reference_basis(vals_ptr, npoint, X_ptr)
        assert np.isclose(np.sum(vals), npoint)
        np.set_printoptions(suppress=True)
        print('X=', X, 'vals = ', vals, np.sum(vals))


def test_evaluate_reference_basis_quad(quadrilateral_element):
    ufl_elements, compiled_elements, module = quadrilateral_element
    for e, compiled_e in zip(ufl_elements, compiled_elements):
        space_dim = compiled_e[0].space_dimension
        X = np.array([[0.0, 0.0], [0.5, 0.5], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
        npoint = X.shape[0]
        X_ptr = module.ffi.cast("const double *", module.ffi.from_buffer(X))
        vals = np.zeros([npoint, space_dim])
        vals_ptr = module.ffi.cast("double *", module.ffi.from_buffer(vals))
        compiled_e[0].evaluate_reference_basis(vals_ptr, npoint, X_ptr)
        assert np.isclose(np.sum(vals), npoint)
        np.set_printoptions(suppress=True)
        print('X=', X, 'vals = ', vals, np.sum(vals))


def test_cmap():
    cell = ufl.triangle
    element = ufl.VectorElement("Lagrange", cell, 1)
    mesh = ufl.Mesh(element)
    compiled_cmap, module = ffc.codegeneration.jit.compile_coordinate_maps([mesh])
    x = np.array([[0.5, 0.5]], dtype=np.float64)
    x_ptr = module.ffi.cast("double *", module.ffi.from_buffer(x))
    X = np.zeros_like(x)
    X_ptr = module.ffi.cast("double *", module.ffi.from_buffer(X))
    coords = np.array([0.0, 0.0, 2.0, 0.0, 0.0, 4.0], dtype=np.float64)
    coords_ptr = module.ffi.cast("double *", module.ffi.from_buffer(coords))
    compiled_cmap[0].compute_reference_coordinates(X_ptr, X.shape[0], x_ptr, coords_ptr, 0)
    assert(np.isclose(X[0, 0], 0.25))
    assert(np.isclose(X[0, 1], 0.125))


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)),
    ("double complex",
     np.array(
         [[1.0 + 0j, -0.5 + 0j, -0.5 + 0j], [-0.5 + 0j, 0.5 + 0j, 0.0 + 0j],
          [-0.5 + 0j, 0.0 + 0j, 0.5 + 0j]],
         dtype=np.complex128)),
])
def test_laplace_bilinear_form_2d(mode, expected_result):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode, 'cache_dir': './compile_cache'})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = cffi.FFI()
    form0 = compiled_forms[0][0]

    assert form0.num_cell_integrals == 1
    ids = np.zeros(form0.num_cell_integrals, dtype=np.int32)
    form0.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    default_integral = form0.create_cell_integral(ids[0])

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)

    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    default_integral.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), 0)

    assert np.allclose(A, expected_result)


@pytest.mark.parametrize("mode,expected_result", [
    ("float",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.float32)),
    ("long double",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.longdouble)),
    ("double",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.float64)),
    ("double complex",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.complex128)),
    ("float complex",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.complex64)),
])
def test_mass_bilinear_form_2d(mode, expected_result):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(u, v) * ufl.dx
    L = ufl.conj(v) * ufl.dx
    forms = [a, L]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)
    form1 = compiled_forms[1][0].create_cell_integral(-1)

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), 0)

    b = np.zeros(3, dtype=np_type)
    form1.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), b.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), 0)

    assert np.allclose(A, expected_result)
    assert np.allclose(b, 1.0 / 6.0)


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)
     - (1.0 / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64)),
    ("double complex",
     np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.complex128)
     - (1.0j / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.complex128)),
])
def test_helmholtz_form_2d(mode, expected_result):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    if mode == "double":
        k = 1.0
    elif mode == "double complex":
        k = ufl.constantvalue.ComplexValue(1j)
    else:
        raise RuntimeError("Unknown mode type")

    a = (ufl.inner(ufl.grad(u), ufl.grad(v)) - ufl.inner(k * u, v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), 0)

    assert np.allclose(A, expected_result)


def test_form_coefficient():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TestFunction(element), ufl.TrialFunction(element)
    g = ufl.Coefficient(element)
    a = g * ufl.inner(u, v) * ufl.dx
    forms = [a]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(forms)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)
    A = np.zeros((3, 3), dtype=np.float64)
    w = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('double  *', A.ctypes.data), ffi.cast('double  *', w.ctypes.data),
        ffi.cast('double  *', coords.ctypes.data), 0)

    A_analytic = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64) / 24.0
    A_diff = (A - A_analytic)
    assert np.isclose(A_diff.max(), 0.0)
    assert np.isclose(A_diff.min(), 0.0)


def test_subdomains():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a0 = ufl.inner(u, v) * ufl.dx + ufl.inner(u, v) * ufl.dx(2)
    a1 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx
    a2 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx(1)
    a3 = ufl.inner(u, v) * ufl.ds(210) + ufl.inner(u, v) * ufl.ds(0)
    forms = [a0, a1, a2, a3]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': 'double'})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = cffi.FFI()

    form0 = compiled_forms[0][0]
    ids = np.zeros(form0.num_cell_integrals, dtype=np.int32)
    form0.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1 and ids[1] == 2

    form1 = compiled_forms[1][0]
    ids = np.zeros(form1.num_cell_integrals, dtype=np.int32)
    form1.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1 and ids[1] == 2

    form2 = compiled_forms[2][0]
    ids = np.zeros(form2.num_cell_integrals, dtype=np.int32)
    form2.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == 1 and ids[1] == 2

    form3 = compiled_forms[3][0]
    ids = np.zeros(form3.num_cell_integrals, dtype=np.int32)
    form3.get_cell_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert len(ids) == 0
    ids = np.zeros(form3.num_exterior_facet_integrals, dtype=np.int32)
    form3.get_exterior_facet_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == 0 and ids[1] == 210


def test_interior_facet_integral():
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a0 = ufl.inner(ufl.jump(u), ufl.jump(v)) * ufl.dS
    forms = [a0]
    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': 'double'})

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = cffi.FFI()

    form0 = compiled_forms[0][0]
    ids = np.zeros(form0.num_interior_facet_integrals, dtype=np.int32)
    form0.get_interior_facet_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    integral0 = form0.create_interior_facet_integral(-1)
    A = np.zeros((6, 6), dtype=np.float64)
    w = np.array([], dtype=np.float64)
    facets = np.array([0, 2], dtype=np.int32)
    orients = np.array([1, 1], dtype=np.int32)
    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
                       1.0, 0.0, 0.0, 1.0, 1.0, 1.0], dtype=np.float64)

    integral0.tabulate_tensor(
        ffi.cast('double  *', A.ctypes.data), ffi.cast('double  *', w.ctypes.data),
        ffi.cast('double  *', coords.ctypes.data), ffi.cast('int *', facets.ctypes.data),
        ffi.cast('int *', orients.ctypes.data))

    print(A)


@pytest.mark.parametrize("mode", ["double", "double complex"])
def test_conditional(mode):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    x = ufl.SpatialCoordinate(cell)
    condition = ufl.Or(ufl.ge(ufl.real(x[0] + x[1]), 0.1),
                       ufl.ge(ufl.real(x[1] + x[1]**2), 0.1))
    c1 = ufl.conditional(condition, 2.0, 1.0)
    a = c1 * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    x1x2 = ufl.real(x[0] + ufl.as_ufl(2) * x[1])
    c2 = ufl.conditional(ufl.ge(x1x2, 0), 6.0, 0.0)
    b = c2 * ufl.conj(v) * ufl.dx

    forms = [a, b]

    compiled_forms, module = ffc.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode})

    form0 = compiled_forms[0][0].create_cell_integral(-1)
    form1 = compiled_forms[1][0].create_cell_integral(-1)

    ffi = cffi.FFI()
    c_type, np_type = float_to_type(mode)

    A1 = np.zeros((3, 3), dtype=np_type)
    w1 = np.array([1.0, 1.0, 1.0], dtype=np_type)
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)

    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A1.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w1.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), 0)

    expected_result = np.array([[2, -1, -1], [-1, 1, 0], [-1, 0, 1]], dtype=np_type)
    assert np.allclose(A1, expected_result)

    A2 = np.zeros(3, dtype=np_type)
    w2 = np.array([1.0, 1.0, 1.0], dtype=np_type)
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)

    form1.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A2.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w2.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), 0)

    expected_result = np.ones(3, dtype=np_type)
    assert np.allclose(A2, expected_result)
