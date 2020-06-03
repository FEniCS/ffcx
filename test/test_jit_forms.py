# Copyright (C) 2018 Garth N. Wells
#
# This file is part of FFCX.(https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np

import cffi
import ffcx.codegeneration.jit
import pytest
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


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)),
    ("double complex",
     np.array(
         [[1.0 + 0j, -0.5 + 0j, -0.5 + 0j], [-0.5 + 0j, 0.5 + 0j, 0.0 + 0j],
          [-0.5 + 0j, 0.0 + 0j, 0.5 + 0j]],
         dtype=np.complex128)),
])
def test_laplace_bilinear_form_2d(mode, expected_result, compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    kappa = ufl.Constant(cell, shape=(2, 2))
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)

    a = ufl.tr(kappa) * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

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

    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=np_type)

    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    default_integral.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    assert np.allclose(A, np.trace(kappa_value) * expected_result)


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
def test_mass_bilinear_form_2d(mode, expected_result, compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(u, v) * ufl.dx
    L = ufl.conj(v) * ufl.dx
    forms = [a, L]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)
    form1 = compiled_forms[1][0].create_cell_integral(-1)

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    b = np.zeros(3, dtype=np_type)
    form1.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), b.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    assert np.allclose(A, expected_result)
    assert np.allclose(b, 1.0 / 6.0)


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)
     - (1.0 / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64)),
    ("double complex",
     np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.complex128)
     - (1.0j / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.complex128)),
])
def test_helmholtz_form_2d(mode, expected_result, compile_args):
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
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)

    c_type, np_type = float_to_type(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    assert np.allclose(A, expected_result)


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
def test_laplace_bilinear_form_3d(mode, expected_result, compile_args):
    cell = ufl.tetrahedron
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

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
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    assert np.allclose(A, expected_result)


def test_form_coefficient(compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TestFunction(element), ufl.TrialFunction(element)
    g = ufl.Coefficient(element)
    a = g * ufl.inner(u, v) * ufl.dx
    forms = [a]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(forms, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0][0].create_cell_integral(-1)
    A = np.zeros((3, 3), dtype=np.float64)
    w = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    c = np.array([], dtype=np.float64)
    perm = np.array([0], dtype=np.int)

    ffi = cffi.FFI()
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)
    form0.tabulate_tensor(
        ffi.cast('double  *', A.ctypes.data),
        ffi.cast('double  *', w.ctypes.data),
        ffi.cast('double  *', c.ctypes.data),
        ffi.cast('double  *', coords.ctypes.data), ffi.NULL,
        ffi.cast('uint8_t *', perm.ctypes.data), 0)

    A_analytic = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64) / 24.0
    A_diff = (A - A_analytic)
    assert np.isclose(A_diff.max(), 0.0)
    assert np.isclose(A_diff.min(), 0.0)


def test_subdomains(compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a0 = ufl.inner(u, v) * ufl.dx + ufl.inner(u, v) * ufl.dx(2)
    a1 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx
    a2 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx(1)
    a3 = ufl.inner(u, v) * ufl.ds(210) + ufl.inner(u, v) * ufl.ds(0)
    forms = [a0, a1, a2, a3]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': 'double'}, cffi_extra_compile_args=compile_args)

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


@pytest.mark.parametrize("mode", ["double", "double complex"])
def test_interior_facet_integral(mode, compile_args):
    cell = ufl.triangle
    element = ufl.FiniteElement("Lagrange", cell, 1)
    u, v = ufl.TrialFunction(element), ufl.TestFunction(element)
    a0 = ufl.inner(ufl.jump(ufl.grad(u)), ufl.jump(ufl.grad(v))) * ufl.dS
    forms = [a0]
    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = cffi.FFI()

    form0 = compiled_forms[0][0]
    ids = np.zeros(form0.num_interior_facet_integrals, dtype=np.int32)
    form0.get_interior_facet_integral_ids(ffi.cast('int *', ids.ctypes.data))
    assert ids[0] == -1

    ffi = cffi.FFI()
    c_type, np_type = float_to_type(mode)

    integral0 = form0.create_interior_facet_integral(-1)
    A = np.zeros((6, 6), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np.float64)

    facets = np.array([0, 2], dtype=np.int32)
    perms = np.array([0, 1], dtype=np.int32)

    coords = np.array([[0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
                       [1.0, 0.0, 0.0, 1.0, 1.0, 1.0]], dtype=np.float64)

    integral0.tabulate_tensor(
        ffi.cast('{}  *'.format(c_type), A.ctypes.data),
        ffi.cast('{}  *'.format(c_type), w.ctypes.data),
        ffi.cast('{}  *'.format(c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.cast('int *', facets.ctypes.data),
        ffi.cast('uint8_t *', perms.ctypes.data), 0)


@pytest.mark.parametrize("mode", ["double", "double complex"])
def test_conditional(mode, compile_args):
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

    compiled_forms, module = ffcx.codegeneration.jit.compile_forms(
        forms, parameters={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    form0 = compiled_forms[0][0].create_cell_integral(-1)
    form1 = compiled_forms[1][0].create_cell_integral(-1)

    ffi = cffi.FFI()
    c_type, np_type = float_to_type(mode)

    A1 = np.zeros((3, 3), dtype=np_type)
    w1 = np.array([1.0, 1.0, 1.0], dtype=np_type)
    c = np.array([], dtype=np.float64)

    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)

    form0.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A1.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w1.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    expected_result = np.array([[2, -1, -1], [-1, 1, 0], [-1, 0, 1]], dtype=np_type)
    assert np.allclose(A1, expected_result)

    A2 = np.zeros(3, dtype=np_type)
    w2 = np.array([1.0, 1.0, 1.0], dtype=np_type)
    coords = np.array([0.0, 0.0, 1.0, 0.0, 0.0, 1.0], dtype=np.float64)

    form1.tabulate_tensor(
        ffi.cast('{type} *'.format(type=c_type), A2.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), w2.ctypes.data),
        ffi.cast('{type} *'.format(type=c_type), c.ctypes.data),
        ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL, 0)

    expected_result = np.ones(3, dtype=np_type)
    assert np.allclose(A2, expected_result)
