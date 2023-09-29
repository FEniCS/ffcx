# Copyright (C) 2018-2020 Garth N. Wells & Matthew Scroggs
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import numpy as np
import pytest
import sympy
from sympy.abc import x, y, z

import basix.ufl
import ffcx.codegeneration.jit
import ufl
from ffcx.codegeneration.utils import cdtype_to_numpy, scalar_to_value_type


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)),
    ("double _Complex",
     np.array(
         [[1.0 + 0j, -0.5 + 0j, -0.5 + 0j], [-0.5 + 0j, 0.5 + 0j, 0.0 + 0j],
          [-0.5 + 0j, 0.0 + 0j, 0.5 + 0j]],
         dtype=np.complex128)),
])
def test_laplace_bilinear_form_2d(mode, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    kappa = ufl.Constant(domain, shape=(2, 2))
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)

    a = ufl.tr(kappa) * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1

    default_integral = form0.form_integrals[offsets[cell]]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)

    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)

    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")

    kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

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
    ("double _Complex",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.complex128)),
    ("float _Complex",
     np.array(
         [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
          [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
         dtype=np.complex64)),
])
def test_mass_bilinear_form_2d(mode, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(u, v) * ufl.dx
    L = ufl.conj(v) * ufl.dx
    forms = [a, L]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]
    form1 = compiled_forms[1].form_integrals[0]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    ffi = module.ffi
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)

    kernel0 = ffi.cast(f"ufcx_tabulate_tensor_{np_type} *", getattr(form0, f"tabulate_tensor_{np_type}"))
    kernel0(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
            ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    b = np.zeros(3, dtype=np_type)
    kernel1 = ffi.cast(f"ufcx_tabulate_tensor_{np_type} *", getattr(form1, f"tabulate_tensor_{np_type}"))
    kernel1(ffi.cast('{type} *'.format(type=mode), b.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
            ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    assert np.allclose(A, expected_result)
    assert np.allclose(b, 1.0 / 6.0)


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)
     - (1.0 / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64)),
    ("double _Complex",
     np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.complex128)
     - (1.0j / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.complex128)),
])
def test_helmholtz_form_2d(mode, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    if mode == "double":
        k = 1.0
    elif mode == "double _Complex":
        k = ufl.constantvalue.ComplexValue(1j)
    else:
        raise RuntimeError("Unknown mode type")

    a = (ufl.inner(ufl.grad(u), ufl.grad(v)) - ufl.inner(k * u, v)) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((3, 3), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    ffi = module.ffi
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)
    kernel = getattr(form0, f"tabulate_tensor_{np_type}")

    kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    assert np.allclose(A, expected_result)


@pytest.mark.parametrize("mode,expected_result", [
    ("double", np.array([[0.5, -1 / 6, -1 / 6, -1 / 6],
                         [-1 / 6, 1 / 6, 0.0, 0.0],
                         [-1 / 6, 0.0, 1 / 6, 0.0],
                         [-1 / 6, 0.0, 0.0, 1 / 6]], dtype=np.float64)),
    ("double _Complex",
     np.array(
         [[0.5 + 0j, -1 / 6 + 0j, -1 / 6 + 0j, -1 / 6 + 0j],
          [-1 / 6 + 0j, 1 / 6 + 0j, 0.0 + 0j, 0.0 + 0j],
          [-1 / 6 + 0j, 0.0 + 0j, 1 / 6 + 0j, 0.0 + 0j],
          [-1 / 6 + 0j, 0.0 + 0j, 0.0 + 0j, 1 / 6 + 0j]],
         dtype=np.complex128)),
])
def test_laplace_bilinear_form_3d(mode, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "tetrahedron", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]

    np_type = cdtype_to_numpy(mode)
    A = np.zeros((4, 4), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    ffi = module.ffi
    coords = np.array([0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0], dtype=np_gtype)

    kernel = getattr(form0, f"tabulate_tensor_{np_type}")
    kernel(ffi.cast('{type} *'.format(type=mode), A.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    assert np.allclose(A, expected_result)


def test_form_coefficient(compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TestFunction(space), ufl.TrialFunction(space)
    g = ufl.Coefficient(space)
    a = g * ufl.inner(u, v) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(forms, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]
    A = np.zeros((3, 3), dtype=np.float64)
    w = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    c = np.array([], dtype=np.float64)
    perm = np.array([0], dtype=np.uint8)

    ffi = module.ffi
    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np.float64)

    kernel = getattr(form0, "tabulate_tensor_float64")
    kernel(ffi.cast('double  *', A.ctypes.data),
           ffi.cast('double  *', w.ctypes.data),
           ffi.cast('double  *', c.ctypes.data),
           ffi.cast('double  *', coords.ctypes.data), ffi.NULL,
           ffi.cast('uint8_t *', perm.ctypes.data))

    A_analytic = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64) / 24.0
    A_diff = (A - A_analytic)
    assert np.isclose(A_diff.max(), 0.0)
    assert np.isclose(A_diff.min(), 0.0)


def test_subdomains(compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a0 = ufl.inner(u, v) * ufl.dx + ufl.inner(u, v) * ufl.dx(2)
    a1 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx
    a2 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx(1)
    a3 = ufl.inner(u, v) * ufl.ds(210) + ufl.inner(u, v) * ufl.ds(0)
    forms = [a0, a1, a2, a3]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': 'double'}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0]
    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    ids = [form0.form_integral_ids[j] for j in range(offsets[cell], offsets[cell + 1])]
    assert ids[0] == -1 and ids[1] == 2

    form1 = compiled_forms[1]
    offsets = form1.form_integral_offsets
    ids = [form1.form_integral_ids[j] for j in range(offsets[cell], offsets[cell + 1])]
    assert ids[0] == -1 and ids[1] == 2

    form2 = compiled_forms[2]
    offsets = form2.form_integral_offsets
    ids = [form2.form_integral_ids[j] for j in range(offsets[cell], offsets[cell + 1])]
    assert ids[0] == 1 and ids[1] == 2

    form3 = compiled_forms[3]
    offsets = form3.form_integral_offsets
    assert offsets[cell + 1] - offsets[cell] == 0
    exf = module.lib.exterior_facet
    ids = [form3.form_integral_ids[j] for j in range(offsets[exf], offsets[exf + 1])]
    assert ids[0] == 0 and ids[1] == 210


@pytest.mark.parametrize("mode", ["double", "double _Complex"])
def test_interior_facet_integral(mode, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a0 = ufl.inner(ufl.jump(ufl.grad(u)), ufl.jump(ufl.grad(v))) * ufl.dS
    forms = [a0]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi

    form0 = compiled_forms[0]

    ffi = module.ffi
    np_type = cdtype_to_numpy(mode)

    integral0 = form0.form_integrals[0]
    A = np.zeros((6, 6), dtype=np_type)
    w = np.array([], dtype=np_type)
    c = np.array([], dtype=np.float64)

    facets = np.array([0, 2], dtype=np.intc)
    perms = np.array([0, 1], dtype=np.uint8)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    coords = np.array([[0.0, 0.0, 0.0,
                        1.0, 0.0, 0.0,
                        0.0, 1.0, 0.0],
                       [1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       1.0, 1.0, 0.0]], dtype=np_gtype)

    kernel = getattr(integral0, f"tabulate_tensor_{np_type}")
    kernel(ffi.cast(f'{mode}  *', A.ctypes.data),
           ffi.cast(f'{mode}  *', w.ctypes.data),
           ffi.cast(f'{mode}  *', c.ctypes.data),
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.cast('int *', facets.ctypes.data),
           ffi.cast('uint8_t *', perms.ctypes.data))


@pytest.mark.parametrize("mode", ["double", "double _Complex"])
def test_conditional(mode, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    x = ufl.SpatialCoordinate(domain)
    condition = ufl.Or(ufl.ge(ufl.real(x[0] + x[1]), 0.1),
                       ufl.ge(ufl.real(x[1] + x[1]**2), 0.1))
    c1 = ufl.conditional(condition, 2.0, 1.0)
    a = c1 * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    x1x2 = ufl.real(x[0] + ufl.as_ufl(2) * x[1])
    c2 = ufl.conditional(ufl.ge(x1x2, 0), 6.0, 0.0)
    b = c2 * ufl.conj(v) * ufl.dx

    forms = [a, b]

    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    form0 = compiled_forms[0].form_integrals[0]
    form1 = compiled_forms[1].form_integrals[0]

    ffi = module.ffi
    np_type = cdtype_to_numpy(mode)

    A1 = np.zeros((3, 3), dtype=np_type)
    w1 = np.array([1.0, 1.0, 1.0], dtype=np_type)
    c = np.array([], dtype=np.float64)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)

    kernel0 = ffi.cast(f"ufcx_tabulate_tensor_{np_type} *", getattr(form0, f"tabulate_tensor_{np_type}"))
    kernel0(ffi.cast('{type} *'.format(type=mode), A1.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), w1.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
            ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    expected_result = np.array([[2, -1, -1], [-1, 1, 0], [-1, 0, 1]], dtype=np_type)
    assert np.allclose(A1, expected_result)

    A2 = np.zeros(3, dtype=np_type)
    w2 = np.array([1.0, 1.0, 1.0], dtype=np_type)

    kernel1 = ffi.cast(f"ufcx_tabulate_tensor_{np_type} *", getattr(form1, f"tabulate_tensor_{np_type}"))
    kernel1(ffi.cast('{type} *'.format(type=mode), A2.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), w2.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
            ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    expected_result = np.ones(3, dtype=np_type)
    assert np.allclose(A2, expected_result)


def test_custom_quadrature(compile_args):
    ve = basix.ufl.element("P", "triangle", 1, shape=(2, ))
    mesh = ufl.Mesh(ve)

    e = basix.ufl.element("P", mesh.ufl_cell().cellname(), 2)
    V = ufl.FunctionSpace(mesh, e)
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.5], [0.0, 0.5], [0.5, 0.0]]
    weights = [1 / 12] * 6
    a = u * v * ufl.dx(metadata={"quadrature_rule": "custom",
                                 "quadrature_points": points, "quadrature_weights": weights})

    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(forms, cffi_extra_compile_args=compile_args)

    ffi = module.ffi
    form = compiled_forms[0]
    default_integral = form.form_integrals[0]

    A = np.zeros((6, 6), dtype=np.float64)
    w = np.array([], dtype=np.float64)
    c = np.array([], dtype=np.float64)

    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np.float64)

    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(ffi.cast("double *", A.ctypes.data),
           ffi.cast("double *", w.ctypes.data),
           ffi.cast("double *", c.ctypes.data),
           ffi.cast("double *", coords.ctypes.data), ffi.NULL, ffi.NULL)

    # Check that A is diagonal
    assert np.count_nonzero(A - np.diag(np.diagonal(A))) == 0


def test_curl_curl(compile_args):
    V = basix.ufl.element("N1curl", "triangle", 2)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, V)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx

    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(forms, cffi_extra_compile_args=compile_args)


def lagrange_triangle_symbolic(order, corners=[(1, 0), (2, 0), (0, 1)], fun=lambda i: i):
    from sympy import S
    poly_basis = [x**i * y**j for i in range(order + 1) for j in range(order + 1 - i)]
    # vertices
    eval_points = [S(c) for c in corners]
    # edges
    for e in [(1, 2), (0, 2), (0, 1)]:
        p0 = corners[e[0]]
        p1 = corners[e[1]]
        if order > 3:
            raise NotImplementedError
        elif order == 3:
            eval_points += [tuple(S(a) + (b - a) * i for a, b in zip(p0, p1))
                            for i in [(1 - 1 / sympy.sqrt(5)) / 2, (1 + 1 / sympy.sqrt(5)) / 2]]
        else:
            eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order) for a, b in zip(p0, p1))
                            for i in range(1, order)]
    # face
    for f in [(0, 1, 2)]:
        p0 = corners[f[0]]
        p1 = corners[f[1]]
        p2 = corners[f[2]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order)
                        + sympy.Rational((c - a) * j, order) for a, b, c in zip(p0, p1, p2))
                        for i in range(1, order) for j in range(1, order - i)]

    dual_mat = [[f.subs(x, p[0]).subs(y, p[1]) for p in eval_points] for f in poly_basis]
    dual_mat = sympy.Matrix(dual_mat)
    mat = dual_mat.inv()
    functions = [sum(i * j for i, j in zip(mat.row(k), poly_basis)) for k in range(mat.rows)]
    results = []
    for f in functions:
        integrand = fun(f)
        results.append(integrand.integrate((x, 1 - y, 2 - 2 * y), (y, 0, 1)))
    return results


@pytest.mark.parametrize("mode", ["double"])
@pytest.mark.parametrize("sym_fun,ufl_fun", [
    (lambda i: i, lambda i: i),
    (lambda i: i.diff(x), lambda i: ufl.grad(i)[0]),
    (lambda i: i.diff(y), lambda i: ufl.grad(i)[1])])
@pytest.mark.parametrize("order", [1, 2, 3])
def test_lagrange_triangle(compile_args, order, mode, sym_fun, ufl_fun):
    sym = lagrange_triangle_symbolic(order, fun=sym_fun)
    element = basix.ufl.element("Lagrange", "triangle", order)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, element)
    v = ufl.TestFunction(space)

    a = ufl_fun(v) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    ffi = module.ffi
    form0 = compiled_forms[0]

    assert form0.form_integral_offsets[module.lib.cell + 1] == 1
    default_integral = form0.form_integrals[0]

    np_type = cdtype_to_numpy(mode)
    b = np.zeros((order + 2) * (order + 1) // 2, dtype=np_type)
    w = np.array([], dtype=np_type)
    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)
    coords = np.array([[1.0, 0.0, 0.0],
                       [2.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")
    kernel(ffi.cast('{type} *'.format(type=mode), b.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.NULL,
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    # Check that the result is the same as for sympy
    assert np.allclose(b, [float(i) for i in sym])


def lagrange_tetrahedron_symbolic(order, corners=[(1, 0, 0), (2, 0, 0), (0, 1, 0), (0, 0, 1)], fun=lambda i: i):
    from sympy import S
    poly_basis = [
        x**i * y**j * z**k for i in range(order + 1) for j in range(order + 1 - i)
        for k in range(order + 1 - i - j)]
    # vertices
    eval_points = [S(c) for c in corners]
    # edges
    for e in [(2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)]:
        p0 = corners[e[0]]
        p1 = corners[e[1]]
        if order > 3:
            raise NotImplementedError
        elif order == 3:
            eval_points += [tuple(S(a) + (b - a) * i for a, b in zip(p0, p1))
                            for i in [(1 - 1 / sympy.sqrt(5)) / 2, (1 + 1 / sympy.sqrt(5)) / 2]]
        else:
            eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order) for a, b in zip(p0, p1))
                            for i in range(1, order)]
    # face
    for f in [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]:
        p0 = corners[f[0]]
        p1 = corners[f[1]]
        p2 = corners[f[2]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order)
                        + sympy.Rational((c - a) * j, order) for a, b, c in zip(p0, p1, p2))
                        for i in range(1, order) for j in range(1, order - i)]
    # interior
    for v in [(0, 1, 2, 3)]:
        p0 = corners[v[0]]
        p1 = corners[v[1]]
        p2 = corners[v[2]]
        p3 = corners[v[3]]
        eval_points += [tuple(S(a) + sympy.Rational((b - a) * i, order) + sympy.Rational((c - a) * j, order)
                        + sympy.Rational((d - a) * k, order) for a, b, c, d in zip(p0, p1, p2, p3))
                        for i in range(1, order) for j in range(1, order - i) for k in range(1, order - i - j)]

    dual_mat = [[f.subs(x, p[0]).subs(y, p[1]).subs(z, p[2]) for p in eval_points] for f in poly_basis]
    dual_mat = sympy.Matrix(dual_mat)
    mat = dual_mat.inv()
    functions = [sum(i * j for i, j in zip(mat.row(k), poly_basis)) for k in range(mat.rows)]
    results = []
    for f in functions:
        integrand = fun(f)
        results.append(integrand.integrate((x, 1 - y - z, 2 - 2 * y - 2 * z), (y, 0, 1 - z), (z, 0, 1)))
    return results


@pytest.mark.parametrize("mode", ["double"])
@pytest.mark.parametrize("sym_fun,ufl_fun", [
    (lambda i: i, lambda i: i),
    (lambda i: i.diff(x), lambda i: ufl.grad(i)[0]),
    (lambda i: i.diff(y), lambda i: ufl.grad(i)[1])])
@pytest.mark.parametrize("order", [1, 2, 3])
def test_lagrange_tetrahedron(compile_args, order, mode, sym_fun, ufl_fun):
    sym = lagrange_tetrahedron_symbolic(order, fun=sym_fun)
    element = basix.ufl.element("Lagrange", "tetrahedron", order)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3, )))
    space = ufl.FunctionSpace(domain, element)
    v = ufl.TestFunction(space)

    a = ufl_fun(v) * ufl.dx
    forms = [a]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    ffi = module.ffi
    form0 = compiled_forms[0]

    assert form0.form_integral_offsets[module.lib.cell + 1] == 1

    default_integral = form0.form_integrals[0]

    np_type = cdtype_to_numpy(mode)
    b = np.zeros((order + 3) * (order + 2) * (order + 1) // 6, dtype=np_type)
    w = np.array([], dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    coords = np.array([1.0, 0.0, 0.0,
                       2.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 1.0], dtype=np_gtype)

    kernel = getattr(default_integral, f"tabulate_tensor_{np_type}")
    kernel(ffi.cast('{type} *'.format(type=mode), b.ctypes.data),
           ffi.cast('{type} *'.format(type=mode), w.ctypes.data),
           ffi.NULL,
           ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    # Check that the result is the same as for sympy
    assert np.allclose(b, [float(i) for i in sym])


def test_prism(compile_args):
    element = basix.ufl.element("Lagrange", "prism", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "prism", 1, shape=(3, )))
    space = ufl.FunctionSpace(domain, element)
    v = ufl.TestFunction(space)

    L = v * ufl.dx
    forms = [L]
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': 'double'}, cffi_extra_compile_args=compile_args)

    ffi = module.ffi
    form0 = compiled_forms[0]
    assert form0.form_integral_offsets[module.lib.cell + 1] == 1

    default_integral = form0.form_integrals[0]
    b = np.zeros(6, dtype=np.float64)
    coords = np.array([1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       0.0, 0.0, 0.0,
                       1.0, 0.0, 1.0,
                       0.0, 1.0, 1.0,
                       0.0, 0.0, 1.0], dtype=np.float64)

    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(ffi.cast('double *', b.ctypes.data),
           ffi.NULL,
           ffi.NULL,
           ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    assert np.isclose(sum(b), 0.5)


def test_complex_operations(compile_args):
    mode = "double _Complex"
    cell = "triangle"
    c_element = basix.ufl.element("Lagrange", cell, 1, shape=(2, ))
    mesh = ufl.Mesh(c_element)
    element = basix.ufl.element("DG", cell, 0, shape=(2, ))
    V = ufl.FunctionSpace(mesh, element)
    u = ufl.Coefficient(V)
    J1 = ufl.real(u)[0] * ufl.imag(u)[1] * ufl.conj(u)[0] * ufl.dx
    J2 = ufl.real(u[0]) * ufl.imag(u[1]) * ufl.conj(u[0]) * ufl.dx
    forms = [J1, J2]

    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, options={'scalar_type': mode}, cffi_extra_compile_args=compile_args)

    form0 = compiled_forms[0].form_integrals[0]
    form1 = compiled_forms[1].form_integrals[0]

    ffi = module.ffi
    np_type = cdtype_to_numpy(mode)
    w1 = np.array([3 + 5j, 8 - 7j], dtype=np_type)
    c = np.array([], dtype=np_type)

    geom_type = scalar_to_value_type(mode)
    np_gtype = cdtype_to_numpy(geom_type)

    coords = np.array([[0.0, 0.0, 0.0],
                       [1.0, 0.0, 0.0],
                       [0.0, 1.0, 0.0]], dtype=np_gtype)
    J_1 = np.zeros((1), dtype=np_type)
    kernel0 = ffi.cast(f"ufcx_tabulate_tensor_{np_type} *", getattr(form0, f"tabulate_tensor_{np_type}"))
    kernel0(ffi.cast('{type} *'.format(type=mode), J_1.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), w1.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
            ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    expected_result = np.array([0.5 * np.real(w1[0]) * np.imag(w1[1])
                               * (np.real(w1[0]) - 1j * np.imag(w1[0]))], dtype=np_type)
    assert np.allclose(J_1, expected_result)

    J_2 = np.zeros((1), dtype=np_type)

    kernel1 = ffi.cast(f"ufcx_tabulate_tensor_{np_type} *", getattr(form1, f"tabulate_tensor_{np_type}"))
    kernel1(ffi.cast('{type} *'.format(type=mode), J_2.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), w1.ctypes.data),
            ffi.cast('{type} *'.format(type=mode), c.ctypes.data),
            ffi.cast(f'{geom_type} *', coords.ctypes.data), ffi.NULL, ffi.NULL)

    assert np.allclose(J_2, expected_result)

    assert np.allclose(J_1, J_2)


def test_invalid_function_name(compile_args):
    # Monkey patch to force invalid name
    old_str = ufl.Coefficient.__str__
    ufl.Coefficient.__str__ = lambda self: "invalid function name"

    V = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2, )))
    space = ufl.FunctionSpace(domain, V)
    u = ufl.Coefficient(space)
    a = ufl.inner(u, u) * ufl.dx

    forms = [a]

    try:
        compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
            forms, cffi_extra_compile_args=compile_args)
    except ValueError:
        pass
    except Exception:
        raise RuntimeError("Compilation should fail with ValueError.")

    # Revert monkey patch for other tests
    ufl.Coefficient.__str__ = old_str


def test_interval_vertex_quadrature(compile_args):

    c_el = basix.ufl.element("Lagrange", "interval", 1, shape=(1, ))
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    dx = ufl.Measure(
        "dx", metadata={"quadrature_rule": "vertex"})
    b = x[0] * dx

    forms = [b]
    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args)

    ffi = module.ffi
    form0 = compiled_forms[0]
    assert form0.form_integral_offsets[module.lib.cell + 1] == 1

    default_integral = form0.form_integrals[0]
    J = np.zeros(1, dtype=np.float64)
    a = np.pi
    b = np.exp(1)
    coords = np.array([a, 0.0, 0.0,

                       b, 0.0, 0.0], dtype=np.float64)

    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(ffi.cast('double *', J.ctypes.data),
           ffi.NULL,
           ffi.NULL,
           ffi.cast('double *', coords.ctypes.data), ffi.NULL, ffi.NULL)
    assert np.isclose(J[0], (0.5 * a + 0.5 * b) * np.abs(b - a))


def test_facet_vertex_quadrature(compile_args):
    """Test facet vertex quadrature"""
    c_el = basix.ufl.element("Lagrange", "quadrilateral", 1, shape=(2,))
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    ds = ufl.Measure(
        "ds", metadata={"quadrature_rule": "vertex"})
    expr = (x[0] + ufl.cos(x[1]))
    b1 = expr * ds
    ds_c = ufl.Measure(
        "ds",
        metadata={
            "quadrature_rule": "custom",
            "quadrature_points": np.array([[0.0], [1.0]]),
            "quadrature_weights": np.array([1.0 / 2.0, 1.0 / 2.0]),
        }
    )
    b2 = expr * ds_c
    forms = [b1, b2]
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args)

    ffi = module.ffi
    assert len(compiled_forms) == 2
    solutions = []
    for form in compiled_forms:
        offsets = form.form_integral_offsets
        exf = module.lib.exterior_facet
        assert offsets[exf + 1] - offsets[exf] == 1

        default_integral = form.form_integrals[offsets[exf]]
        J = np.zeros(1, dtype=np.float64)
        a = np.pi
        b = np.exp(1)
        coords = np.array([a, 0.1, 0.0,
                           a + b, 0.0, 0.0,
                           a, a, 0.,
                           a + 2 * b, a, 0.], dtype=np.float64)
        # First facet is between vertex 0 and 1 in coords
        facets = np.array([0], dtype=np.intc)

        kernel = getattr(default_integral, "tabulate_tensor_float64")
        kernel(ffi.cast('double *', J.ctypes.data),
               ffi.NULL,
               ffi.NULL,
               ffi.cast('double *', coords.ctypes.data),
               ffi.cast('int *', facets.ctypes.data),
               ffi.NULL)
        solutions.append(J[0])
        # Test against exact result
        assert np.isclose(J[0], (0.5 * (a + np.cos(0.1)) + 0.5 * (a + b + np.cos(0))) * np.sqrt(b**2 + 0.1**2))

    # Compare custom quadrature with vertex quadrature
    assert np.isclose(solutions[0], solutions[1])


def test_manifold_derivatives(compile_args):
    """Test higher order derivatives on manifolds"""
    c_el = basix.ufl.element("Lagrange", "interval", 1, shape=(2,), gdim=2)
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    dx = ufl.Measure("dx", domain=mesh)
    order = 4
    el = basix.ufl.element("Lagrange", "interval", order, gdim=2)
    V = ufl.FunctionSpace(mesh, el)

    u = ufl.Coefficient(V)
    d = 5.3
    f_ex = d * order * (order - 1) * x[1]**(order - 2)
    expr = u.dx(1).dx(1) - f_ex
    J = expr * expr * dx

    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        [J], cffi_extra_compile_args=compile_args)

    default_integral = compiled_forms[0].form_integrals[0]
    scale = 2.5
    coords = np.array([0.0, 0.0, 0.0, 0.0, scale, 0.0], dtype=np.float64)
    dof_coords = scale * el.element.points.reshape(-1)

    w = np.array([d * d_c**order for d_c in dof_coords], dtype=np.float64)
    c = np.array([], dtype=np.float64)
    perm = np.array([0], dtype=np.uint8)

    ffi = module.ffi
    J = np.zeros(1, dtype=np.float64)
    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(ffi.cast('double *', J.ctypes.data),
           ffi.cast('double  *', w.ctypes.data),
           ffi.cast('double  *', c.ctypes.data),
           ffi.cast('double  *', coords.ctypes.data), ffi.NULL,
           ffi.cast('uint8_t *', perm.ctypes.data))

    assert np.isclose(J[0], 0.0)


def test_integral_grouping(compile_args):
    """We group integrals with common integrands to avoid duplicated
    integration kernels. This means that `inner(u, v)*dx((1,2,3))  +
    inner(grad(u), grad(v))*dx(2) + inner(u,v)*dx` is grouped as
    1. `inner(u,v)*dx(("everywhere", 1, 3))`
    2. `(inner(grad(u), grad(v)) + inner(u, v))*dx(2)`
    Each of the forms has one generated `tabulate_tensor_*` function,
    which is referred to multiple times in `integrals_` and
    `integral_ids_`

    """
    mesh = ufl.Mesh(ufl.VectorElement("Lagrange", ufl.triangle, 1))
    V = ufl.FunctionSpace(mesh, ufl.FiniteElement("Lagrange", ufl.triangle, 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = ufl.inner(u, v) * ufl.dx((1, 2, 3)) + ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        [a], cffi_extra_compile_args=compile_args)
    # NOTE: This assumes that the first integral type is cell integrals, see UFCx.h
    cell = module.lib.cell
    num_integrals = compiled_forms[0].form_integral_offsets[cell + 1] - compiled_forms[0].form_integral_offsets[cell]
    assert num_integrals == 4
    unique_integrals = set([compiled_forms[0].form_integrals[compiled_forms[0].form_integral_offsets[cell] + i]
                            for i in range(num_integrals)])
    assert len(unique_integrals) == 2
