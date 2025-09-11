# Copyright (C) 2018-2025 Garth N. Wells, Matthew Scroggs, Paul T. Kühner and Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import os
import sys

import basix.ufl
import numpy as np
import pytest
import sympy
import ufl
from sympy.abc import x, y, z

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


@pytest.mark.parametrize(
    "dtype,expected_result",
    [
        (
            "float64",
            np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64),
        ),
        pytest.param(
            "complex128",
            np.array(
                [
                    [1.0 + 0j, -0.5 + 0j, -0.5 + 0j],
                    [-0.5 + 0j, 0.5 + 0j, 0.0 + 0j],
                    [-0.5 + 0j, 0.0 + 0j, 0.5 + 0j],
                ],
                dtype=np.complex128,
            ),
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_laplace_bilinear_form_2d(dtype, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    kappa = ufl.Constant(domain, shape=(2, 2))
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)

    a = ufl.tr(kappa) * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

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

    assert domain.ufl_coordinate_element().basix_hash() == default_integral.coordinate_element_hash

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)

    kappa_value = np.array([[1.0, 2.0], [3.0, 4.0]])
    c = np.array(kappa_value.flatten(), dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    assert np.allclose(A, np.trace(kappa_value) * expected_result)


@pytest.mark.parametrize(
    "dtype,expected_result",
    [
        (
            np.float32,
            np.array(
                [
                    [1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0],
                ],
                dtype=np.float32,
            ),
        ),
        # ("longdouble",
        #  np.array(
        #      [[1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0], [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
        #       [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0]],
        #      dtype=np.longdouble)),
        (
            np.float64,
            np.array(
                [
                    [1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0],
                ],
                dtype=np.float64,
            ),
        ),
        pytest.param(
            np.complex128,
            np.array(
                [
                    [1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0],
                ],
                dtype=np.complex128,
            ),
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
        pytest.param(
            np.complex64,
            np.array(
                [
                    [1.0 / 12.0, 1.0 / 24.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 12.0, 1.0 / 24.0],
                    [1.0 / 24.0, 1.0 / 24.0, 1.0 / 12.0],
                ],
                dtype=np.complex64,
            ),
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_mass_bilinear_form_2d(dtype, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(u, v) * ufl.dx
    L = ufl.conj(v) * ufl.dx
    forms = [a, L]
    _compiled_forms, _module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )


@pytest.mark.parametrize(
    "dtype,expected_result",
    [
        (
            "float64",
            np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.float64)
            - (1.0 / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64),
        ),
        pytest.param(
            "complex128",
            np.array([[1.0, -0.5, -0.5], [-0.5, 0.5, 0.0], [-0.5, 0.0, 0.5]], dtype=np.complex128)
            - (1.0j / 24.0) * np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.complex128),
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_helmholtz_form_2d(dtype, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    if np.issubdtype(dtype, np.complexfloating):
        k = ufl.constantvalue.ComplexValue(1j)
    elif np.issubdtype(dtype, np.floating):
        k = 1.0
    else:
        raise RuntimeError(
            "Unknown mode type",
        )

    a = (ufl.inner(ufl.grad(u), ufl.grad(v)) - ufl.inner(k * u, v)) * ufl.dx
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    ffi = module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(form0, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    np.testing.assert_allclose(A, expected_result)


@pytest.mark.parametrize(
    "dtype,expected_result",
    [
        (
            "float64",
            np.array(
                [
                    [0.5, -1 / 6, -1 / 6, -1 / 6],
                    [-1 / 6, 1 / 6, 0.0, 0.0],
                    [-1 / 6, 0.0, 1 / 6, 0.0],
                    [-1 / 6, 0.0, 0.0, 1 / 6],
                ],
                dtype=np.float64,
            ),
        ),
        pytest.param(
            "complex128",
            np.array(
                [
                    [0.5 + 0j, -1 / 6 + 0j, -1 / 6 + 0j, -1 / 6 + 0j],
                    [-1 / 6 + 0j, 1 / 6 + 0j, 0.0 + 0j, 0.0 + 0j],
                    [-1 / 6 + 0j, 0.0 + 0j, 1 / 6 + 0j, 0.0 + 0j],
                    [-1 / 6 + 0j, 0.0 + 0j, 0.0 + 0j, 1 / 6 + 0j],
                ],
                dtype=np.complex128,
            ),
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_laplace_bilinear_form_3d(dtype, expected_result, compile_args):
    element = basix.ufl.element("Lagrange", "tetrahedron", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]

    A = np.zeros((4, 4), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    ffi = module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(form0, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    assert np.allclose(A, expected_result)


def test_form_coefficient(compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TestFunction(space), ufl.TrialFunction(space)
    g = ufl.Coefficient(space)
    a = g * ufl.inner(u, v) * ufl.dx
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    form0 = compiled_forms[0].form_integrals[0]
    A = np.zeros((3, 3), dtype=np.float64)
    w = np.array([1.0, 1.0, 1.0], dtype=np.float64)
    c = np.array([], dtype=np.float64)
    perm = np.array([0], dtype=np.uint8)

    ffi = module.ffi
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64)

    kernel = getattr(form0, "tabulate_tensor_float64")
    kernel(
        ffi.cast("double  *", A.ctypes.data),
        ffi.cast("double  *", w.ctypes.data),
        ffi.cast("double  *", c.ctypes.data),
        ffi.cast("double  *", coords.ctypes.data),
        ffi.NULL,
        ffi.cast("uint8_t *", perm.ctypes.data),
        ffi.NULL,
    )

    A_analytic = np.array([[2, 1, 1], [1, 2, 1], [1, 1, 2]], dtype=np.float64) / 24.0
    A_diff = A - A_analytic
    assert np.isclose(A_diff.max(), 0.0)
    assert np.isclose(A_diff.min(), 0.0)


def test_subdomains(compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a0 = ufl.inner(u, v) * ufl.dx + ufl.inner(u, v) * ufl.dx(2)
    a1 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx
    a2 = ufl.inner(u, v) * ufl.dx(2) + ufl.inner(u, v) * ufl.dx(1)
    a3 = ufl.inner(u, v) * ufl.ds(210) + ufl.inner(u, v) * ufl.ds(0)
    forms = [a0, a1, a2, a3]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
    )

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


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_interior_facet_integral(dtype, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a0 = ufl.inner(ufl.jump(ufl.grad(u)), ufl.jump(ufl.grad(v))) * ufl.dS
    forms = [a0]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi

    form0 = compiled_forms[0]

    ffi = module.ffi

    integral0 = form0.form_integrals[0]
    A = np.zeros((6, 6), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    facets = np.array([0, 2], dtype=np.intc)
    perms = np.array([0, 1], dtype=np.uint8)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array(
        [
            [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
            [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0],
        ],
        dtype=xdtype,
    )

    c_type = dtype_to_c_type(dtype)
    c_xtype = dtype_to_c_type(xdtype)
    kernel = getattr(integral0, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type}  *", A.ctypes.data),
        ffi.cast(f"{c_type}  *", w.ctypes.data),
        ffi.cast(f"{c_type}  *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", facets.ctypes.data),
        ffi.cast("uint8_t *", perms.ctypes.data),
        ffi.NULL,
    )


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_conditional(dtype, compile_args):
    element = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    x = ufl.SpatialCoordinate(domain)
    condition = ufl.Or(ufl.ge(ufl.real(x[0] + x[1]), 0.1), ufl.ge(ufl.real(x[1] + x[1] ** 2), 0.1))
    c1 = ufl.conditional(condition, 2.0, 1.0)
    a = c1 * ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx

    x1x2 = ufl.real(x[0] + ufl.as_ufl(2) * x[1])
    c2 = ufl.conditional(ufl.ge(x1x2, 0), 6.0, 0.0)
    b = c2 * ufl.conj(v) * ufl.dx

    forms = [a, b]

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    form0 = compiled_forms[0].form_integrals[0]
    form1 = compiled_forms[1].form_integrals[0]

    ffi = module.ffi

    A1 = np.zeros((3, 3), dtype=dtype)
    w1 = np.array([1.0, 1.0, 1.0], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel0 = ffi.cast(
        f"ufcx_tabulate_tensor_{dtype} *", getattr(form0, f"tabulate_tensor_{dtype}")
    )
    kernel0(
        ffi.cast(f"{c_type} *", A1.ctypes.data),
        ffi.cast(f"{c_type} *", w1.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    expected_result = np.array([[2, -1, -1], [-1, 1, 0], [-1, 0, 1]], dtype=dtype)
    assert np.allclose(A1, expected_result)

    A2 = np.zeros(3, dtype=dtype)
    w2 = np.array([1.0, 1.0, 1.0], dtype=dtype)

    kernel1 = ffi.cast(
        f"ufcx_tabulate_tensor_{dtype} *", getattr(form1, f"tabulate_tensor_{dtype}")
    )
    kernel1(
        ffi.cast(f"{c_type} *", A2.ctypes.data),
        ffi.cast(f"{c_type} *", w2.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    expected_result = np.ones(3, dtype=dtype)
    assert np.allclose(A2, expected_result)


def test_custom_quadrature(compile_args):
    ve = basix.ufl.element("P", "triangle", 1, shape=(2,))
    mesh = ufl.Mesh(ve)

    e = basix.ufl.element("P", mesh.ufl_cell().cellname(), 2)
    V = ufl.FunctionSpace(mesh, e)
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.5, 0.5], [0.0, 0.5], [0.5, 0.0]]
    weights = [1 / 12] * 6
    a = (
        u
        * v
        * ufl.dx(
            metadata={
                "quadrature_rule": "custom",
                "quadrature_points": points,
                "quadrature_weights": weights,
            }
        )
    )

    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form = compiled_forms[0]
    default_integral = form.form_integrals[0]

    A = np.zeros((6, 6), dtype=np.float64)
    w = np.array([], dtype=np.float64)
    c = np.array([], dtype=np.float64)

    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=np.float64)

    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(
        ffi.cast("double *", A.ctypes.data),
        ffi.cast("double *", w.ctypes.data),
        ffi.cast("double *", c.ctypes.data),
        ffi.cast("double *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    # Check that A is diagonal
    assert np.count_nonzero(A - np.diag(np.diagonal(A))) == 0


def test_curl_curl(compile_args):
    V = basix.ufl.element("N1curl", "triangle", 2)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, V)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = ufl.inner(ufl.curl(u), ufl.curl(v)) * ufl.dx

    forms = [a]
    _compiled_forms, _module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )


def lagrange_triangle_symbolic(order, corners=((1, 0), (2, 0), (0, 1)), fun=lambda i: i):
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
            eval_points += [
                tuple(S(a) + (b - a) * i for a, b in zip(p0, p1))
                for i in [(1 - 1 / sympy.sqrt(5)) / 2, (1 + 1 / sympy.sqrt(5)) / 2]
            ]
        else:
            eval_points += [
                tuple(S(a) + sympy.Rational((b - a) * i, order) for a, b in zip(p0, p1))
                for i in range(1, order)
            ]
    # face
    for f in [(0, 1, 2)]:
        p0 = corners[f[0]]
        p1 = corners[f[1]]
        p2 = corners[f[2]]
        eval_points += [
            tuple(
                S(a) + sympy.Rational((b - a) * i, order) + sympy.Rational((c - a) * j, order)
                for a, b, c in zip(p0, p1, p2)
            )
            for i in range(1, order)
            for j in range(1, order - i)
        ]

    dual_mat = [[f.subs(x, p[0]).subs(y, p[1]) for p in eval_points] for f in poly_basis]
    dual_mat = sympy.Matrix(dual_mat)
    mat = dual_mat.inv()
    functions = [sum(i * j for i, j in zip(mat.row(k), poly_basis)) for k in range(mat.rows)]
    results = []
    for f in functions:
        integrand = fun(f)
        results.append(integrand.integrate((x, 1 - y, 2 - 2 * y), (y, 0, 1)))
    return results


@pytest.mark.parametrize("dtype", ["float64"])
@pytest.mark.parametrize(
    "sym_fun,ufl_fun",
    [
        (lambda i: i, lambda i: i),
        (lambda i: i.diff(x), lambda i: ufl.grad(i)[0]),
        (lambda i: i.diff(y), lambda i: ufl.grad(i)[1]),
    ],
)
@pytest.mark.parametrize("order", [1, 2, 3])
def test_lagrange_triangle(compile_args, order, dtype, sym_fun, ufl_fun):
    sym = lagrange_triangle_symbolic(order, fun=sym_fun)
    element = basix.ufl.element("Lagrange", "triangle", order)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, element)
    v = ufl.TestFunction(space)

    a = ufl_fun(v) * ufl.dx
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]

    assert form0.form_integral_offsets[module.lib.cell + 1] == 1
    default_integral = form0.form_integrals[0]

    b = np.zeros((order + 2) * (order + 1) // 2, dtype=dtype)
    w = np.array([], dtype=dtype)
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[1.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", b.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.NULL,
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    # Check that the result is the same as for sympy
    assert np.allclose(b, [float(i) for i in sym])


def lagrange_tetrahedron_symbolic(
    order, corners=((1, 0, 0), (2, 0, 0), (0, 1, 0), (0, 0, 1)), fun=lambda i: i
):
    from sympy import S

    poly_basis = [
        x**i * y**j * z**k
        for i in range(order + 1)
        for j in range(order + 1 - i)
        for k in range(order + 1 - i - j)
    ]
    # vertices
    eval_points = [S(c) for c in corners]
    # edges
    for e in [(2, 3), (1, 3), (1, 2), (0, 3), (0, 2), (0, 1)]:
        p0 = corners[e[0]]
        p1 = corners[e[1]]
        if order > 3:
            raise NotImplementedError
        elif order == 3:
            eval_points += [
                tuple(S(a) + (b - a) * i for a, b in zip(p0, p1))
                for i in [(1 - 1 / sympy.sqrt(5)) / 2, (1 + 1 / sympy.sqrt(5)) / 2]
            ]
        else:
            eval_points += [
                tuple(S(a) + sympy.Rational((b - a) * i, order) for a, b in zip(p0, p1))
                for i in range(1, order)
            ]
    # face
    for f in [(1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)]:
        p0 = corners[f[0]]
        p1 = corners[f[1]]
        p2 = corners[f[2]]
        eval_points += [
            tuple(
                S(a) + sympy.Rational((b - a) * i, order) + sympy.Rational((c - a) * j, order)
                for a, b, c in zip(p0, p1, p2)
            )
            for i in range(1, order)
            for j in range(1, order - i)
        ]
    # interior
    for v in [(0, 1, 2, 3)]:
        p0 = corners[v[0]]
        p1 = corners[v[1]]
        p2 = corners[v[2]]
        p3 = corners[v[3]]
        eval_points += [
            tuple(
                S(a)
                + sympy.Rational((b - a) * i, order)
                + sympy.Rational((c - a) * j, order)
                + sympy.Rational((d - a) * k, order)
                for a, b, c, d in zip(p0, p1, p2, p3)
            )
            for i in range(1, order)
            for j in range(1, order - i)
            for k in range(1, order - i - j)
        ]

    dual_mat = [
        [f.subs(x, p[0]).subs(y, p[1]).subs(z, p[2]) for p in eval_points] for f in poly_basis
    ]
    dual_mat = sympy.Matrix(dual_mat)
    mat = dual_mat.inv()
    functions = [sum(i * j for i, j in zip(mat.row(k), poly_basis)) for k in range(mat.rows)]
    results = []
    for f in functions:
        integrand = fun(f)
        results.append(
            integrand.integrate((x, 1 - y - z, 2 - 2 * y - 2 * z), (y, 0, 1 - z), (z, 0, 1))
        )
    return results


@pytest.mark.parametrize("dtype", ["float64"])
@pytest.mark.parametrize(
    "sym_fun,ufl_fun",
    [
        (lambda i: i, lambda i: i),
        (lambda i: i.diff(x), lambda i: ufl.grad(i)[0]),
        (lambda i: i.diff(y), lambda i: ufl.grad(i)[1]),
    ],
)
@pytest.mark.parametrize("order", [1, 2, 3])
def test_lagrange_tetrahedron(compile_args, order, dtype, sym_fun, ufl_fun):
    sym = lagrange_tetrahedron_symbolic(order, fun=sym_fun)
    element = basix.ufl.element("Lagrange", "tetrahedron", order)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    v = ufl.TestFunction(space)

    a = ufl_fun(v) * ufl.dx
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]

    assert form0.form_integral_offsets[module.lib.cell + 1] == 1

    default_integral = form0.form_integrals[0]

    b = np.zeros((order + 3) * (order + 2) * (order + 1) // 6, dtype=dtype)
    w = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([1.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", b.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.NULL,
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    # Check that the result is the same as for sympy
    assert np.allclose(b, [float(i) for i in sym])


def test_prism(compile_args):
    element = basix.ufl.element("Lagrange", "prism", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "prism", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    v = ufl.TestFunction(space)
    L = v * ufl.dx
    forms = [L]
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": "float64"}, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]
    assert form0.form_integral_offsets[module.lib.cell + 1] == 1

    default_integral = form0.form_integrals[0]
    b = np.zeros(6, dtype=np.float64)
    coords = np.array(
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0],
        dtype=np.float64,
    )
    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(
        ffi.cast("double *", b.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.cast("double *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    assert np.isclose(sum(b), 0.5)


@pytest.mark.xfail(
    sys.platform.startswith("win32"), raises=NotImplementedError, reason="missing _Complex"
)
def test_complex_operations(compile_args):
    dtype = "complex128"
    cell = "triangle"
    c_element = basix.ufl.element("Lagrange", cell, 1, shape=(2,))
    mesh = ufl.Mesh(c_element)
    element = basix.ufl.element("DG", cell, 0, shape=(2,))
    V = ufl.FunctionSpace(mesh, element)
    u = ufl.Coefficient(V)
    J1 = ufl.real(u)[0] * ufl.imag(u)[1] * ufl.conj(u)[0] * ufl.dx
    J2 = ufl.real(u[0]) * ufl.imag(u[1]) * ufl.conj(u[0]) * ufl.dx
    forms = [J1, J2]

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    form0 = compiled_forms[0].form_integrals[0]
    form1 = compiled_forms[1].form_integrals[0]

    ffi = module.ffi
    w1 = np.array([3 + 5j, 8 - 7j], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)
    J_1 = np.zeros((1), dtype=dtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel0 = ffi.cast(
        f"ufcx_tabulate_tensor_{dtype} *", getattr(form0, f"tabulate_tensor_{dtype}")
    )
    kernel0(
        ffi.cast(f"{c_type} *", J_1.ctypes.data),
        ffi.cast(f"{c_type} *", w1.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    expected_result = np.array(
        [0.5 * np.real(w1[0]) * np.imag(w1[1]) * (np.real(w1[0]) - 1j * np.imag(w1[0]))],
        dtype=dtype,
    )
    assert np.allclose(J_1, expected_result)

    J_2 = np.zeros((1), dtype=dtype)

    kernel1 = ffi.cast(
        f"ufcx_tabulate_tensor_{dtype} *", getattr(form1, f"tabulate_tensor_{dtype}")
    )
    kernel1(
        ffi.cast(f"{c_type} *", J_2.ctypes.data),
        ffi.cast(f"{c_type} *", w1.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    assert np.allclose(J_2, expected_result)

    assert np.allclose(J_1, J_2)


def test_invalid_function_name(compile_args):
    # Monkey patch to force invalid name
    old_str = ufl.Coefficient.__str__
    ufl.Coefficient.__str__ = lambda self: "invalid function name"

    V = basix.ufl.element("Lagrange", "triangle", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    space = ufl.FunctionSpace(domain, V)
    u = ufl.Coefficient(space)
    a = ufl.inner(u, u) * ufl.dx
    forms = [a]
    try:
        _compiled_forms, _module, _code = ffcx.codegeneration.jit.compile_forms(
            forms, cffi_extra_compile_args=compile_args
        )
    except ValueError:
        pass
    except Exception:
        raise RuntimeError("Compilation should fail with ValueError.")

    # Revert monkey patch for other tests
    ufl.Coefficient.__str__ = old_str


def test_interval_vertex_quadrature(compile_args):
    c_el = basix.ufl.element("Lagrange", "interval", 1, shape=(1,))
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    dx = ufl.Measure("dx", metadata={"quadrature_rule": "vertex"})
    b = x[0] * dx

    forms = [b]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]
    assert form0.form_integral_offsets[module.lib.cell + 1] == 1

    default_integral = form0.form_integrals[0]
    J = np.zeros(1, dtype=np.float64)
    a = np.pi
    b = np.exp(1)
    coords = np.array([a, 0.0, 0.0, b, 0.0, 0.0], dtype=np.float64)

    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(
        ffi.cast("double *", J.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.cast("double *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    assert np.isclose(J[0], (0.5 * a + 0.5 * b) * np.abs(b - a))


@pytest.mark.filterwarnings("ignore:Explicitly selected vertex")
def test_facet_vertex_quadrature(compile_args):
    """Test facet vertex quadrature"""
    c_el = basix.ufl.element("Lagrange", "quadrilateral", 1, shape=(2,))
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    ds = ufl.Measure("ds", metadata={"quadrature_rule": "vertex"})
    expr = x[0] + ufl.cos(x[1])
    b1 = expr * ds
    ds_c = ufl.Measure(
        "ds",
        metadata={
            "quadrature_rule": "custom",
            "quadrature_points": np.array([[0.0], [1.0]]),
            "quadrature_weights": np.array([1.0 / 2.0, 1.0 / 2.0]),
        },
    )
    b2 = expr * ds_c
    forms = [b1, b2]
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cffi_extra_compile_args=compile_args
    )

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
        coords = np.array(
            [a, 0.1, 0.0, a + b, 0.0, 0.0, a, a, 0.0, a + 2 * b, a, 0.0], dtype=np.float64
        )
        # First facet is between vertex 0 and 1 in coords
        facets = np.array([0], dtype=np.intc)

        kernel = getattr(default_integral, "tabulate_tensor_float64")
        kernel(
            ffi.cast("double *", J.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.cast("double *", coords.ctypes.data),
            ffi.cast("int *", facets.ctypes.data),
            ffi.NULL,
            ffi.NULL,
        )
        solutions.append(J[0])
        # Test against exact result
        assert np.isclose(
            J[0], (0.5 * (a + np.cos(0.1)) + 0.5 * (a + b + np.cos(0))) * np.sqrt(b**2 + 0.1**2)
        )

    # Compare custom quadrature with vertex quadrature
    assert np.isclose(solutions[0], solutions[1])


def test_manifold_derivatives(compile_args):
    """Test higher order derivatives on manifolds"""
    c_el = basix.ufl.element("Lagrange", "interval", 1, shape=(2,))
    mesh = ufl.Mesh(c_el)

    x = ufl.SpatialCoordinate(mesh)
    dx = ufl.Measure("dx", domain=mesh)
    order = 4
    el = basix.ufl.element("Lagrange", "interval", order)
    V = ufl.FunctionSpace(mesh, el)

    u = ufl.Coefficient(V)
    d = 5.3
    f_ex = d * order * (order - 1) * x[1] ** (order - 2)
    expr = u.dx(1).dx(1) - f_ex
    J = expr * expr * dx

    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        [J], cffi_extra_compile_args=compile_args
    )

    default_integral = compiled_forms[0].form_integrals[0]
    scale = 2.5
    coords = np.array([0.0, 0.0, 0.0, 0.0, scale, 0.0], dtype=np.float64)
    dof_coords = scale * el._element.points.reshape(-1)

    w = np.array([d * d_c**order for d_c in dof_coords], dtype=np.float64)
    c = np.array([], dtype=np.float64)
    perm = np.array([0], dtype=np.uint8)

    ffi = module.ffi
    J = np.zeros(1, dtype=np.float64)
    kernel = getattr(default_integral, "tabulate_tensor_float64")
    kernel(
        ffi.cast("double *", J.ctypes.data),
        ffi.cast("double  *", w.ctypes.data),
        ffi.cast("double  *", c.ctypes.data),
        ffi.cast("double  *", coords.ctypes.data),
        ffi.NULL,
        ffi.cast("uint8_t *", perm.ctypes.data),
        ffi.NULL,
    )

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
    mesh = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    V = ufl.FunctionSpace(mesh, basix.ufl.element("Lagrange", "triangle", 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    a = (
        ufl.inner(u, v) * ufl.dx((1, 2, 3))
        + ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx(2)
        + ufl.inner(u, v) * ufl.dx
    )
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        [a], cffi_extra_compile_args=compile_args
    )
    # NOTE: This assumes that the first integral type is cell integrals, see UFCx.h
    cell = module.lib.cell
    num_integrals = (
        compiled_forms[0].form_integral_offsets[cell + 1]
        - compiled_forms[0].form_integral_offsets[cell]
    )
    assert num_integrals == 4
    unique_integrals = set(
        [
            compiled_forms[0].form_integrals[compiled_forms[0].form_integral_offsets[cell] + i]
            for i in range(num_integrals)
        ]
    )
    assert len(unique_integrals) == 2


def test_derivative_domains(compile_args):
    """Test a form with derivatives on two different domains will generate valid code."""

    V_ele = basix.ufl.element("Lagrange", "triangle", 2)
    W_ele = basix.ufl.element("Lagrange", "interval", 1)

    gdim = 2
    V_domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(gdim,)))
    W_domain = ufl.Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(gdim,)))

    V = ufl.FunctionSpace(V_domain, V_ele)
    W = ufl.FunctionSpace(W_domain, W_ele)

    u = ufl.TrialFunction(V)
    q = ufl.TestFunction(W)

    ds = ufl.Measure("ds", domain=V_domain)

    forms = [ufl.inner(u.dx(0), q.dx(0)) * ds]
    _compiled_forms, _module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": np.float64}, cffi_extra_compile_args=compile_args
    )


@pytest.mark.parametrize("dtype", ["float64"])
@pytest.mark.parametrize("permutation", [[0], [1]])
def test_mixed_dim_form(compile_args, dtype, permutation):
    """Test that the local element tensor corresponding to a mixed-dimensional form is correct.
    The form involves an integral over a facet of the cell. The trial function and a coefficient f
    are of codim 0. The test function and a coefficient g are of codim 1. We compare against another
    form where the test function and g are codim 0 but have the same trace on the facet.
    """

    def tabulate_tensor(ele_type, V_cell_type, W_cell_type, coeffs):
        "Helper function to create a form and compute the local element tensor"
        V_ele = basix.ufl.element(ele_type, V_cell_type, 2)
        W_ele = basix.ufl.element(ele_type, W_cell_type, 1)

        gdim = 2
        V_domain = ufl.Mesh(basix.ufl.element("Lagrange", V_cell_type, 1, shape=(gdim,)))
        W_domain = ufl.Mesh(basix.ufl.element("Lagrange", W_cell_type, 1, shape=(gdim,)))

        V = ufl.FunctionSpace(V_domain, V_ele)
        W = ufl.FunctionSpace(W_domain, W_ele)

        u = ufl.TrialFunction(V)
        q = ufl.TestFunction(W)

        f = ufl.Coefficient(V)
        g = ufl.Coefficient(W)

        ds = ufl.Measure("ds", domain=V_domain)

        n = ufl.FacetNormal(V_domain)
        forms = [ufl.inner(f * g * ufl.grad(u), n * q) * ds]
        compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
            forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
        )
        form0 = compiled_forms[0]
        default_integral = form0.form_integrals[0]
        kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

        A = np.zeros((W_ele.dim, V_ele.dim), dtype=dtype)
        w = np.array(coeffs, dtype=dtype)
        c = np.array([], dtype=dtype)
        facet = np.array([0], dtype=np.intc)
        perm = np.array(permutation, dtype=np.uint8)

        xdtype = dtype_to_scalar_dtype(dtype)
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]], dtype=xdtype)

        c_type = dtype_to_c_type(dtype)
        c_xtype = dtype_to_c_type(xdtype)

        ffi = module.ffi
        kernel(
            ffi.cast(f"{c_type}  *", A.ctypes.data),
            ffi.cast(f"{c_type}  *", w.ctypes.data),
            ffi.cast(f"{c_type}  *", c.ctypes.data),
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.cast("int *", facet.ctypes.data),
            ffi.cast("uint8_t *", perm.ctypes.data),
            ffi.NULL,
        )

        return A

    # Define the element type
    ele_type = "Lagrange"
    # Define the cell type for each space
    V_cell_type = "triangle"
    Vbar_cell_type = "interval"

    # Coefficient data
    # f is a quadratic on each edge that is 0 at the vertices and 1 at the midpoint
    f_data = [0, 0, 0, 1, 1, 1]
    # g is a linear function along the edge that is 0 at one vertex and 1 at the other
    g_data = [0, 1]
    # Collect coefficient data
    coeffs = f_data + g_data

    # Tabulate the tensor for the mixed-dimensional form
    A = tabulate_tensor(ele_type, V_cell_type, Vbar_cell_type, coeffs)

    # Compare to a reference result. Here, we compare against the same kernel but with
    # the interval element replaced with a triangle.
    # We create some data for g on the triangle whose trace coincides with g on the interval
    g_data = [0, 0, 1]
    coeffs_ref = f_data + g_data
    A_ref = tabulate_tensor(ele_type, V_cell_type, V_cell_type, coeffs_ref)
    # Remove the entries for the extra test DOF on the triangle element
    A_ref = A_ref[1:][:]

    # The orientation of the interval is assumed to be fixed (since it is given uniquely
    # by its global orientation in the mesh). Let us focus on local facet 0. It can have
    # two possible orientations depending on the cell topology:
    #
    # 2   1          1   1
    # | \  \         | \  \
    # |  \  \   or   |  \  \
    # 0---1  0       0---2  0
    #
    # In the second case, local facet 0 is flipped relative to the interval. If we look at
    # the degrees of freedom second-order Lagrange on the triangle and first-order Lagrange
    # on the interval, we have
    #
    # 2    1           1    1
    # | \   \          | \   \
    # 4  3   \         5  3   \
    # |    \  \    or  |    \  \
    # 0--5--1  0       0--4--2  0
    #
    # Since the trial function is defined on the triangle, the second case (where the
    # permutation is 1) is thus the same as swapping cols 1 and 2 and cols 4 and 5 of the
    # first case result.
    # NOTE: Although we choose to fix the orientation of the interval, the same result
    # can be obtained by fixing the triangle and considering flipping the interval.
    if permutation[0] == 1:
        A_ref[:, [1, 2]] = A_ref[:, [2, 1]]
        A_ref[:, [4, 5]] = A_ref[:, [5, 4]]

    assert np.allclose(A, A_ref)


@pytest.mark.parametrize("dtype", ["float64"])
def test_ds_prism(compile_args, dtype):
    element = basix.ufl.element("Lagrange", "prism", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "prism", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)

    a = ufl.inner(u, v) * ufl.ds
    forms = [a]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())

    ffi = module.ffi
    form0 = compiled_forms[0]

    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    exterior_facet = module.lib.exterior_facet
    interior_facet = module.lib.interior_facet
    assert offsets[cell + 1] - offsets[cell] == 0
    assert offsets[exterior_facet + 1] - offsets[exterior_facet] == 2
    assert offsets[interior_facet + 1] - offsets[interior_facet] == 0

    integral_id0 = form0.form_integral_ids[offsets[exterior_facet]]
    integral_id1 = form0.form_integral_ids[offsets[exterior_facet] + 1]
    assert integral_id0 == integral_id1 == -1

    integral0 = form0.form_integrals[offsets[exterior_facet]]
    integral1 = form0.form_integrals[offsets[exterior_facet] + 1]

    if basix.CellType(integral0.domain) == basix.CellType.triangle:
        assert basix.CellType(integral1.domain) == basix.CellType.quadrilateral
        integral_tri = integral0
        integral_quad = integral1
    else:
        assert basix.CellType(integral0.domain) == basix.CellType.quadrilateral
        assert basix.CellType(integral1.domain) == basix.CellType.triangle
        integral_tri = integral1
        integral_quad = integral0

    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)
    entity_perm = np.array([0], dtype=np.uint8)

    # Test integral over triangle (facet 0)
    A = np.zeros((6, 6), dtype=dtype)
    entity_index = np.array([0], dtype=int)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=xdtype,
    )

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)

    kernel = getattr(integral_tri, f"tabulate_tensor_{dtype}")

    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", entity_index.ctypes.data),
        ffi.cast("uint8_t *", entity_perm.ctypes.data),
        ffi.NULL,
    )

    assert np.allclose(
        A,
        np.array(
            [
                [1 / 12, 1 / 24, 1 / 24, 0, 0, 0],
                [1 / 24, 1 / 12, 1 / 24, 0, 0, 0],
                [1 / 24, 1 / 24, 1 / 12, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
                [0, 0, 0, 0, 0, 0],
            ]
        ),
    )

    # Test integral over quadrilateral (facet 1)
    A = np.zeros((6, 6), dtype=dtype)
    entity_index = np.array([1], dtype=np.int64)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
        ],
        dtype=xdtype,
    )

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)

    kernel = getattr(integral_quad, f"tabulate_tensor_{dtype}")

    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", entity_index.ctypes.data),
        ffi.cast("uint8_t *", entity_perm.ctypes.data),
        ffi.NULL,
    )

    assert np.allclose(
        A,
        np.array(
            [
                [1 / 9, 1 / 18, 0, 1 / 18, 1 / 36, 0],
                [1 / 18, 1 / 9, 0, 1 / 36, 1 / 18, 0],
                [0, 0, 0, 0, 0, 0],
                [1 / 18, 1 / 36, 0, 1 / 9, 1 / 18, 0],
                [1 / 36, 1 / 18, 0, 1 / 18, 1 / 9, 0],
                [0, 0, 0, 0, 0, 0],
            ]
        ),
    )


@pytest.mark.parametrize("geometry", [("interval", 1), ("triangle", 2), ("tetrahedron", 3)])
@pytest.mark.parametrize("rank", [0, 1, 2])
@pytest.mark.parametrize(
    "dtype",
    [
        np.float32,
        np.float64,
        pytest.param(
            np.complex64,
            marks=pytest.mark.skipif(
                os.name == "nt", reason="win32 platform does not support C99 _Complex numbers"
            ),
        ),
        pytest.param(
            np.complex128,
            marks=pytest.mark.skipif(
                os.name == "nt", reason="win32 platform does not support C99 _Complex numbers"
            ),
        ),
    ],
)
@pytest.mark.parametrize("element_type", ["Lagrange", "Discontinuous Lagrange"])
def test_vertex_integral(compile_args, geometry, rank, dtype, element_type):
    cell, gdim = geometry
    rdtype = np.real(dtype(0)).dtype
    domain = ufl.Mesh(basix.ufl.element("Lagrange", cell, 1, shape=(gdim,), dtype=rdtype))
    element = basix.ufl.element(element_type, cell, 1)
    dP = ufl.Measure("dP")
    x = ufl.SpatialCoordinate(domain)
    V = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    # With a Lagrange space for test and trial functions, all these forms should evaluate as the
    # x-coordinate times an indicator for alignment between the argument(s) and the vertices.
    if rank == 0:
        F = x[0] * dP
    elif rank == 1:
        F = x[0] * ufl.conj(v) * dP
    else:
        F = x[0] * u * ufl.conj(v) * dP

    if element.discontinuous:
        if rank == 0:
            # No elements involved in assembly.
            return

        with pytest.raises(TypeError):
            ffcx.codegeneration.jit.compile_forms(
                [F], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
            )
        return

    (form0,), module, _ = ffcx.codegeneration.jit.compile_forms(
        [F], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    assert form0.rank == rank

    ffi = module.ffi
    kernel = getattr(
        form0.form_integrals[0],
        f"tabulate_tensor_{dtype().dtype.name}",
    )

    vertex_count = domain.ufl_coordinate_element().basix_element.points.shape[0]

    for a, b in np.array([(0, 1), (1, 0), (2, 0), (5, -2)], dtype=rdtype):
        # General geometry for simplices of gdim 1,2 and 3.
        # gdim 1: creates the interval (a, b)
        # gdim 2: creates the triangle with vertices (a, 0), (b, 0), (0,0)
        # gdim 3: creates the tetrahedron with vertices (a, 0, 0), (b, 0, 0), (0, 0, 0), (0, 0, 1)
        coords = np.array([a, 0.0, 0.0, b, 0.0, 0.0, 0, 0, 0, 0, 0, 1], dtype=rdtype)

        for vertex in range(vertex_count):
            J = np.zeros(vertex_count**2, dtype=dtype)
            e = np.array([vertex], dtype=np.int32)
            kernel(
                ffi.cast(f"{dtype_to_c_type(dtype)} *", J.ctypes.data),
                ffi.NULL,
                ffi.NULL,
                ffi.cast(f"{dtype_to_c_type(rdtype)} *", coords.ctypes.data),
                ffi.cast("int *", e.ctypes.data),
                ffi.NULL,
                ffi.NULL,
            )
            idx = (rank > 0) * vertex + (rank > 1) * vertex * vertex_count
            assert np.isclose(coords[vertex * 3], J[idx], atol=1e2 * np.finfo(dtype).eps)
            # Check all other entries are not touched
            assert np.allclose(np.delete(J, [idx]), 0)


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_vertex_mesh(compile_args, dtype):
    """Test that a mesh with only vertices can be created and used."""

    xdtype = dtype_to_scalar_dtype(dtype)

    cell = basix.CellType.point
    c_el = basix.ufl.element("Lagrange", cell, 0, shape=(2,), discontinuous=True, dtype=xdtype)
    msh = ufl.Mesh(c_el)

    el = basix.ufl.element("Lagrange", cell, 0, discontinuous=True, dtype=xdtype)
    V = ufl.FunctionSpace(msh, el)
    u = ufl.Coefficient(V)
    dx = ufl.Measure("dx", domain=msh)
    Jh = u * dx

    forms = [Jh]
    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]
    assert form0.form_integral_offsets[module.lib.cell + 1] == 1
    default_integral = form0.form_integrals[0]
    J = np.zeros(1, dtype=dtype)
    coords = np.array([2.0, 3.0, 0.0], dtype=xdtype)
    coeffs = np.array([-5.2], dtype=dtype)
    w = np.array(coeffs, dtype=dtype)
    c = np.array([], dtype=dtype)

    c_type = dtype_to_c_type(dtype)
    c_xtype = dtype_to_c_type(xdtype)
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type}  *", J.ctypes.data),
        ffi.cast(f"{c_type}  *", w.ctypes.data),
        ffi.cast(f"{c_type}  *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    assert np.isclose(J[0], coeffs[0])


@pytest.mark.parametrize("dtype", ["float64"])
@pytest.mark.parametrize("permutation", [[0], [1]])
@pytest.mark.parametrize("local_entity_index", [0, 1, 2, 3, 4, 5])
def test_mixed_dim_form_codim2(compile_args, dtype, permutation, local_entity_index):
    """Test that the local element tensor corresponding to a mixed-dimensional
    form of codim 2 is correct. The form involves an integral over an edge of
    the cell. The trial function and a coefficient `f` are of codim 0.
    The test function is of codim 2.
    """

    def tabulate_tensor(ele_type, V_cell_type, W_cell_type, coeffs, l_i):
        V_ele = basix.ufl.element(ele_type, V_cell_type, 2)
        W_ele = basix.ufl.element(ele_type, W_cell_type, 1)

        gdim = 3
        V_domain = ufl.Mesh(basix.ufl.element("Lagrange", V_cell_type, 1, shape=(gdim,)))
        W_domain = ufl.Mesh(basix.ufl.element("Lagrange", W_cell_type, 1, shape=(gdim,)))

        V = ufl.FunctionSpace(V_domain, V_ele)
        W = ufl.FunctionSpace(W_domain, W_ele)

        u = ufl.TrialFunction(V)
        q = ufl.TestFunction(W)

        f = ufl.Coefficient(V)

        dr = ufl.Measure("dr", domain=V_domain)
        forms = [u.dx(0) * f * q * dr]

        compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
            forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
        )
        form0 = compiled_forms[0]
        default_integral = form0.form_integrals[0]
        kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

        A = np.zeros((W_ele.dim, V_ele.dim), dtype=dtype)
        w = np.array(coeffs, dtype=dtype)
        c = np.array([], dtype=dtype)
        edge = np.array([l_i], dtype=np.intc)
        perm = np.array(permutation, dtype=np.uint8)

        xdtype = dtype_to_scalar_dtype(dtype)
        coords = np.array(
            [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 1.3]], dtype=xdtype
        )

        c_type = dtype_to_c_type(dtype)
        c_xtype = dtype_to_c_type(xdtype)

        ffi = module.ffi
        kernel(
            ffi.cast(f"{c_type}  *", A.ctypes.data),
            ffi.cast(f"{c_type}  *", w.ctypes.data),
            ffi.cast(f"{c_type}  *", c.ctypes.data),
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.cast("int *", edge.ctypes.data),
            ffi.cast("uint8_t *", perm.ctypes.data),
            ffi.NULL,
        )

        return A

    # Define the element type
    ele_type = "Lagrange"
    # Define the cell type for each space
    V_cell_type = "tetrahedron"
    Vbar_cell_type = "interval"

    # Coefficient data
    f_data = [5.0, 6.0, 7.0, -1.0, 2.0, 3.0, 5.0, 5.0, 6.0, 2.0]
    coeff_data = f_data

    # Tabulate the tensor for the mixed-dimensional form
    A = tabulate_tensor(ele_type, V_cell_type, Vbar_cell_type, coeff_data, local_entity_index)

    # Compare to a reference result. Here, we compare against the same form but with
    # the interval element replaced with a triangle.
    A_ref = tabulate_tensor(ele_type, V_cell_type, V_cell_type, coeff_data, local_entity_index)

    # Map from local edge index to (local) DOF indices on that edge
    local_index_to_slice = {0: [2, 3], 1: [1, 3], 2: [1, 2], 3: [0, 3], 4: [0, 2], 5: [0, 1]}

    A_ref = A_ref[local_index_to_slice[local_entity_index]]

    # See test_mixed_dim_form for an explanation of the permutation logic. Here, for simplicity,
    # we choose to fix the orientation of the tet and flip the interval.
    if permutation[0] == 1:
        A_ref[[0, 1], :] = A_ref[[1, 0], :]

    assert np.allclose(A, A_ref)


@pytest.mark.parametrize("dtype", ["float64"])
@pytest.mark.parametrize("local_entity_index", [0, 1, 2])
def test_mixed_dim_form_codim2_2D(compile_args, dtype, local_entity_index):
    """Test that mixed assembly between 2D and 0D gives a consistent result"""

    # Define the element type
    ele_type = "Lagrange"
    # Define the cell type for each space
    V_cell_type = "triangle"
    W_cell_type = "point"

    # Coefficient data
    coeffs = [3, 4, 7, 8, 5, 2]

    # Tabulate the tensor for the mixed-dimensional form
    V_ele = basix.ufl.element(ele_type, V_cell_type, 2)
    W_ele = basix.ufl.element(ele_type, W_cell_type, 0, discontinuous=True)

    gdim = 2
    V_domain = ufl.Mesh(basix.ufl.element("Lagrange", V_cell_type, 1, shape=(gdim,)))
    W_domain = ufl.Mesh(
        basix.ufl.element("Lagrange", W_cell_type, 0, shape=(gdim,), discontinuous=True)
    )

    V = ufl.FunctionSpace(V_domain, V_ele)
    W = ufl.FunctionSpace(W_domain, W_ele)

    u = ufl.TrialFunction(V)
    q = ufl.TestFunction(W)

    f = ufl.Coefficient(V)

    dr = ufl.Measure("dr", domain=V_domain)
    x = ufl.SpatialCoordinate(V_domain)
    forms = [u.dx(0) * f * q * x[1] * dr]

    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    form0 = compiled_forms[0]
    default_integral = form0.form_integrals[0]
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

    A = np.zeros((W_ele.dim, V_ele.dim), dtype=dtype)
    w = np.array(coeffs, dtype=dtype)
    c = np.array([], dtype=dtype)
    edge = np.array([local_entity_index], dtype=np.intc)
    perm = np.zeros(1, dtype=np.uint8)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[-0.1, 0.2, 0.0], [2.5, 0.0, 0.0], [0.0, 2.0, 0.0]], dtype=xdtype)

    c_type = dtype_to_c_type(dtype)
    c_xtype = dtype_to_c_type(xdtype)

    ffi = module.ffi
    kernel(
        ffi.cast(f"{c_type}  *", A.ctypes.data),
        ffi.cast(f"{c_type}  *", w.ctypes.data),
        ffi.cast(f"{c_type}  *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", edge.ctypes.data),
        ffi.cast("uint8_t *", perm.ctypes.data),
        ffi.NULL,
    )

    # Tabulate reference solution, equivalent of assembling a vector with the same spaces from V
    # and vertex quadrature
    vertices = basix.topology(basix.cell.CellType.triangle)[0]
    vertex_coords = basix.geometry(basix.cell.CellType.triangle)
    qp = np.array(vertex_coords[vertices[local_entity_index]])
    qw = np.ones(qp.shape[0], dtype=xdtype)
    assert len(qw) == 1
    dx = ufl.Measure(
        "dx",
        domain=V_domain,
        metadata={"quadrature_rule": "custom", "quadrature_points": qp, "quadrature_weights": qw},
    )
    forms = [u.dx(0) * f * x[1] / abs(ufl.JacobianDeterminant(V_domain)) * dx]
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    form0 = compiled_forms[0]
    default_integral = form0.form_integrals[0]
    kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

    b = np.zeros(V_ele.dim, dtype=dtype)
    ffi = module.ffi
    kernel(
        ffi.cast(f"{c_type}  *", b.ctypes.data),
        ffi.cast(f"{c_type}  *", w.ctypes.data),
        ffi.cast(f"{c_type}  *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.cast("int *", edge.ctypes.data),
        ffi.cast("uint8_t *", perm.ctypes.data),
        ffi.NULL,
    )
    np.testing.assert_allclose(A[0, :], b)


@pytest.mark.parametrize("dtype", ["float64"])
@pytest.mark.parametrize("permutation", [[0], [1]])
def test_ridge_integral(compile_args, dtype, permutation):
    """Test that one can assemble a ridge integral on a single domain."""
    y_ext = 1.7
    z_ext = 1.3

    def tabulate_tensor(ele_type, V_cell_type, l_i):
        V_ele = basix.ufl.element(ele_type, V_cell_type, 1)

        gdim = 3
        domain = ufl.Mesh(basix.ufl.element("Lagrange", V_cell_type, 1, shape=(gdim,)))

        V = ufl.FunctionSpace(domain, V_ele)
        x = ufl.SpatialCoordinate(domain)
        u = ufl.TestFunction(V)

        dr = ufl.Measure("dr", domain=domain)
        forms = [u * x[2] * dr]

        compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
            forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
        )
        form0 = compiled_forms[0]
        default_integral = form0.form_integrals[0]
        kernel = getattr(default_integral, f"tabulate_tensor_{dtype}")

        A = np.zeros((V_ele.dim,), dtype=dtype)
        w = np.array([], dtype=dtype)
        c = np.array([], dtype=dtype)
        edge = np.array([l_i], dtype=np.intc)
        perm = np.array(permutation, dtype=np.uint8)

        xdtype = dtype_to_scalar_dtype(dtype)
        coords = np.array(
            [[0.0, 0.0, 0.0], [2.5, 0.0, 0.0], [0.0, y_ext, 0.0], [0.0, 0.0, z_ext]], dtype=xdtype
        )

        c_type = dtype_to_c_type(dtype)
        c_xtype = dtype_to_c_type(xdtype)

        ffi = module.ffi
        kernel(
            ffi.cast(f"{c_type}  *", A.ctypes.data),
            ffi.cast(f"{c_type}  *", w.ctypes.data),
            ffi.cast(f"{c_type}  *", c.ctypes.data),
            ffi.cast(f"{c_xtype} *", coords.ctypes.data),
            ffi.cast("int *", edge.ctypes.data),
            ffi.cast("uint8_t *", perm.ctypes.data),
            ffi.NULL,
        )

        return A

    # Define the element type
    ele_type = "Lagrange"
    # Define the cell type for each space
    V_cell_type = "tetrahedron"

    # Tabulate the tensor for the mixed-dimensional form
    b = tabulate_tensor(ele_type, V_cell_type, 0)

    # Ref first ridge integral length
    rl = np.sqrt(y_ext**2 + z_ext**2) / y_ext

    # Compute z as a function of y along the ridge
    y = sympy.Symbol("y")
    z = z_ext - z_ext / y_ext * y
    # Define the test functions along the ridge
    phi_2 = y / y_ext
    phi_3 = 1 - y / y_ext

    # Use sympy for numerical integration
    b_2 = sympy.integrate(phi_2 * z * rl, (y, 0, y_ext))
    b_3 = sympy.integrate(phi_3 * z * rl, (y, 0, y_ext))
    b_ref = np.array([0, 0, b_2, b_3], dtype=b.dtype)
    np.testing.assert_allclose(b_ref, b, atol=1e-10)


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_diagonal_form(dtype, compile_args):
    element = basix.ufl.element("Lagrange", "tetrahedron", 1)
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    space = ufl.FunctionSpace(domain, element)
    u, v = ufl.TrialFunction(space), ufl.TestFunction(space)
    a = (
        ufl.inner(u.dx(0), v.dx(2)) * ufl.dx
        - ufl.inner(u, v) * ufl.dx
        + ufl.inner(u.dx(0), v) * ufl.dx
    )
    forms = [a]
    compiled_diag_forms, diag_module, _ = ffcx.codegeneration.jit.compile_forms(
        forms,
        options={"scalar_type": dtype, "part": "diagonal"},
        cffi_extra_compile_args=compile_args,
    )

    for _, compiled_f in zip(forms, compiled_diag_forms):
        assert compiled_f.rank == 1
    diag_form0 = compiled_diag_forms[0].form_integrals[0]

    A_diag = np.zeros((4,), dtype=dtype)
    A = np.zeros((4, 4), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    ffi = diag_module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    diag_kernel = getattr(diag_form0, f"tabulate_tensor_{dtype}")
    diag_kernel(
        ffi.cast(f"{c_type} *", A_diag.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    ffi = module.ffi

    form0 = compiled_forms[0].form_integrals[0]
    kernel = getattr(form0, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    np.testing.assert_allclose(A_diag, np.diag(A))


@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_diagonal_mixed_form(dtype, compile_args):
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))

    element_u = basix.ufl.element("Lagrange", "tetrahedron", 3, shape=(3,))
    element_p = basix.ufl.element("Lagrange", "tetrahedron", 2)
    element = basix.ufl.mixed_element([element_u, element_p])
    space = ufl.FunctionSpace(domain, element)
    u, p = ufl.TrialFunctions(space)
    v, q = ufl.TestFunctions(space)
    a = (
        ufl.inner(ufl.grad(u), ufl.grad(v))
        - ufl.inner(p, ufl.div(v))
        + ufl.inner(ufl.div(u), q)
        + ufl.inner(p.dx(0), q.dx(1))
    ) * ufl.dx
    forms = [a]
    compiled_diag_forms, diag_module, _ = ffcx.codegeneration.jit.compile_forms(
        forms,
        options={"scalar_type": dtype, "part": "diagonal"},
        cffi_extra_compile_args=compile_args,
    )

    for _, compiled_f in zip(forms, compiled_diag_forms):
        assert compiled_f.rank == 1
    diag_form0 = compiled_diag_forms[0].form_integrals[0]

    A_diag = np.zeros((element.dim,), dtype=dtype)
    A = np.zeros((element.dim, element.dim), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    ffi = diag_module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    diag_kernel = getattr(diag_form0, f"tabulate_tensor_{dtype}")
    diag_kernel(
        ffi.cast(f"{c_type} *", A_diag.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    ffi = module.ffi

    form0 = compiled_forms[0].form_integrals[0]
    kernel = getattr(form0, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    np.testing.assert_allclose(A_diag, np.diag(A))


@pytest.mark.parametrize("shape", [(2,), (3,), (2, 2), (3, 3), (), (1,), (1, 1)])
@pytest.mark.parametrize(
    "dtype",
    [
        "float64",
        pytest.param(
            "complex128",
            marks=pytest.mark.xfail(
                sys.platform.startswith("win32"),
                raises=NotImplementedError,
                reason="missing _Complex",
            ),
        ),
    ],
)
def test_diagonal_mixed_block(dtype, compile_args, shape):
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))

    element_u = basix.ufl.element("Lagrange", "tetrahedron", 2, shape=shape)
    element_v = basix.ufl.element("Lagrange", "tetrahedron", 2, shape=shape)
    element = basix.ufl.mixed_element([element_u, element_v])
    space = ufl.FunctionSpace(domain, element)

    u, p = ufl.TrialFunctions(space)
    v, q = ufl.TestFunctions(space)
    a = (ufl.inner(ufl.grad(u), ufl.grad(v)) - ufl.inner(p, v) + ufl.inner(u, q)) * ufl.dx
    forms = [a]
    compiled_diag_forms, diag_module, _ = ffcx.codegeneration.jit.compile_forms(
        forms,
        options={"scalar_type": dtype, "part": "diagonal"},
        cffi_extra_compile_args=compile_args,
    )

    for _, compiled_f in zip(forms, compiled_diag_forms):
        assert compiled_f.rank == 1
    diag_form0 = compiled_diag_forms[0].form_integrals[0]

    A_diag = np.zeros((element.dim,), dtype=dtype)
    A = np.zeros((element.dim, element.dim), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    ffi = diag_module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0], dtype=xdtype)

    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    diag_kernel = getattr(diag_form0, f"tabulate_tensor_{dtype}")
    diag_kernel(
        ffi.cast(f"{c_type} *", A_diag.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    compiled_forms, module, _code = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    ffi = module.ffi

    form0 = compiled_forms[0].form_integrals[0]
    kernel = getattr(form0, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    np.testing.assert_allclose(A_diag, np.diag(A))
