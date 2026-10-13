# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import math

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
import ffcx.compiler
import ffcx.options
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


def generate_kernel(forms, dtype, options, compile_args=[], cache_dir=None):
    """Generate kernel for given forms."""
    if cache_dir is None:
        # use a different cache directory for each option
        sf = options.get("sum_factorization", False)
        cache_dir = f"./ffcx-cache-{sf}"

    compiled_forms, module, code = ffcx.codegeneration.jit.compile_forms(
        forms, cache_dir=cache_dir, options=options, cffi_extra_compile_args=compile_args
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
    kernel, _code, module = generate_kernel([a], dtype, options={"scalar_type": dtype})
    ffi = module.ffi
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )

    # Use sum factorization
    A1 = np.zeros((ndofs, ndofs), dtype=dtype)
    kernel, _code, module = generate_kernel(
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
        ffi.NULL,
    )

    np.testing.assert_allclose(A, A1, rtol=1e-6, atol=1e-6)


def distorted_cell(cell_type, dtype):
    """Vertices of a non-affine cell."""
    if cell_type == basix.CellType.quadrilateral:
        coords = [[0.0, 0.0, 0.0], [1.2, 0.1, 0.0], [0.1, 0.9, 0.0], [1.3, 1.4, 0.0]]
    else:
        coords = [
            [x + 0.1 * y * z, y + 0.05 * x, z + 0.1 * x * y]
            for z in (0, 1)
            for y in (0, 1)
            for x in (0, 1)
        ]
    return np.array(coords, dtype=dtype)


def sum_factorisation_forms(cell_type, P):
    """Forms covering the cases of a sum-factorised contraction.

    Returns a map from name to form and the number of coefficient values.
    """
    gdim = cell_to_gdim(cell_type)
    variant = basix.LagrangeVariant.gll_warped
    mesh = ufl.Mesh(create_tensor_product_element(cell_type, 1, variant, shape=(gdim,)))
    V = ufl.FunctionSpace(mesh, create_tensor_product_element(cell_type, P, variant))
    W = ufl.FunctionSpace(mesh, create_tensor_product_element(cell_type, P, variant, (gdim,)))
    Q = ufl.FunctionSpace(mesh, create_tensor_product_element(cell_type, P + 1, variant))
    DG0 = ufl.FunctionSpace(mesh, basix.ufl.element("DG", cell_type.name, 0))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    uw, vw = ufl.TrialFunction(W), ufl.TestFunction(W)
    f = ufl.Coefficient(V)
    x = ufl.SpatialCoordinate(mesh)
    nf = V.ufl_element().dim
    return {
        "mass": (u * v * ufl.dx, 0),
        "laplace": (ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx, 0),
        "mixed_derivative": (u.dx(0) * v.dx(gdim - 1) * ufl.dx, 0),
        "coefficient_mass": (ufl.inner(f * u, v) * ufl.dx, nf),
        "sin_coefficient": (ufl.sin(f) * u.dx(0) * v * ufl.dx, nf),
        "vector_laplace": (ufl.inner(ufl.grad(uw), ufl.grad(vw)) * ufl.dx, 0),
        "elasticity": (
            ufl.inner(ufl.sym(ufl.grad(uw)), ufl.sym(ufl.grad(vw))) * ufl.dx
            + ufl.div(uw) * ufl.div(vw) * ufl.dx,
            0,
        ),
        "mixed_spaces": (ufl.TrialFunction(Q) * v * ufl.dx, 0),
        "linear": (ufl.sin(x[0]) * f * v * ufl.dx, nf),
        "linear_gradient": (ufl.inner(ufl.grad(f), ufl.grad(v)) * ufl.dx, nf),
        "non_tensor_product_argument": (ufl.TrialFunction(DG0) * v * ufl.dx, 0),
    }


def tabulate(form, num_coefficients, cell_type, dtype, options, compile_args, cache_dir):
    """Tabulate the cell integral of a form on a distorted cell."""
    kernel, _code, module = generate_kernel(
        [form], dtype, {"scalar_type": dtype, **options}, compile_args, cache_dir
    )
    xdtype = dtype_to_scalar_dtype(dtype)
    coords = distorted_cell(cell_type, xdtype)
    A = np.zeros(
        math.prod(a.ufl_function_space().ufl_element().dim for a in form.arguments()), dtype
    )
    w = np.linspace(0.3, 1.7, max(num_coefficients, 1)).astype(dtype)
    c = np.zeros(1, dtype=dtype)
    ffi = module.ffi
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", w.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coords.ctypes.data),
        ffi.NULL,
        ffi.NULL,
        ffi.NULL,
    )
    return A


def uses_sum_factorised_contraction(form, dtype):
    """Check if the generated code contracts one direction at a time."""
    options = ffcx.options.get_options(
        {"scalar_type": dtype, "sum_factorization": True, "sum_factorization_contraction": True}
    )
    code, _ = ffcx.compiler.compile_ufl_objects([form], options)
    return "_g1_" in "".join(code)


@pytest.mark.parametrize(
    "name",
    [
        "mass",
        "laplace",
        "mixed_derivative",
        "coefficient_mass",
        "sin_coefficient",
        "vector_laplace",
        "elasticity",
        "mixed_spaces",
        "linear",
        "linear_gradient",
        "non_tensor_product_argument",
    ],
)
@pytest.mark.parametrize("P", [1, 3])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorised_contraction(name, P, cell_type, compile_args, tmp_path):
    dtype = "float64"
    form, num_coefficients = sum_factorisation_forms(cell_type, P)[name]
    args = (form, num_coefficients, cell_type, dtype)
    A_ref = tabulate(*args, {}, compile_args, tmp_path / "ref")
    A = tabulate(
        *args,
        {"sum_factorization": True, "sum_factorization_contraction": True},
        compile_args,
        tmp_path / "sf",
    )
    np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12 * np.max(np.abs(A_ref)))

    # A table without a tensor factorisation falls back to the fused loop
    assert uses_sum_factorised_contraction(form, dtype) == (name != "non_tensor_product_argument")


@pytest.mark.parametrize("name", ["coefficient_mass", "vector_laplace"])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorised_contraction_complex(name, cell_type, compile_args, tmp_path):
    dtype = "complex128"
    form, num_coefficients = sum_factorisation_forms(cell_type, 2)[name]
    args = (form, num_coefficients, cell_type, dtype)
    A_ref = tabulate(*args, {}, compile_args, tmp_path / "ref")
    A = tabulate(
        *args,
        {"sum_factorization": True, "sum_factorization_contraction": True},
        compile_args,
        tmp_path / "sf",
    )
    np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12 * np.max(np.abs(A_ref)))
    assert uses_sum_factorised_contraction(form, dtype)


@pytest.mark.parametrize("quadrature_degree,expected", [(8, True), (30, False)])
def test_sum_factorised_contraction_stack_budget(quadrature_degree, expected):
    """Rules whose cached integrand doesn't fit the stack keep the fused loop."""
    form, _ = sum_factorisation_forms(basix.CellType.hexahedron, 1)["laplace"]
    integrand = form.integrals()[0].integrand()
    a = integrand * ufl.dx(metadata={"quadrature_degree": quadrature_degree})
    assert uses_sum_factorised_contraction(a, "float64") == expected


@pytest.mark.parametrize("name", ["two_quadrature_degrees", "piecewise_coefficient", "linear"])
@pytest.mark.parametrize(
    "options",
    [
        {"sum_factorization": True},
        {"sum_factorization": True, "sum_factorization_contraction": True},
    ],
)
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_quadrature_rules(name, options, cell_type, compile_args, tmp_path):
    """Integrals with several or one point rules.

    Each rule has its own 1D tables, and a one point rule makes the
    tables of coefficients piecewise.
    """
    gdim = cell_to_gdim(cell_type)
    variant = basix.LagrangeVariant.gll_warped
    mesh = ufl.Mesh(create_tensor_product_element(cell_type, 1, variant, shape=(gdim,)))
    V = ufl.FunctionSpace(mesh, create_tensor_product_element(cell_type, 2, variant))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)
    f = ufl.Coefficient(V)
    one_point = ufl.dx(metadata={"quadrature_degree": 1})
    form, num_coefficients = {
        "two_quadrature_degrees": (
            u * v * ufl.dx(degree=2) + ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx(degree=5),
            0,
        ),
        "piecewise_coefficient": (f * u * v * one_point, V.ufl_element().dim),
        "linear": (f * v * one_point, V.ufl_element().dim),
    }[name]
    args = (form, num_coefficients, cell_type, "float64")
    A_ref = tabulate(*args, {}, compile_args, tmp_path / "ref")
    A = tabulate(*args, options, compile_args, tmp_path / "sf")
    np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12 * np.max(np.abs(A_ref)))
