# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import math
import sys

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.analysis
import ffcx.codegeneration.jit
import ffcx.ir.representation
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


def jit_compile(forms, options, compile_args, cache_dir):
    """Compile forms into a module."""
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, cache_dir=cache_dir, options=options, cffi_extra_compile_args=compile_args or []
    )
    for f, compiled_f in zip(forms, compiled_forms):
        assert compiled_f.rank == len(f.arguments())
    return compiled_forms, module


def generate_kernel(forms, dtype, options, compile_args=None, cache_dir=None):
    """Generate the kernel of the first of some forms."""
    if cache_dir is None:
        # use a different cache directory for each option
        sf = ffcx.options.get_options(options)["sum_factorization"]
        cache_dir = f"./ffcx-cache-{sf}"
    compiled_forms, module = jit_compile(forms, options, compile_args, cache_dir)
    return cell_kernel(compiled_forms[0], module, dtype), module


def cell_kernel(compiled_form, module, dtype):
    """Return the kernel of the default cell integral of a compiled form."""
    offsets = compiled_form.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    assert compiled_form.form_integral_ids[offsets[cell]] == -1
    return getattr(compiled_form.form_integrals[offsets[cell]], f"tabulate_tensor_{dtype}")


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
    kernel, module = generate_kernel(
        [a], dtype, options={"scalar_type": dtype, "sum_factorization": False}
    )
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
    kernel, module = generate_kernel(
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
    fw = ufl.Coefficient(W)
    g = ufl.Coefficient(Q)
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
        "vector_coefficient_gradient": (
            ufl.inner(ufl.grad(fw), ufl.grad(vw)) * ufl.dx,
            W.ufl_element().dim,
        ),
        "two_coefficients": (ufl.inner(f * g, v) * ufl.dx, nf + Q.ufl_element().dim),
        # The value and derivative of f share a stage of their evaluation
        "coefficient_and_derivative": ((f + f.dx(0)) * v * ufl.dx, nf),
        "non_tensor_product_argument": (ufl.TrialFunction(DG0) * v * ufl.dx, 0),
    }


def tabulate_kernel(kernel, module, form, num_coefficients, cell_type, dtype):
    """Tabulate the cell integral kernel of a form on a distorted cell."""
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


def assert_close(A, A_ref):
    """Check an element tensor matches a reference."""
    np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12 * np.max(np.abs(A_ref)))


def assert_kernels_match(kernels, form, num_coefficients, cell_type, dtype):
    """Check the kernels of a cell integral, by the fused loop and by sum factorisation, match."""
    A_ref, A = (
        tabulate_kernel(kernel, module, form, num_coefficients, cell_type, dtype)
        for kernel, module in kernels
    )
    assert_close(A, A_ref)


def assert_matches_fused_loop(args, compile_args, tmp_path, options=None):
    """Check a cell integral tabulated with sum factorisation matches the fused loop."""
    form, num_coefficients, cell_type, dtype = args
    kernels = [
        generate_kernel(
            [form], dtype, {"scalar_type": dtype, **options}, compile_args, tmp_path / cache_dir
        )
        for options, cache_dir in (({"sum_factorization": False}, "ref"), (options or {}, "sf"))
    ]
    assert_kernels_match(kernels, form, num_coefficients, cell_type, dtype)


def integrands(form, dtype, **options):
    """Return the quadrature rules and IR of the integrands of a form."""
    options = ffcx.options.get_options({"scalar_type": dtype, **options})
    analysis = ffcx.analysis.analyze_ufl_objects(
        [form], dtype, do_cancel_jacobian_products=options["do_cancel_jacobian_products"]
    )
    ir = ffcx.ir.representation.compute_ir(analysis, {}, "", options, None)
    return [
        (rule, integrand)
        for integral in ir.integrals
        for (_, rule), integrand in integral.expression.integrand.items()
    ]


def sum_factorisation_plans(form, dtype, **options):
    """Return the plans of the integrands of a form contracted one direction at a time."""
    return [
        integrand["sum_factorisation"]
        for _, integrand in integrands(form, dtype, **options)
        if integrand["sum_factorisation"] is not None
    ]


def sum_factorised(form, dtype, **options):
    """Check if the integrals of a form are contracted one direction at a time."""
    return len(sum_factorisation_plans(form, dtype, **options)) > 0


# Forms of `sum_factorisation_forms` compiled together into one module
CONTRACTION_FORMS = (
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
    "vector_coefficient_gradient",
    "two_coefficients",
    "coefficient_and_derivative",
    "non_tensor_product_argument",
)
COMPLEX_FORMS = (
    "coefficient_mass",
    "vector_laplace",
    "vector_coefficient_gradient",
    "two_coefficients",
)
VECTOR_FORMS = ("vector_laplace", "elasticity")


def cell_degrees(degrees):
    """Parametrise over cells and degrees, grouping tests that share compiled forms.

    With ``pytest -n <workers> --dist loadgroup``, one worker compiles each module.
    """
    return [
        pytest.param(cell_type, P, marks=pytest.mark.xdist_group(f"{cell_type.name}_P{P}"))
        for cell_type in (basix.CellType.quadrilateral, basix.CellType.hexahedron)
        for P in degrees
    ]


@pytest.fixture(scope="module")
def compiled_forms(tmp_path_factory, compile_args):
    """Compile forms of `sum_factorisation_forms` into one module, once for all tests using them."""
    modules = {}

    def compile(names, cell_type, P, dtype, sum_factorization):
        key = (names, cell_type, P, dtype, sum_factorization)
        if key not in modules:
            forms = sum_factorisation_forms(cell_type, P)
            options = {"scalar_type": dtype, "sum_factorization": sum_factorization}
            compiled, module = jit_compile(
                [forms[name][0] for name in names],
                options,
                compile_args,
                tmp_path_factory.mktemp("forms"),
            )
            modules[key] = forms, compiled, module
        return modules[key]

    return compile


def check_form(compiled_forms, names, name, P, cell_type, dtype):
    """Check a form of `sum_factorisation_forms` matches the fused loop.

    Returns whether it is contracted one direction at a time.
    """
    kernels = []
    for sum_factorization in (False, True):
        forms, compiled, module = compiled_forms(names, cell_type, P, dtype, sum_factorization)
        kernels.append((cell_kernel(compiled[names.index(name)], module, dtype), module))
    form, num_coefficients = forms[name]
    assert_kernels_match(kernels, form, num_coefficients, cell_type, dtype)
    return sum_factorised(form, dtype)


@pytest.mark.parametrize("sum_factorization,expected", [(None, True), (True, True), (False, False)])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_default_and_opt_out(sum_factorization, expected, cell_type):
    form, _ = sum_factorisation_forms(cell_type, 2)["mass"]
    options = {} if sum_factorization is None else {"sum_factorization": sum_factorization}
    assert sum_factorised(form, "float64", **options) == expected


@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_diagonal(cell_type):
    """Assembling only the diagonal uses ordinary quadrature."""
    form, _ = sum_factorisation_forms(cell_type, 2)["laplace"]
    assert not sum_factorised(form, "float64", part="diagonal")


@pytest.mark.parametrize("rule", ["custom", "vertex"])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_unfactorised_rule(rule, cell_type, compile_args, tmp_path):
    """A rule without tensor factors uses ordinary quadrature by default."""
    form, _ = sum_factorisation_forms(cell_type, 2)["mass"]
    metadata = {"quadrature_rule": rule, "quadrature_degree": 1}
    if rule == "custom":
        points, weights = basix.make_quadrature(cell_type, 4)
        metadata.update(quadrature_points=points, quadrature_weights=weights)
    form = form.integrals()[0].integrand() * ufl.dx(metadata=metadata)
    assert_matches_fused_loop((form, 0, cell_type, "float64"), compile_args, tmp_path)
    assert not sum_factorised(form, "float64")


@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_unfactorised_coefficient(cell_type, compile_args, tmp_path):
    """Ordinary coefficient and geometry tables work with tensor-product argument tables."""
    gdim = cell_to_gdim(cell_type)
    mesh = ufl.Mesh(basix.ufl.element("P", cell_type, 1, shape=(gdim,)))
    V = ufl.FunctionSpace(
        mesh, create_tensor_product_element(cell_type, 2, basix.LagrangeVariant.gll_warped)
    )
    Q = ufl.FunctionSpace(mesh, basix.ufl.element("P", cell_type, 2))
    w, v = ufl.Coefficient(Q), ufl.TestFunction(V)
    form = (w * v + ufl.inner(ufl.grad(w), ufl.grad(v))) * ufl.dx
    args = (form, Q.ufl_element().dim, cell_type, "float64")
    assert_matches_fused_loop(args, compile_args, tmp_path)
    assert sum_factorised(form, "float64")


@pytest.mark.parametrize("name", CONTRACTION_FORMS)
@pytest.mark.parametrize("cell_type,P", cell_degrees([1, 2, 3]))
def test_sum_factorised_contraction(name, P, cell_type, compiled_forms):
    contracted = check_form(compiled_forms, CONTRACTION_FORMS, name, P, cell_type, "float64")
    # A table without a tensor factorisation falls back to the fused loop
    assert contracted == (name != "non_tensor_product_argument")


@pytest.mark.xfail(
    sys.platform.startswith("win32"), raises=NotImplementedError, reason="missing _Complex"
)
@pytest.mark.parametrize("name", COMPLEX_FORMS)
@pytest.mark.parametrize("cell_type,P", cell_degrees([1, 2]))
def test_sum_factorised_contraction_complex(name, cell_type, P, compiled_forms):
    assert check_form(compiled_forms, COMPLEX_FORMS, name, P, cell_type, "complex128")


@pytest.mark.parametrize("name", VECTOR_FORMS)
@pytest.mark.parametrize("cell_type,P", cell_degrees([2, 4]))
def test_sum_factorised_contraction_vector(name, P, cell_type, compiled_forms):
    """Vector-valued forms, whose blocks share contributions, at even degree.

    An even degree gives an odd number of dofs and points per direction.
    """
    assert check_form(compiled_forms, VECTOR_FORMS, name, P, cell_type, "float64")


@pytest.mark.parametrize(
    "quadrature_degree,dtype",
    [
        (8, "float64"),
        (22, "float64"),
        (22, "complex128"),
        (34, "float64"),
    ],
)
def test_sum_factorised_contraction_large_rule(quadrature_degree, dtype):
    """Without a stack limit, large rules keep directional contractions."""
    form, _ = sum_factorisation_forms(basix.CellType.hexahedron, 1)["laplace"]
    integrand = form.integrals()[0].integrand()
    a = integrand * ufl.dx(metadata={"quadrature_degree": quadrature_degree})
    assert sum_factorised(a, dtype, sum_factorization_stack_limit=0)


@pytest.mark.parametrize(
    "cell_type,name,limit",
    [
        (basix.CellType.quadrilateral, "elasticity", 2048),
        (basix.CellType.quadrilateral, "elasticity", 1024),
        (basix.CellType.quadrilateral, "vector_laplace", 1024),
        (basix.CellType.hexahedron, "elasticity", 16384),
    ],
)
def test_sum_factorisation_batches(cell_type, name, limit, compile_args, tmp_path):
    """A stack limit splits the integrand factors into batches, contracted in turn."""
    form, num_coefficients = sum_factorisation_forms(cell_type, 2)[name]
    options = {"sum_factorization_stack_limit": limit}
    (plan,) = sum_factorisation_plans(form, "float64", **options)
    assert len(plan.batches) > 1
    args = (form, num_coefficients, cell_type, "float64")
    assert_matches_fused_loop(args, compile_args, tmp_path, options)


@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_stack_limit_fallback(cell_type, compile_args, tmp_path):
    """Ordinary quadrature is used if one integrand factor doesn't fit the stack limit."""
    form, num_coefficients = sum_factorisation_forms(cell_type, 2)["elasticity"]
    options = {"sum_factorization_stack_limit": 256}
    assert not sum_factorised(form, "float64", **options)
    # The tensor product rule is used as an ordinary rule
    assert all(not rule.has_tensor_factors for rule, _ in integrands(form, "float64", **options))
    args = (form, num_coefficients, cell_type, "float64")
    assert_matches_fused_loop(args, compile_args, tmp_path, options)


@pytest.mark.parametrize("P", [7, 8])
def test_sum_factorised_contraction_high_degree_elasticity(P):
    """Hex elasticity keeps directional contractions beyond the former buffer limit."""
    form, _ = sum_factorisation_forms(basix.CellType.hexahedron, P)["elasticity"]
    assert sum_factorised(form, "float64")


@pytest.mark.parametrize("name", ["two_quadrature_degrees", "piecewise_coefficient", "linear"])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_sum_factorisation_quadrature_rules(name, cell_type, compile_args, tmp_path):
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
    assert_matches_fused_loop(args, compile_args, tmp_path)
