# Copyright (C) 2023 Igor A. Baratta
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

import itertools

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype
from ffcx.element_interface import basix_index
from ffcx.ir.elementtables import (
    TensorProductTablePlaceholder,
    analyse_table_type,
    analyse_tensor_product_table_type,
    clamp_table_small_numbers,
)


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


def distorted_cell(cell_type):
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
    return np.array(coords, dtype=np.float64)


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
    form = {
        "two_quadrature_degrees": u * v * ufl.dx(degree=2)
        + ufl.inner(ufl.grad(u), ufl.grad(v)) * ufl.dx(degree=5),
        "piecewise_coefficient": f * u * v * one_point,
        "linear": f * v * one_point,
    }[name]

    def tabulate(options, cache_dir):
        compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
            [form], options=options, cache_dir=cache_dir, cffi_extra_compile_args=compile_args
        )
        form0 = compiled_forms[0]
        integral = form0.form_integrals[form0.form_integral_offsets[module.lib.cell]]
        ndofs = V.ufl_element().dim
        A = np.zeros(ndofs ** len(form.arguments()))
        w = np.linspace(0.3, 1.7, ndofs)
        c = np.zeros(1)
        coords = distorted_cell(cell_type)
        ffi = module.ffi
        integral.tabulate_tensor_float64(
            ffi.cast("double *", A.ctypes.data),
            ffi.cast("double *", w.ctypes.data),
            ffi.cast("double *", c.ctypes.data),
            ffi.cast("double *", coords.ctypes.data),
            ffi.NULL,
            ffi.NULL,
            ffi.NULL,
        )
        return A

    A_ref = tabulate({}, tmp_path / "ref")
    A = tabulate({"sum_factorization": True}, tmp_path / "sf")
    np.testing.assert_allclose(A, A_ref, rtol=1e-12, atol=1e-12 * np.max(np.abs(A_ref)))


@pytest.mark.parametrize("P", [1, 2, 3])
@pytest.mark.parametrize("cell_type", [basix.CellType.quadrilateral, basix.CellType.hexahedron])
def test_tensor_product_table_type(cell_type, P):
    """The table type derived from the 1D factors matches the one of the full table."""
    d = cell_to_gdim(cell_type)
    element = basix.create_tp_element(
        basix.ElementFamily.P, cell_type, P, basix.LagrangeVariant.gll_warped
    )
    factor = element.get_tensor_product_representation()[0][0]
    rules = [basix.make_quadrature(basix.CellType.interval, q)[0] for q in (0, 1, 2 * P)]
    # Collocated with the dofs, giving quadrature tables
    rules.append(basix.create_lattice(basix.CellType.interval, P, basix.LatticeType.gll, True))
    for x in rules:
        points = np.array(list(itertools.product(*[x[:, 0]] * d)))
        for derivatives in itertools.product(range(3), repeat=d):
            if sum(derivatives) > 2:
                continue
            table = element.tabulate(sum(derivatives), points)[basix_index(derivatives)][..., 0]
            table = clamp_table_small_numbers(table.reshape(1, 1, *table.shape))
            factors = [
                clamp_table_small_numbers(
                    factor.tabulate(n, x)[n][..., 0].reshape(1, 1, x.shape[0], -1)
                )
                for n in derivatives
            ]
            assert analyse_tensor_product_table_type(factors) == analyse_table_type(table)


def test_tensor_product_table_type_clamped_to_zero():
    """Factors whose products all clamp to zero give a zero table."""
    rng = np.random.default_rng(0)
    factors = [1e-5 * (1 + rng.random((1, 1, 3, 2))) for _ in range(2)]
    full = np.einsum("pi,qj->pqij", factors[0][0, 0], factors[1][0, 0]).reshape(1, 1, 9, 4)
    assert analyse_table_type(clamp_table_small_numbers(full)) == "zeros"
    assert analyse_tensor_product_table_type(factors) == "zeros"


def test_tensor_product_table_placeholder():
    """The placeholder of a table that isn't formed has a shape but no values."""
    table = TensorProductTablePlaceholder((1, 1, 9, 4))
    assert table.shape == (1, 1, 9, 4)
    with pytest.raises(RuntimeError):
        np.asarray(table)
    with pytest.raises(RuntimeError):
        table[0]
