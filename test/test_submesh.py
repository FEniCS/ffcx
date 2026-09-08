# Copyright (C) 2024 Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from __future__ import annotations

import sys

import basix.ufl
import numpy as np
import pytest
import ufl

import ffcx.codegeneration.jit
from ffcx.codegeneration.utils import dtype_to_c_type, dtype_to_scalar_dtype


def compute_tensor(forms: list[ufl.form.Form], dtype: str, compile_args: list[str]):
    """Helper-function to compute matrix for a P1-Lagrange problem"""
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        forms, options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )

    ffi = module.ffi
    form0 = compiled_forms[0]
    offsets = form0.form_integral_offsets
    cell = module.lib.cell
    assert offsets[cell + 1] - offsets[cell] == 1
    integral_id = form0.form_integral_ids[offsets[cell]]
    assert integral_id == -1

    default_integral = form0.form_integrals[offsets[cell]]

    A = np.zeros((3, 3), dtype=dtype)
    w = np.array([], dtype=dtype)
    c = np.array([], dtype=dtype)

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array([[1.0, 2.0, 0.0], [1.5, 2.3, 0.0], [6.0, 1.8, 0.0]], dtype=xdtype)

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
    return A


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
def test_multiple_mesh_codim0(dtype, compile_args):
    # Define coordinate element and element used in parent and sub-mesh
    element = basix.ufl.element("Lagrange", "triangle", 1)
    coordinate_element = basix.ufl.element("Lagrange", "triangle", 1, shape=(2,))

    domain = ufl.Mesh(coordinate_element)
    space = ufl.FunctionSpace(domain, element)
    u_parent = ufl.TrialFunction(space)

    # Create submesh and functionspace on submesh
    sub_domain = ufl.Mesh(coordinate_element)
    subspace = ufl.FunctionSpace(sub_domain, element)
    v_sub = ufl.TestFunction(subspace)

    #
    a = ufl.inner(u_parent.dx(0), v_sub.dx(0)) * ufl.dx(domain=domain)

    A = compute_tensor([a], dtype, compile_args)

    # Compute reference solution on with test and trial function from same mesh
    v_parent = ufl.TestFunction(space)
    a_org = ufl.inner(u_parent.dx(0), v_parent.dx(0)) * ufl.dx(domain=domain)
    A_org = compute_tensor([a_org], dtype, compile_args)

    np.testing.assert_allclose(A, A_org)


def _compile_scalar_form(form, dtype: str, compile_args: list[str]):
    """Compile a single-integral scalar functional and return its kernel."""
    compiled_forms, module, _ = ffcx.codegeneration.jit.compile_forms(
        [form], options={"scalar_type": dtype}, cffi_extra_compile_args=compile_args
    )
    return compiled_forms[0].form_integrals[0], module


def _call_scalar_kernel(
    integral,
    module,
    dtype: str,
    coordinate_dofs: np.ndarray,
    coefficients: np.ndarray,
    entity_local_index: int | None = None,
    quadrature_permutation: int = 0,
) -> np.ndarray:
    """Call a compiled scalar-valued tabulate_tensor kernel and return its output."""
    ffi = module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)

    A = np.zeros(1, dtype=dtype)
    c = np.array([], dtype=dtype)

    if entity_local_index is None:
        entity_ptr = ffi.NULL
        perm_ptr = ffi.NULL
    else:
        local_index = np.array([entity_local_index], dtype=np.int32)
        perm = np.array([quadrature_permutation], dtype=np.uint8)
        entity_ptr = ffi.cast("int *", local_index.ctypes.data)
        perm_ptr = ffi.cast("uint8_t *", perm.ctypes.data)

    kernel = getattr(integral, f"tabulate_tensor_{dtype}")
    kernel(
        ffi.cast(f"{c_type} *", A.ctypes.data),
        ffi.cast(f"{c_type} *", coefficients.ctypes.data),
        ffi.cast(f"{c_type} *", c.ctypes.data),
        ffi.cast(f"{c_xtype} *", coordinate_dofs.ctypes.data),
        entity_ptr,
        perm_ptr,
        ffi.NULL,
    )
    return A


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
def test_multiple_mesh_codim1_facet_gradient(dtype, compile_args):
    """Grad of a coefficient on a codim-1 (facet) companion domain.

    Integrating it via a facet measure on the parent mesh must reproduce
    the same value, for every facet of the parent cell, as integrating
    directly on the companion domain's own mesh. This is a regression
    test for a bug where the companion domain's own (unavailable)
    Jacobian was used instead of the parent's FacetJacobian.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "quadrilateral", 1, shape=(2,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(2,)))
    element = basix.ufl.element("Lagrange", "interval", 1)
    space = ufl.FunctionSpace(codomain, element)
    u = ufl.Coefficient(space)
    expr = ufl.Dx(u, 0)

    parent_integral, parent_module = _compile_scalar_form(
        expr * ufl.Measure("ds", domain=domain), dtype, compile_args
    )
    manifold_integral, manifold_module = _compile_scalar_form(
        expr * ufl.Measure("dx", domain=codomain), dtype, compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    # An irregular (non-parallelogram) quadrilateral, so a facet's
    # tangent genuinely differs from facet to facet.
    coords = np.array([0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 4.0, 4.0, 0.0], dtype=xdtype)
    w = np.array([3.2, 4.0], dtype=dtype)

    facet_vertices = basix.topology(basix.CellType.quadrilateral)[1]
    for i, verts in enumerate(facet_vertices):
        A_parent = _call_scalar_kernel(
            parent_integral, parent_module, dtype, coords, w, entity_local_index=i
        )
        facet_coords = coords.reshape(-1, 3)[list(verts)].copy()
        A_manifold = _call_scalar_kernel(manifold_integral, manifold_module, dtype, facet_coords, w)
        np.testing.assert_allclose(A_parent, A_manifold, atol=1e-10)


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
def test_multiple_mesh_codim2_ridge_gradient(dtype, compile_args):
    """Grad of a coefficient on a codim-2 (ridge) companion domain.

    Integrating it via a ridge measure on the parent mesh must reproduce
    the same value, for every edge of the parent cell, as integrating
    directly on the companion domain's own mesh.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "interval", 1, shape=(3,)))
    element = basix.ufl.element("Lagrange", "interval", 1)
    space = ufl.FunctionSpace(codomain, element)
    u = ufl.Coefficient(space)
    expr = ufl.Dx(u, 0)

    parent_integral, parent_module = _compile_scalar_form(
        expr * ufl.Measure("ridge", domain=domain), dtype, compile_args
    )
    manifold_integral, manifold_module = _compile_scalar_form(
        expr * ufl.Measure("dx", domain=codomain), dtype, compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    # An irregular tetrahedron, so an edge's tangent genuinely differs
    # from edge to edge.
    coords = np.array([0.0, 0.0, 0.0, 2.0, 0.3, 0.1, 0.2, 3.0, 0.4, 0.5, 0.6, 4.0], dtype=xdtype)
    w = np.array([3.2, 4.0], dtype=dtype)

    ridge_vertices = basix.topology(basix.CellType.tetrahedron)[1]
    for i, verts in enumerate(ridge_vertices):
        A_parent = _call_scalar_kernel(
            parent_integral, parent_module, dtype, coords, w, entity_local_index=i
        )
        ridge_coords = coords.reshape(-1, 3)[list(verts)].copy()
        A_manifold = _call_scalar_kernel(manifold_integral, manifold_module, dtype, ridge_coords, w)
        np.testing.assert_allclose(A_parent, A_manifold, atol=1e-10)
