# Copyright (C) 2024 Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later

from __future__ import annotations

import sys
from collections.abc import Sequence

import basix.ufl
import numpy as np
import pytest
import ufl
from utils import reference_vertex_permutations

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
    entity_local_index: int | Sequence[int] | None = None,
    quadrature_permutation: int | Sequence[int] = 0,
) -> np.ndarray:
    """Call a compiled scalar-valued tabulate_tensor kernel and return its output.

    Pass a pair for `entity_local_index`/`quadrature_permutation` to call
    an interior-facet kernel, which sees both sides of the facet.
    """
    ffi = module.ffi
    xdtype = dtype_to_scalar_dtype(dtype)
    c_type, c_xtype = dtype_to_c_type(dtype), dtype_to_c_type(xdtype)

    A = np.zeros(1, dtype=dtype)
    c = np.array([], dtype=dtype)

    if entity_local_index is None:
        entity_ptr = ffi.NULL
        perm_ptr = ffi.NULL
    else:
        local_index = np.atleast_1d(np.asarray(entity_local_index, dtype=np.int32))
        perm = np.atleast_1d(np.asarray(quadrature_permutation, dtype=np.uint8))
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
    """Grad of a coefficient on a codim-1 (facet) mixed-dimensional submesh.

    Integrating it via a facet measure on the parent mesh must reproduce
    the same value, for every facet of the parent cell, as integrating
    directly on the submesh's own mesh. This is a regression test for a
    bug where the submesh's own dofs were read from the wrong offset
    into the parent's `coordinate_dofs` buffer.
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
def test_multiple_mesh_codim1_facet_rt(dtype, compile_args):
    """A Raviart-Thomas coefficient on a codim-1 (facet) mixed-dimensional submesh.

    An RT coefficient's contravariant Piola pullback needs the submesh's own Jacobian
    directly (not just Grad(SpatialCoordinate)), and `Grad(w)` of that
    Piola-mapped field additionally needs its JacobianInverse -- so this
    exercises the coordinate-dofs gather through two different UFL code
    paths than the plain-Lagrange Grad tests above use. Integrating via
    a facet measure on the tetrahedron parent must reproduce the same
    value, for every facet, as integrating directly on the submesh's
    own triangle mesh.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
    element = basix.ufl.element("RT", "triangle", 1)
    space = ufl.FunctionSpace(codomain, element)
    w_coeff = ufl.Coefficient(space)
    expr = ufl.inner(w_coeff, w_coeff) + ufl.inner(ufl.grad(w_coeff), ufl.grad(w_coeff))

    parent_integral, parent_module = _compile_scalar_form(
        expr * ufl.Measure("ds", domain=domain), dtype, compile_args
    )
    manifold_integral, manifold_module = _compile_scalar_form(
        expr * ufl.Measure("dx", domain=codomain), dtype, compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    # A scalene, asymmetric tetrahedron, so a facet's shape genuinely
    # differs from facet to facet.
    coords = np.array([0.1, 0.2, 0.0, 2.3, 0.1, 0.4, 0.2, 3.1, 0.3, 0.5, 0.6, 4.7], dtype=xdtype)
    w = np.array([1.3, -2.1, 0.7], dtype=dtype)

    facet_vertices = basix.topology(basix.CellType.tetrahedron)[2]
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
    """Grad of a coefficient on a codim-2 (ridge) mixed-dimensional submesh.

    Integrating it via a ridge measure on the parent mesh must reproduce
    the same value, for every edge of the parent cell, as integrating
    directly on the submesh's own mesh.
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


def test_multiple_mesh_codim1_facet_gradient_quadrature_permutation(compile_args):
    """Directly exercise the `quadrature_permutation` input for a codim-1
    facet mixed-dimensional-submesh gradient, independent of any real
    (MPI-partitioned) mesh.

    The submesh's own coordinate dofs are gathered out of the parent's
    `coordinate_dofs` buffer via a permutation-aware closure-dofs table
    (see `ffcx/codegeneration/geometry.py`): flipping
    `quadrature_permutation` for a fixed facet/coordinates must change
    the result, and must match swapping which physical vertex the
    submesh's dof 0 vs. dof 1 corresponds to.
    """
    dtype = "float64"
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
    coords = np.array([0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0, 0.0, 4.0, 4.0, 0.0], dtype=xdtype)
    w = np.array([3.2, 4.0], dtype=dtype)
    facet_local_index = 3  # connects parent-local vertices 2, 3 (basix quadrilateral facet 3)

    A_perm0 = _call_scalar_kernel(
        parent_integral,
        parent_module,
        dtype,
        coords,
        w,
        entity_local_index=facet_local_index,
        quadrature_permutation=0,
    )
    A_perm1 = _call_scalar_kernel(
        parent_integral,
        parent_module,
        dtype,
        coords,
        w,
        entity_local_index=facet_local_index,
        quadrature_permutation=1,
    )

    # Ground truth for each of the two possible submesh-dof <-> vertex
    # correspondences, computed independently by direct manifold assembly.
    verts = coords.reshape(-1, 3)
    coords_forward = np.concatenate([verts[2], verts[3]]).astype(xdtype)
    coords_reversed = np.concatenate([verts[3], verts[2]]).astype(xdtype)
    A_manifold_forward = _call_scalar_kernel(
        manifold_integral, manifold_module, dtype, coords_forward, w
    )
    A_manifold_reversed = _call_scalar_kernel(
        manifold_integral, manifold_module, dtype, coords_reversed, w
    )

    # This facet's tangent is not axis-aligned, so the two permutations
    # must genuinely give different results (proving the kernel is
    # permutation-sensitive), each matching its own vertex ordering
    # (proving it's sensitive in the *correct* way).
    assert not np.isclose(A_perm0[0], A_perm1[0])
    np.testing.assert_allclose(A_perm0, A_manifold_forward, atol=1e-10)
    np.testing.assert_allclose(A_perm1, A_manifold_reversed, atol=1e-10)


def test_multiple_mesh_codim2_ridge_gradient_quadrature_permutation(compile_args):
    """Directly exercise the `quadrature_permutation` input for a codim-2
    ridge mixed-dimensional-submesh gradient (see
    test_multiple_mesh_codim1_facet_gradient_quadrature_permutation).
    """
    dtype = "float64"
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
    coords = np.array([0.0, 0.0, 0.0, 2.0, 0.3, 0.1, 0.2, 3.0, 0.4, 0.5, 0.6, 4.0], dtype=xdtype)
    w = np.array([3.2, 4.0], dtype=dtype)
    ridge_local_index = 0  # connects parent-local vertices 2, 3 (basix tetrahedron ridge 0)

    A_perm0 = _call_scalar_kernel(
        parent_integral,
        parent_module,
        dtype,
        coords,
        w,
        entity_local_index=ridge_local_index,
        quadrature_permutation=0,
    )
    A_perm1 = _call_scalar_kernel(
        parent_integral,
        parent_module,
        dtype,
        coords,
        w,
        entity_local_index=ridge_local_index,
        quadrature_permutation=1,
    )

    verts = coords.reshape(-1, 3)
    coords_forward = np.concatenate([verts[2], verts[3]]).astype(xdtype)
    coords_reversed = np.concatenate([verts[3], verts[2]]).astype(xdtype)
    A_manifold_forward = _call_scalar_kernel(
        manifold_integral, manifold_module, dtype, coords_forward, w
    )
    A_manifold_reversed = _call_scalar_kernel(
        manifold_integral, manifold_module, dtype, coords_reversed, w
    )

    assert not np.isclose(A_perm0[0], A_perm1[0])
    np.testing.assert_allclose(A_perm0, A_manifold_forward, atol=1e-10)
    np.testing.assert_allclose(A_perm1, A_manifold_reversed, atol=1e-10)


def test_multiple_mesh_codim1_2d_facet_gradient_quadrature_permutation(compile_args):
    """Directly exercise all 6 `quadrature_permutation` values for a
    triangle (2D) facet mixed-dimensional-submesh gradient.

    Regression test: the closure-dofs gather (see
    `ffcx/codegeneration/geometry.py`) originally reordered a triangle
    facet's dofs using basix's own `permute_subentity_closure`
    numbering, which disagrees with the `quadrature_permutation`
    ordering DOLFINx/FFCx actually use (`ffcx/ir/elementtables.py`'s
    rotation-major, reflection-minor enumeration) for permutation
    values other than 0 and 1 -- caught by a real 3D DOLFINx assembly
    test, not by this file's original (interval-facet-only) permutation
    test.
    """
    dtype = "float64"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
    element = basix.ufl.element("Lagrange", "triangle", 1)
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
    # An irregular tetrahedron, so a facet's tangent genuinely differs
    # with orientation.
    coords = np.array([0.0, 0.0, 0.0, 2.0, 0.3, 0.1, 0.2, 3.0, 0.4, 0.5, 0.6, 4.0], dtype=xdtype)
    w = np.array([3.2, 4.0, 1.5], dtype=dtype)
    facet_local_index = 0  # connects parent-local vertices 1, 2, 3 (basix tetrahedron facet 0)
    facet_vertices = [1, 2, 3]

    vertex_permutations = reference_vertex_permutations(basix.CellType.triangle)
    results = []
    for perm, vertex_permutation in enumerate(vertex_permutations):
        A_perm = _call_scalar_kernel(
            parent_integral,
            parent_module,
            dtype,
            coords,
            w,
            entity_local_index=facet_local_index,
            quadrature_permutation=perm,
        )
        permuted_vertices = [facet_vertices[v] for v in vertex_permutation]
        manifold_coords = coords.reshape(-1, 3)[permuted_vertices].copy()
        A_manifold = _call_scalar_kernel(
            manifold_integral, manifold_module, dtype, manifold_coords, w
        )
        np.testing.assert_allclose(A_perm, A_manifold, atol=1e-10)
        results.append(A_perm[0])

    # A generic (non-degenerate) facet must give 6 genuinely different
    # results across its 6 possible orientations.
    assert len(set(np.round(results, 8))) == 6


def _point_submesh_forms(dtype: str, compile_args: list[str], measure: str):
    """Compile a codim-2 point-submesh functional and its manifold reference.

    `measure` is the parent measure to integrate over: "dP" (a vertex
    integral) or "ridge" -- on a 2D cell a ridge *is* a vertex, so both
    gather the submesh's coordinate dofs over the parent's vertices.
    """
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(2,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "point", 0, shape=(2,), discontinuous=True))
    # The function space's element must be continuous: `ffcx/analysis.py`
    # rejects discontinuous elements in vertex integrals. The *coordinate*
    # element above stays discontinuous, matching what DOLFINx builds.
    element = basix.ufl.element("Lagrange", "point", 0)
    u = ufl.Coefficient(ufl.FunctionSpace(codomain, element))
    x = ufl.SpatialCoordinate(codomain)
    expr = u * (x[0] + 2 * x[1])

    return (
        _compile_scalar_form(expr * ufl.Measure(measure, domain=domain), dtype, compile_args),
        _compile_scalar_form(expr * ufl.Measure("dx", domain=codomain), dtype, compile_args),
    )


# An irregular triangle, so the three vertices give genuinely different values.
_POINT_SUBMESH_COORDS = [0.3, 0.2, 0.0, 2.0, 0.7, 0.0, 0.6, 3.0, 0.0]


@pytest.mark.parametrize("measure", ["dP", "ridge"])
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
def test_multiple_mesh_codim2_point_spatial_coordinate(dtype, compile_args, measure):
    """Spatial coordinates of a codim-2 (point) mixed-dimensional submesh.

    Integrating them over the parent's vertices must reproduce, for every
    vertex, the value of integrating directly on the submesh's own mesh.
    Covers the vertex-closure-dofs gather for both spellings of a 2D
    codim-2 measure.
    """
    (parent_integral, parent_module), (manifold_integral, manifold_module) = _point_submesh_forms(
        dtype, compile_args, measure
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    coords = np.array(_POINT_SUBMESH_COORDS, dtype=xdtype)
    w = np.array([3.2], dtype=dtype)

    results = []
    for i, verts in enumerate(basix.topology(basix.CellType.triangle)[0]):
        A_parent = _call_scalar_kernel(
            parent_integral, parent_module, dtype, coords, w, entity_local_index=i
        )
        vertex_coords = coords.reshape(-1, 3)[list(verts)].copy()
        A_manifold = _call_scalar_kernel(
            manifold_integral, manifold_module, dtype, vertex_coords, w
        )
        np.testing.assert_allclose(A_parent, A_manifold, atol=1e-10)
        results.append(A_parent[0])

    # A generic triangle must give 3 genuinely different results, i.e. the
    # gather really does follow `entity_local_index`.
    assert len(set(np.round(results, 8))) == 3


@pytest.mark.parametrize("measure", ["dP", "ridge"])
def test_multiple_mesh_codim2_point_quadrature_permutation_invariance(compile_args, measure):
    """A vertex has no orientation, so its gather must ignore the permutation.

    The inverse of the facet/ridge permutation sweeps above, which assert
    the kernel *is* permutation-sensitive. Here a single-row closure table
    is indexed with a literal 0, so a varying (and, for a one-row table,
    out-of-range) `quadrature_permutation` must change nothing.
    """
    dtype = "float64"
    (parent_integral, parent_module), _ = _point_submesh_forms(dtype, compile_args, measure)

    coords = np.array(_POINT_SUBMESH_COORDS, dtype=dtype_to_scalar_dtype(dtype))
    w = np.array([3.2], dtype=dtype)

    for i in range(len(basix.topology(basix.CellType.triangle)[0])):
        values = [
            _call_scalar_kernel(
                parent_integral,
                parent_module,
                dtype,
                coords,
                w,
                entity_local_index=i,
                quadrature_permutation=perm,
            )[0]
            for perm in (0, 1, 5)
        ]
        assert len(set(values)) == 1


def test_geometry_only_submesh_needs_facet_permutations(compile_args):
    """A submesh that enters only through geometry is still mixed-dimensional.

    With no coefficient or argument of its own, the submesh was invisible
    to the `is_mixed_dim` scan, so the integral under-reported
    `needs_facet_permutations` and DOLFINx would skip computing the
    permutation the closure-dofs gather indexes with.

    Note the *values* below are the same for every orientation, and were
    so before the fix too: an integral of a function of the physical
    coordinate over a fixed facet cannot depend on how that facet is
    parameterised, so a geometry-only integrand is orientation-invariant
    by construction. The sweep is a correctness guard on the gather; the
    flag is the regression.
    """
    dtype = "float64"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "tetrahedron", 1, shape=(3,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 1, shape=(3,)))
    x = ufl.SpatialCoordinate(codomain)
    expr = x[0] + 2 * x[1] + 3 * x[2]

    parent_integral, parent_module = _compile_scalar_form(
        expr * ufl.Measure("ds", domain=domain), dtype, compile_args
    )
    manifold_integral, manifold_module = _compile_scalar_form(
        expr * ufl.Measure("dx", domain=codomain), dtype, compile_args
    )

    assert parent_integral.needs_facet_permutations

    # An irregular tetrahedron, so each facet genuinely differs.
    coords = np.array(
        [0.0, 0.0, 0.0, 2.0, 0.3, 0.1, 0.2, 3.0, 0.4, 0.5, 0.6, 4.0],
        dtype=dtype_to_scalar_dtype(dtype),
    )
    w = np.array([], dtype=dtype)
    facet_local_index = 0
    facet_vertices = basix.topology(basix.CellType.tetrahedron)[2][facet_local_index]

    vertex_permutations = reference_vertex_permutations(basix.CellType.triangle)
    for perm, vertex_permutation in enumerate(vertex_permutations):
        A_perm = _call_scalar_kernel(
            parent_integral,
            parent_module,
            dtype,
            coords,
            w,
            entity_local_index=facet_local_index,
            quadrature_permutation=perm,
        )
        permuted_vertices = [facet_vertices[v] for v in vertex_permutation]
        manifold_coords = coords.reshape(-1, 3)[permuted_vertices].copy()
        A_manifold = _call_scalar_kernel(
            manifold_integral, manifold_module, dtype, manifold_coords, w
        )
        np.testing.assert_allclose(A_perm, A_manifold, atol=1e-10)


def test_multiple_mesh_codim1_interior_facet_restriction_offset(compile_args):
    """The "-" side of an interior facet starts one *parent* cell in.

    `coordinate_dofs` holds the two cells of the integration domain back
    to back, so the "-" restriction's offset is the parent's scalar dof
    count -- not the submesh's, which is smaller. A degree-2 triangle
    parent with a degree-2 interval submesh makes the two strides (18 vs
    9) far enough apart that the wrong one reads a mix of both cells.

    Ground truth needs no second mesh: feed the *same* cell on both sides
    of the facet, and the "-" restriction must then agree exactly with
    the "+" restriction.
    """
    dtype = "float64"
    domain = ufl.Mesh(basix.ufl.element("Lagrange", "triangle", 2, shape=(2,)))
    codomain = ufl.Mesh(basix.ufl.element("Lagrange", "interval", 2, shape=(2,)))
    x = ufl.SpatialCoordinate(codomain)
    dS = ufl.Measure("dS", domain=domain)

    minus_integral, minus_module = _compile_scalar_form(
        (x[0]("-") + 2 * x[1]("-")) * dS, dtype, compile_args
    )
    plus_integral, plus_module = _compile_scalar_form(
        (x[0]("+") + 2 * x[1]("+")) * dS, dtype, compile_args
    )

    xdtype = dtype_to_scalar_dtype(dtype)
    rng = np.random.default_rng(7)
    cell = rng.uniform(0.0, 1.0, (6, 3))
    cell[:, 2] = 0.0
    coords = np.concatenate([cell, cell]).ravel().astype(xdtype)
    w = np.array([], dtype=dtype)

    for facet in range(3):
        for perm in (0, 1):
            A_minus = _call_scalar_kernel(
                minus_integral,
                minus_module,
                dtype,
                coords,
                w,
                entity_local_index=(facet, facet),
                quadrature_permutation=(perm, perm),
            )
            A_plus = _call_scalar_kernel(
                plus_integral,
                plus_module,
                dtype,
                coords,
                w,
                entity_local_index=(facet, facet),
                quadrature_permutation=(perm, perm),
            )
            np.testing.assert_allclose(A_minus, A_plus, atol=1e-10)
            assert abs(A_plus[0]) > 1e-3
