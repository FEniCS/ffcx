# Copyright (C) 2026 Jørgen S. Dokken
#
# This file is part of FFCx. (https://www.fenicsproject.org)
#
# SPDX-License-Identifier:    LGPL-3.0-or-later
"""Unit tests for ffcx.codegeneration.geometry's static reference tables."""

import basix
import basix.ufl
import pytest
from conftest import _reference_vertex_permutations

from ffcx.codegeneration import geometry


def _coordinate_element(cellname, gdim, degree=1):
    # equispaced: basix's default (GLL) variant doesn't support degree > 2
    # on pyramids; equispaced avoids that unrelated limitation uniformly.
    return basix.ufl.element(
        "Lagrange",
        cellname,
        degree,
        shape=(gdim,),
        lagrange_variant=basix.LagrangeVariant.equispaced,
    )


@pytest.mark.parametrize(
    "cellname,gdim,nperm,nfacets,ndofs",
    [
        ("triangle", 2, 2, 3, 2),
        ("quadrilateral", 2, 2, 4, 2),
        ("tetrahedron", 3, 6, 4, 3),
        ("hexahedron", 3, 8, 6, 4),
    ],
)
def test_facet_closure_dofs_shape(cellname, gdim, nperm, nfacets, ndofs):
    """The table has one row per possible facet orientation, one column
    per facet, and one entry per scalar dof on that facet.
    """
    element = _coordinate_element(cellname, gdim)
    decl = geometry.write_table("facet_closure_dofs", cellname, element)
    assert decl.symbol.name == f"{cellname}_facet_closure_dofs"
    assert decl.values.shape == (nperm, nfacets, ndofs)


@pytest.mark.parametrize("cellname,gdim", [("triangle", 2), ("quadrilateral", 2)])
def test_facet_closure_dofs_p1_forward_and_reversed(cellname, gdim):
    """For a P1 (vertex-only) coordinate element, permutation 0 must be
    the identity and permutation 1 the reverse -- this is the exact
    orientation-correction the mixed-domain gradient bug needed.
    """
    element = _coordinate_element(cellname, gdim)
    table = geometry.write_table("facet_closure_dofs", cellname, element).values
    for facet in range(table.shape[1]):
        forward = table[0, facet]
        reversed_ = table[1, facet]
        assert list(reversed_) == list(reversed(forward))


@pytest.mark.parametrize(
    "cellname,gdim,facet_celltype",
    [
        ("triangle", 2, basix.CellType.interval),
        ("quadrilateral", 2, basix.CellType.interval),
        ("tetrahedron", 3, basix.CellType.triangle),
        ("hexahedron", 3, basix.CellType.quadrilateral),
    ],
)
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_facet_closure_dofs_vertex_part_matches_reference_permutation(
    cellname, gdim, facet_celltype, degree
):
    """The vertex-dof sub-part of every row must match the independently
    validated (via `permute_quadrature_interval`/`triangle`/`quadrilateral`)
    reference-vertex permutation, regardless of the coordinate element's
    degree -- this is what distinguishes basix's general
    `permute_subentity_closure_inv`-based reordering (which also reorders
    edge-/face-interior dofs) from a naive vertex-only one, which cannot
    handle higher-degree coordinate elements.
    """
    element = _coordinate_element(cellname, gdim, degree=degree)
    be = element.basix_element
    entity_dim = gdim - 1
    expected_vertex_perms = _reference_vertex_permutations(facet_celltype)
    num_vertices = len(basix.geometry(facet_celltype))
    table = geometry.write_table("facet_closure_dofs", cellname, element).values
    for facet in range(table.shape[1]):
        base = be.entity_closure_dofs[entity_dim][facet]
        for perm, vertex_permutation in enumerate(expected_vertex_perms):
            expected_vertices = [base[v] for v in vertex_permutation]
            assert list(table[perm, facet, :num_vertices]) == expected_vertices


@pytest.mark.parametrize(
    "cellname,gdim,nridges",
    [
        ("tetrahedron", 3, 6),
        ("hexahedron", 3, 12),
        ("prism", 3, 9),
        ("pyramid", 3, 8),
    ],
)
def test_ridge_closure_dofs_shape(cellname, gdim, nridges):
    """Same shape contract as facet_closure_dofs, one dimension down.

    Ridges (codim-2 entities) of any 3D cell are edges -- i.e. intervals
    -- regardless of the cell's own shape, so every standard 3D cell type
    has the same two possible orientations (forward/reversed).
    """
    element = _coordinate_element(cellname, gdim)
    decl = geometry.write_table("ridge_closure_dofs", cellname, element)
    assert decl.symbol.name == f"{cellname}_ridge_closure_dofs"
    assert decl.values.shape == (2, nridges, 2)


@pytest.mark.parametrize("cellname", ["tetrahedron", "hexahedron", "prism", "pyramid"])
def test_ridge_closure_dofs_p1_forward_and_reversed(cellname):
    """Same forward/reversed contract as the facet case, for ridges."""
    element = _coordinate_element(cellname, 3)
    table = geometry.write_table("ridge_closure_dofs", cellname, element).values
    for ridge in range(table.shape[1]):
        forward = table[0, ridge]
        reversed_ = table[1, ridge]
        assert list(reversed_) == list(reversed(forward))


@pytest.mark.parametrize("cellname", ["tetrahedron", "hexahedron", "prism", "pyramid"])
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_ridge_closure_dofs_vertex_part_matches_reference_permutation(cellname, degree):
    """Same higher-degree cross-check as the facet case, for ridges."""
    element = _coordinate_element(cellname, 3, degree=degree)
    be = element.basix_element
    interval = basix.CellType.interval
    expected_vertex_perms = _reference_vertex_permutations(interval)
    table = geometry.write_table("ridge_closure_dofs", cellname, element).values
    for ridge in range(table.shape[1]):
        base = be.entity_closure_dofs[1][ridge]
        for perm, vertex_permutation in enumerate(expected_vertex_perms):
            expected_vertices = [base[v] for v in vertex_permutation]
            assert list(table[perm, ridge, :2]) == expected_vertices


@pytest.mark.parametrize(
    "cellname,gdim,entity_dim",
    [
        ("triangle", 2, 1),
        ("quadrilateral", 2, 1),
        ("tetrahedron", 3, 2),
        ("hexahedron", 3, 2),
    ],
)
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_facet_closure_dofs_every_permutation_is_a_reordering(cellname, gdim, entity_dim, degree):
    """Regardless of cell type, orientation, or coordinate element degree,
    every permutation row must contain exactly the same set of dof
    indices as the unpermuted (basix reference) closure -- only their
    order may change.
    """
    element = _coordinate_element(cellname, gdim, degree=degree)
    be = element.basix_element
    table = geometry.write_table("facet_closure_dofs", cellname, element).values
    for facet in range(table.shape[1]):
        expected = set(be.entity_closure_dofs[entity_dim][facet])
        for perm in range(table.shape[0]):
            assert set(table[perm, facet]) == expected


@pytest.mark.parametrize("cellname", ["tetrahedron", "hexahedron", "prism", "pyramid"])
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_ridge_closure_dofs_every_permutation_is_a_reordering(cellname, degree):
    """Same reordering-only invariant as the facet case, for ridges."""
    element = _coordinate_element(cellname, 3, degree=degree)
    be = element.basix_element
    table = geometry.write_table("ridge_closure_dofs", cellname, element).values
    for ridge in range(table.shape[1]):
        expected = set(be.entity_closure_dofs[1][ridge])
        for perm in range(table.shape[0]):
            assert set(table[perm, ridge]) == expected


@pytest.mark.parametrize("cellname", ["triangle", "quadrilateral", "tetrahedron", "hexahedron"])
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_peak_closure_dofs_shape_and_value(cellname, degree):
    """A vertex ("peak", the codim-3 entity) closure is always its own
    single dof, regardless of cell type or coordinate element degree --
    there is only ever one possible orientation (nothing to permute).
    """
    gdim = 2 if cellname in ("triangle", "quadrilateral") else 3
    element = _coordinate_element(cellname, gdim, degree=degree)
    be = element.basix_element
    celltype = getattr(basix.CellType, cellname)
    num_vertices = len(basix.topology(celltype)[0])
    decl = geometry.write_table("peak_closure_dofs", cellname, element)
    assert decl.symbol.name == f"{cellname}_peak_closure_dofs"
    assert decl.values.shape == (1, num_vertices, 1)
    for v in range(num_vertices):
        assert list(decl.values[0, v]) == be.entity_closure_dofs[0][v]


@pytest.mark.xfail(
    raises=NotImplementedError,
    reason="Prism/pyramid facets are not uniform (mixed triangle/quadrilateral), matching the "
    "same restriction `ffcx/ir/elementtables.py`'s mixed-dimensional-submesh element "
    "permutation already has.",
)
def test_facet_closure_dofs_unsupported_cell_type():
    element = _coordinate_element("prism", 3)
    geometry.write_table("facet_closure_dofs", "prism", element)


@pytest.mark.parametrize("cellname", ["triangle", "quadrilateral"])
@pytest.mark.parametrize("degree", [1, 2, 3])
def test_ridge_closure_dofs_2d_is_the_vertex_closure(cellname, degree):
    """A 2D cell's ridges are its vertices, so its ridge table is the vertex table.

    It therefore has a single orientation row (a vertex has nothing to
    permute), unlike the 3D ridge tables above whose ridges are edges.
    """
    gdim = 2
    element = _coordinate_element(cellname, gdim, degree=degree)
    be = element.basix_element
    celltype = getattr(basix.CellType, cellname)
    num_vertices = len(basix.topology(celltype)[0])

    decl = geometry.write_table("ridge_closure_dofs", cellname, element)
    assert decl.symbol.name == f"{cellname}_ridge_closure_dofs"
    assert decl.values.shape == (1, num_vertices, 1)
    for v in range(num_vertices):
        assert list(decl.values[0, v]) == be.entity_closure_dofs[0][v]

    # The peak table is the same values under the other name.
    peak = geometry.write_table("peak_closure_dofs", cellname, element)
    assert (peak.values == decl.values).all()
